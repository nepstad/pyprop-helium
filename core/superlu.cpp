#include <core/potential/dynamicpotentialevaluator.h>
#include <slu_zdefs.h>

using namespace blitz;

template<int Rank>
class SuperLUSolver
{
public:
	typedef blitz::Array<cplx, 1> VectorType;
	typedef blitz::Array<int, 1> VectorTypeInt;
	typedef blitz::Array<double, 1> VectorTypeReal;

private:
	bool Initialized;
	bool UseNaturalOrdering;
	VectorType MatrixData;
	blitz::Array<cplx, Rank> TempData;
	VectorTypeInt RowIndices;
	VectorTypeInt ColStartIndices;
	VectorTypeInt PermutationRow;
	VectorTypeInt PermutationCol;
	
	SuperMatrix A;
	SuperMatrix L; 
	SuperMatrix U;
	VectorTypeReal R;
	VectorTypeReal C;
	VectorTypeInt ETree;
	superlu_options_t Options;
    SuperLUStat_t Statistics;
	char Equed[2];
	
public:
	SuperLUSolver()
	{
		Initialized = false;
	}

	~SuperLUSolver()
	{
		if (Initialized)
		{
			SUPERLU_FREE(A.Store);
			Destroy_SuperNode_Matrix(&L);
			Destroy_CompCol_Matrix(&U);
		}
	}

	void Setup(int matrixSize, VectorType matrixData, VectorTypeInt rowIndices, VectorTypeInt colStartIndices)
	{
		/* Set the default input options:
		options.Fact = DOFACT;
		options.Equil = YES;
		options.ColPerm = COLAMD;
		options.DiagPivotThresh = 1.0;
		options.Trans = NOTRANS;
		options.IterRefine = NOREFINE;
		options.SymmetricMode = NO;
		options.PivotGrowth = NO;
		options.ConditionNumber = NO;
		options.PrintStat = YES;
		*/
		cout << "Starting SuperLU Setup" << endl;
		set_default_options(&Options);

		//Options.ColPerm = NATURAL;
		Options.ColPerm = MMD_ATA;
		Options.DiagPivotThresh = 0.1;
		Options.SymmetricMode = YES;


		//Copy matrix data so it won't go out of scope
		MatrixData.reference( matrixData.copy() );
		RowIndices.reference( rowIndices.copy() );
		ColStartIndices.reference( colStartIndices.copy() );
		//We need permutations of rows and cols
		PermutationRow.resize( matrixSize*2 );
		PermutationCol.resize( matrixSize*2 );
		R.resize( matrixSize );
		C.resize( matrixSize );
		ETree.resize( matrixSize );

		int N = matrixSize;
		int nnz = MatrixData.size();

		//Create compressed col matrix from matrix data
		doublecomplex* data = (doublecomplex*) MatrixData.data();
    	zCreate_CompCol_Matrix(&A, N, N, nnz, data, RowIndices.data(), ColStartIndices.data(), SLU_NC, SLU_Z, SLU_GE);

   		//Create an empty right hand side to create factorization without solving anything
		int nrhs = 0;
		doublecomplex* rhs = 0;
		SuperMatrix B, X;
    	zCreate_Dense_Matrix(&B, N, nrhs, rhs, N, SLU_DN, SLU_Z, SLU_GE);
    	zCreate_Dense_Matrix(&X, N, nrhs, rhs, N, SLU_DN, SLU_Z, SLU_GE);

		//Factorize A
		int lwork = 0;
		doublecomplex *work = 0;

		double rpg; //reciprocal growth
		double ferr; //err[nrhs]
		double berr; //err[nrhs]
		double rcond;
		mem_usage_t memUsage;

		int info;
    	StatInit(&Statistics);
		Options.Fact = DOFACT;
		zgssvx(&Options, &A, PermutationCol.data(), PermutationRow.data(), ETree.data(), Equed, 
			R.data(), C.data(), &L, &U, work, lwork, &B, &X, &rpg, &rcond, &ferr, &berr,
           &memUsage, &Statistics, &info);

		if (info != 0)
		{
			cout << "Error from SuperLU: " << info << endl;
			throw std::runtime_error("Could not factorize matrix with superlu");
		}

		SCformat* Lstore = (SCformat *) L.Store;
		NCformat* Ustore = (NCformat *) U.Store;
		cout << "SuperLU-factorization of matrix returned " << info << endl;
		cout << "    Nonzero elemts in L = " << Lstore->nnz << endl;
		cout << "    Nonzero elemts in U = " << Ustore->nnz << endl;
		cout << "    Fill factor         = " << (double)(Lstore->nnz + Ustore->nnz - N) / (double) nnz << endl;

    	//Destroy temporary matrix B
		SUPERLU_FREE(B.Store);
		SUPERLU_FREE(X.Store);

		cout << "Finishing SuperLU Setup" << endl;

		Initialized = true;
	}

	void Solve(blitz::Array<cplx, Rank> rhs)
	{
		//allocate Temp-array if we don't already have it
		if (rhs.size() != TempData.size())
		{
			TempData.resize(rhs.shape());
		}

		//Setup right hand side in correct format
   		int nrhs = 1;
		doublecomplex* rhsData = (doublecomplex*)rhs.data();
		doublecomplex* tempData = (doublecomplex*)TempData.data();
		SuperMatrix B, X;
	 	zCreate_Dense_Matrix(&B, rhs.size(), nrhs, rhsData, rhs.size(), SLU_DN, SLU_Z, SLU_GE);
	 	zCreate_Dense_Matrix(&X, rhs.size(), nrhs, tempData, rhs.size(), SLU_DN, SLU_Z, SLU_GE);

		//Solve
		int lwork = 0;
		doublecomplex *work = 0;
		double rpg; //reciprocal growth
		double ferr; //err[nrhs]
		double berr; //err[nrhs]
		double rcond;
		mem_usage_t memUsage;

		int info;
		Options.Fact = FACTORED;
		zgssvx(&Options, &A, PermutationCol.data(), PermutationRow.data(), ETree.data(), Equed, 
			R.data(), C.data(), &L, &U, work, lwork, &B, &X, &rpg, &rcond, &ferr, &berr,
           &memUsage, &Statistics, &info);


		if (info != 0)
		{
			cout << "Error from SuperLU: " << info << endl;
			throw std::runtime_error("Could not factorize matrix with superlu");
		}

		//Copy data from temp data storage to overwrite input data
		rhs = TempData;

		//free temp datastructures
		SUPERLU_FREE(B.Store);
		SUPERLU_FREE(X.Store);
	}

	void PrintStatistics()
	{
		StatPrint(&Statistics);
	}
};

extern "C"
{

/*! \brief

 <pre>
    Purpose   
    =======   

    sp_ienv() is inquired to choose machine-dependent parameters for the
    local environment. See ISPEC for a description of the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    Arguments   
    =========   

    ISPEC   (input) int
            Specifies the parameter to be returned as the value of SP_IENV.   
            = 1: the panel size w; a panel consists of w consecutive
	         columns of matrix A in the process of Gaussian elimination.
		 The best value depends on machine's cache characters.
            = 2: the relaxation parameter relax; if the number of
	         nodes (columns) in a subtree of the elimination tree is less
		 than relax, this subtree is considered as one supernode,
		 regardless of their row structures.
            = 3: the maximum size for a supernode;
	    = 4: the minimum row dimension for 2-D blocking to be used;
	    = 5: the minimum column dimension for 2-D blocking to be used;
	    = 6: the estimated fills factor for L and U, compared with A;
	    
   (SP_IENV) (output) int
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if SP_IENV = -k, the k-th argument had an illegal value. 
  
    ===================================================================== 
</pre>
*/
int sp_ienv(int ispec)
{
    int i;
	
    switch (ispec) {
	case 1: return (10);
	case 2: return (5);
	case 3: return (100);
	case 4: return (200);
	case 5: return (40);
        case 6: return (5);
    }

    /* Invalid value for ISPEC */
    i = 1;
    xerbla_("sp_ienv", &i);
    return 0;

} /* sp_ienv_ */

}



