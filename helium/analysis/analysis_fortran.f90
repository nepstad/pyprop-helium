subroutine SumAngularProjection(angularData, radialProj, angularFactor, adExtent0, adExtent1)
	implicit none

	integer, parameter                             :: dbl = selected_real_kind(p=13) 

	complex (kind = dbl), intent(inout)            :: angularFactor
	integer, intent(in)                            :: adExtent0, adExtent1

	complex (kind = dbl), dimension(adExtent0, adExtent1), intent(inout)    :: angularData
	complex (kind = dbl), dimension(adExtent0, adExtent1), intent(inout)    :: radialProj

	integer                                        :: i, j

	do j = 1, adExtent1
!		do i = 1, adExtent0
		angularData(:,j) = angularData(:,j) + angularFactor * radialProj(:,j)
!		enddo
	enddo

end subroutine SumAngularProjection