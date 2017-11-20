*     cmplxd.f (FORTRAN 77)
*     Demonstration of COMPLEX numbers
*
*     Prints the values of e ** (j * i * pi / 4) for i = 0, 1, 2, ..., 7
*         where j is the imaginary number sqrt(-1)

      subroutine square_cube(i,isquare,icube)
        integer, intent(in)  :: i             ! input
        integer, intent(out) :: isquare,icube ! output
        isquare = i**2
        icube   = i**3
      end subroutine square_cube

      program xx
       implicit none
       integer :: i,isq,icub
       i = 4
       call square_cube(i,isq,icub)
       print*,"i,i^2,i^3=",i,isq,icub
      end program xx
	  