*     cmplxd.f (FORTRAN 77)
*     Demonstration of COMPLEX numbers
*
*     Prints the values of e ** (j * i * pi / 4) for i = 0, 1, 2, ..., 7
*         where j is the imaginary number sqrt(-1)

      pure function square(x)
        real, intent(in) :: x
        real :: square
        square = x * x
      end function

      program main
        real :: a, b, square
        a = 2.0
        b = square(a)
		print*,b
        ! After invoking the square(.) pure function, we can be sure that 
        ! besides assigning the output value of square(a) to b,
        ! nothing else has been changed.
      end program main
	  