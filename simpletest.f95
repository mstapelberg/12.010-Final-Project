        Program Monte
            implicit none

            integer, parameter :: seed = 1
            integer :: i,j,k
            real(kind = 8), dimension(10) :: phi, theta, r
            real(kind = 8), dimension(10) :: x, y, z
            real(kind = 8) :: PI = 3.14159
            real(kind = 8) :: ran
            call srand(seed)
            print *, rand(), rand(), rand(), rand()
            print *, rand(seed), rand(), rand(), rand()
            !Here we fill the arrays
            do i= 1,10
                phi(i) = 2*PI*ran(-4)
            end do
            do j = 1,10
                write(*,*) phi(j)
            end do
        End Program Monte

        FUNCTION ran(idum)
            IMPLICIT NONE
            INTEGER, PARAMETER :: K4B=selected_int_kind(9)
            INTEGER(K4B), INTENT(INOUT) :: idum
            REAL :: ran
                !“Minimal” random number generator of Park and Miller
                !combined with a Marsaglia shift
                !sequence. Returns a uniform random deviate between 0.0 and
                !1.0 (exclusive of the endpoint
                !values). This fully portable, scalar generator has the
                !“traditional” (not Fortran 90) calling
                !sequence with a random deviate as the returned function
                !value: call with idum a negative
                !integer to initialize; thereafter, do not alter idum except
                !to reinitialize. The period of this
                !generator is about 3.1 × 10 18 .
            INTEGER(K4B), PARAMETER ::IA=16807,IM=2147483647,IQ=127773
            INTEGER(K4B), PARAMETER ::IR = 2836
            REAL, SAVE :: am
            INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
            if (idum <= 0 .or. iy < 0) then
                !Initialize.
                am=nearest(1.0,-1.0)/IM
                iy=ior(ieor(888889999,abs(idum)),1)
                ix=ieor(777755555,abs(idum))
                idum=abs(idum)+1
                !Set idum positive.
            end if
            ix=ieor(ix,ishft(ix,13))
            !Marsaglia shift sequence with period 2 32 − 1.
            ix=ieor(ix,ishft(ix,-17))
            ix=ieor(ix,ishft(ix,5))
            k=iy/IQ
            !Park-Miller sequence by Schrage’s method,
            !period 2 31 − 2.
            iy=IA*(iy-k*IQ)-IR*k
            if (iy < 0) iy=iy+IM
            ran=am*ior(iand(IM,ieor(ix,iy)),1)
            !Combine the two generators with masking to
            !ensure nonzero value.
        END FUNCTION ran
