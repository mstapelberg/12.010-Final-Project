module MCdriver
    implicit none
    real, parameter :: PI = 3.1415926535
contains
subroutine MC(maxneutrons, maxinteractions, atomic_number, radius, a, phi, theta, r, x, y, z, neutron_abs, neutron_alpha,&
    neutron_esc_energy, neutron_scat)
    real, external :: sigma_t, sigma_gamma, sigma_naplha, sigma_s
    integer :: i, j !looping variables
    integer(kind = 4), intent(in) :: maxneutrons, maxinteractions, atomic_number 
    !We are assuming that we are only dealing with Nickel 58
    real(kind = 4), intent(in) :: radius, a 
    real(kind = 4) :: alpha 
!   real(kind = 4) :: E Energy Value for each particle
    real(kind = 4), dimension(maxneutrons, maxinteractions), intent(out) :: neutron_abs, neutron_alpha, neutron_esc_energy, x, y, z
    !The above arrays are empty, and then given a value of one (for neutron_abs and neutron_alpha) to 
    !denote at which interaction the neutron was absorbed or released an alpha, this will then 
    !let me find the corresponding x, y, z value 
    integer,dimension(maxneutrons), intent(out) :: neutron_scat
    real(kind = 4) :: interaction 
    integer :: n, m, seed = 1
    real(kind = 4), dimension(maxneutrons, maxinteractions) :: r, theta, phi, energy
    logical :: baby, life !used to see if the neutron is first born or not, or if the neutron is alive 
        !a = 1.1 !1.1meters
        alpha = ((atomic_number - 1.)/(atomic_number + 1.))**2.
        do i = 1, maxneutrons 
        !Initialize the Arrays
            baby = .True.
            r(i,1) = 0 
            theta(i,1) = 0
            phi(i,1) = 0
            energy(i,1) = 14E6
            do while((Life.EQV. .TRUE.) .AND. (j .LE. maxinteractions))
                if ( (radius .LT. a) .AND. (baby.EQV. .TRUE.)) then 
                !When neutrons are in the vacuum
                !Use SampleI for all of the neutrons to get
                !their initial direction. 
                    r(i,2) = 1.1
                    theta(i,2) = 2.*PI*rand(seed)
                    phi(i,2) = acos(2.*rand(seed) - 1.)
                    energy(i, 2) = 14E6
                else if ((radius .GE. 1.1) .AND. (radius .LT. 1.2)) then !Distance in meters from center
                !For when in the Nickel. Call the Scatter, absorb,
                !n-alpha
                    baby = .False.
                    !do j = 3, maxinteractions 
                        interaction = rand()
                        !SCATTERING 
                        if (interaction .LE. sigma_s(energy(i,j))/sigma_t(energy(i,j))) then 
                            energy(i,j) = energy(i,j-1) - energy(i,j-1)*(1.-alpha)
                            r(i,j) = -log(rand())/sigma_t(energy(i,j))
                            theta(i,j) = 2.*PI*rand()
                            phi(i,j) = acos(2.*rand()-1.)
                            neutron_scat(i) = neutron_scat(i) + 1. !Adds to the scatter counter for each neutron 
                            j = j+1
                        !NOW WE HAVE THE SIGMA_GAMMA    
                        else if (interaction .LE. sigma_gamma(energy(i,j))/sigma_t(energy(i,j))) then 
                            neutron_abs(i,j) = 1 !Sets a placeholder value for the neutron, to let us know which neutron was
                                                 !absorbed and at what energy
                            Life = .FALSE.
                        !NOW WE HAVE THE SIGMA_NALPHA
                        else if (interaction .LE. sigma_nalpha(energy(i,j))/sigma_t(energy(i,j))) then 
                            neutron_alpha(i,j) = 1 !Sets a placeholder value for the neutron, to let us know which neutron was
                                                   !absorbed and at what energy 
                            Life = .FALSE. !breaks out of the do while loop
                        end if
                    !end do 
                !HERE THE NEUTRON HAS LEFT THE SPHERE
                else !(r .GT. 1.2) or is outside the sphere 
                !Here we kill the neutron and save its value to an array     
                    neutron_esc_energy(i,j) = energy(i,j) !saves the current energy of that neutron to the neutron_esc_energy
                                                          !All other entries in that column will be NAN, which will make data
                                                          !analysis easier in numpy
                    Life = .FALSE.           
                end if 
            end do
        end do
        do n = 1, maxneutrons
            do m = 1, maxinteractions
                x(n,m) = r(n,m)*sin(theta(n,m))*cos(phi(n,m)) 
                y(n,m) = r(n,m)*sin(theta(n,m))*sin(phi(n,m))
                z(n,m) = r(n,m)*cos(theta(n,m))
            end do
        end do
        !print out arrays for debugging 
!write(*,90) (n, n=1,maxneutrons)
!90 format(' N', 10(1x, I2.2,' '))
!do m = 1, maxinteractions
!    write(*,100) m, (phi(n,m), n = 1, maxneutrons)
!    100 format(I4, 20(1x, F20.15))
!end do
end subroutine MC
!Here we convert the finished r, theta, phi arrays into x, y, z
!Will be implemented soon

!END OF MAIN PROGRAM
!**************************************************************************************************************
!THE FOLLOWING SUBROUTINES HAVE AN ISSUE, I CANNOT LOAD THE CSV FILES WITH THEIR ENERGY DEPENDENT CROSS SECTIONS
!MY PLAN IS TO LOAD IN THE CROSS SECTIONS, AND THEN HAVE THE INPUT BE THE ENERGY VALUE OF THE NEUTRON
!THE NEEDED SUBROUTINE WILL THEN SUBTRACT THE INPUT ENERGY FROM THE ABSOLUTE VALUE OF EACH ELEMENT 
!IN THE FIRST COLUMN OF THE CROSS SECTION TABLES (THE ENERGY VALUES CORRESPONDING TO THE CROSS SECTIONS)
!WHEN THE ABS(SIGMAX(1,i) - energy(i,j)) IS MINIMIZED, THE SUBROUTINE WILL OUTPUT THE CROSS SECTION THAT
!CORRESPONDS TO THAT ENERGY LEVEL. THIS MAKES SURE THAT WE HAVE AN ENERGY DEPENDENT CROSS SECTION AFTER
!ELASTIC SCATTERING.  \

            !SUBROUTINE TO OPEN THE SCATTERING CROSS SECTIONS 
            !These have been changed to functions for the time being cause they only output one value, 
            !And the name of their output is the name of the function
         real function sigma_t(E)
            real(kind = 4), intent(in) :: E !current energy of the neutron
            integer :: i, j ,k 
            !char(125) :: filename1 = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMAT.csv'

            !open(unit = 13, file = filename,status = 'old', form = 'formatted')
            !temp
            if (E.EQ. 0) then 
                sigma_t = 25
            else if ((E.GT.0) .AND. (E.LT.1E6)) then 
                sigma_t = 10
            else 
                sigma_t = 5
            end if
        end function sigma_t

            
        real function sigma_gamma(E)
            real(kind = 4), intent(in) :: E !Current energy of the neutron 
            integer :: i,j,k !looping variables
            !char(125) :: filename = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMAGAMMA.csv'
            !open(unit = 11, file = filename, )
            if (E.EQ. 0) then 
                sigma_gamma = 10
            else if ((E.GT.0) .AND. (E.LT.1E6)) then 
                sigma_gamma = 3
            else 
                sigma_gamma = 2
            end if
        end function sigma_gamma

            !SUBROUTINE TO OPEN THE SCATTERING CROSS SECTIONS 
        real function sigma_s(E)
            real(kind=4), intent(in) :: E !Current energy of the neutron
            integer :: i,j,k
            !char(125) :: filename = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMAS.csv'

            !open(unit = 10, file = filename, )
            if (E.EQ. 0) then 
                sigma_s = 10
            else if ((E.GT.0) .AND. (E.LT.1E6)) then 
                sigma_s = 3
            else 
                sigma_s = 2
            end if
        end function sigma_s 

        !SUBROUTINE TO OPEN THE NALPA CROSS SECTIONS
        real function sigma_nalpha(E)
            implicit none
            real(kind = 4), intent(in) :: E!Current energy of then neutron
            integer :: i,j,k
            !char(125) :: filename = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMANALPHA.csv'
            
            !open(unit = 12, file = filename, )
            if (E.EQ. 0) then 
                sigma_nalpha = 5
            else if ((E.GT.0) .AND. (E.LT.1E6)) then 
                sigma_nalpha = 2
            else 
                sigma_nalpha = 1
            end if
        end function sigma_nalpha
end module MCdriver
