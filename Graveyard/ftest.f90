program ftest
    implicit none
    real, parameter :: PI = 3.14159
    integer(kind = 4), parameter :: maxneutrons = 10, maxinteractions = 10, atomic_number = 58 
    !We are assuming that we are only dealing with Nickel 58
    real(kind = 4) :: radius, a, alpha 
!   real(kind = 4) :: E Energy Value for each particle
    real(kind = 4), dimension(maxneutrons, maxinteractions) :: neutron_abs, neutron_alpha, neutron_esc_energy
    integer,dimension(maxneutrons) :: neutron_scat
    real(kind = 4), dimension(maxneutrons, maxinteractions) :: r, theta, phi, energy
    logical :: baby !used to see if the neutron is first born or not 
        a = 1.1 !1.1meters
        alpha = ((atomic_number - 1)/(atomic_number + 1 ))**2
        do i = 1, maxneutrons 
        !Initialize the Arrays
            baby = .True.
            r(i,1) = 0 
            theta(i,1) = 0
            phi(i,1) = 0
            energy(i,1) = 14E6
            do while((Life = .TRUE.) .AND. (j .LE. maxinteractions))
                if ( (radius .LT. a) .AND. (baby = .TRUE.)) then 
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
                            phi(i,j) = acos(2.*PI-1.)
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
        end program ftest




        contains 
            !SUBROUTINE TO OPEN THE SCATTERING CROSS SECTIONS 
         subroutine sigma_t(E)
            implicit none
            real(kind = 4), intent(in) :: E !current energy of the neutron
            integer :: i, j ,k 
            char(125) :: filename1 = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMAT.csv'

            open(unit = 13, file = filename,status = 'old', form = 'formatted')

        end subroutine sigma_t(E)

            
        subroutine sigma_gamma(E)
            implicit none 
            real(kind = 4), intent(in) :: E !Current energy of the neutron 
            integer :: i,j,k !looping variables
            char(125) :: filename = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMAGAMMA.csv'
            open(unit = 11, file = filename, )

        end subroutine sigma_gamma
            !SUBROUTINE TO OPEN THE SCATTERING CROSS SECTIONS 
        subroutine sigma_s(E)
            implicit none 
            real(kind=4), intent(in) :: E !Current energy of the neutron
            integer :: i,j,k
            char(125) :: filename = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMAS.csv'

            open(unit = 10, file = filename, )

        end subroutine sigma_s 

        !SUBROUTINE TO OPEN THE NALPA CROSS SECTIONS
        subroutine sigma_nalpha(E)
            implicit none
            real(kind = 4), intent(in) :: E!Current energy of then neutron
            integer :: i,j,k
            char(125) :: filename = '/home/myless/12.010-Final-Project/12.010-Final-Project/NI58SIGMANALPHA.csv'
            
            open(unit = 12, file = filename, )
        end subroutine sigma_nalpha

        

            
!**********************************************************************
!       SUBROUTINE TO SAMPLE R, THETA, AND PHI VALUES
        subroutine sample(neutrons, ninteractions)
            implicit none
            intent(in) :: neutrons, ninteractions
            !This should be a subroutine 
            integer :: i,j,k !looping varibles
            real,dimension(neutrons) :: phi, theta, r !sample arrays
            integer :: maxinteractions, maxneutrons
            integer, parameter :: seed = 1
            real :: sigma_t = 4
            
            !maxneutrons = 10
            !maxinteractions = 10

            do i = 1,neutrons
                do  j = 1, ninteractions
                    theta(j) = 2*PI*rand()
                    phi(j) =  acos(2*rand()-1)
                    r(j) = -log(rand())/sigma_t
                    !write(*,*) i, theta(j)
                end do 
            end do 
                
        end subroutine sample
!********************************************************************
!       SUBROUTINE THAT SAMPLES INITIAL THETA, AND PHI VALUES
        subroutine sample_i(neutrons)
            implicit none
            intent(in) :: neutrons
            integer :: i, j, k !looping variables
            real, dimension(neutrons) :: phi, theta
            integer, parameter :: seed = 1



        subroutine scatter(energy, )
        
        end subroutine scatter
        end program ftest
