program ftest
    implicit none
    real, parameter :: PI = 3.14159
    integer(kind = 4), parameter :: neutrons = 10, maxinteractions = 10
    real(kind = 4) :: radius, a
    real(kind = 4) :: E !Energy Value for each particle
    real(kind = 4), dimension(neutrons) :: neutron_abs, neutron_alpha, neutron_esc_energy
    real(kind = 4), dimension(maxinteractions) :: r, theta, phi
    logical :: baby !used to see if the neutron is first born or not 
        a = 1.1 !1.1meters
        do i = 1, neutrons 
            baby = .True.
            r(1) = 0 
            theta(1) = 0
            phi(1) = 0
            if (radius < a && baby = .TRUE.) then 
            !When neutrons are in the vacuum
            !Use SampleI for all of the neutrons to get
            !their initial direction. 
                r(2) = 1.1
                theta(2) = 2*PI*rand(seed)
                phi(2) = acos(2*rand(seed) - 1)

            else if (radius .geq. 1.1 && radius < 1.2) then !Distance in meters from center
            !For when in the Nickel. Call the Scatter, absorb,
            !n-alpha
                baby = .False.
                for j = 1, maxinteractions 
                    interaction = rand()
                    !SCATTERING 
                    if (interaction .leq. sigma_s(E)/sigma_t(E)) then 
                        neutron_char(E, )
                    






            else
            !Here we kill the neutron and save its value to an array     
            

            end if 
        end do





        contains 

            
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
