        program mcint
            implicit none
            int*4 N = 10000
            real*8 DEN = 1
            real*8 SW = 0
            real*8 SWX = 0
            real*8 SWY = 0
            real*8 SWZ = 0
            real*8 VARW = 0
            real*8 VARX = 0
            real*8 VARY = 0
            real*8 VARZ = 0 
            real*8 VOL = 3.*7.*2. !Volume of the sampled region 

            int*4 J !Counting Variable 

            do J = 1, N
                X = 1. + 3. *RAN2(IDUM)
                T = -3. + 7. *RAN2(IDUM)
