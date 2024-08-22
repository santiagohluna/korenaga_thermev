program korenaga_thermev

    ! Main

    use korenaga_thermev_subs
    use integrators

    implicit none

    integer :: l,m,p,q,ok
    integer, parameter :: nvar = 1
    real :: t,dt,Ur0,Tpot0,Tpot,Ur,vel,ucal 
    real, parameter :: TW = 1.e-12      ! Terawatt
    real,dimension(1) :: y,dydx
    integer, parameter :: idunit = 20
    real, parameter :: tol = 1.e-6
    
    ! Get data necessary to evaluate the maximum plate thickness.
    call get_hdTC_data()

    ! Initialize
    call set_up_integration(idunit,ur0,Tpot0,dt)

    ! Set up summation arrays for evaluation of tidal heating
    nterms = 0
    do l=2,lmax
        do m = 0,l
            do p = 0,l
                do q = -qmax,qmax
                    nterms = nterms + 1
                end do
            end do
        end do
    end do
    allocate(asuma(nterms),stat=ok)
    if (ok.ne.0) then
        print *,'Error en la asignación de memoria de asuma'
        stop
    endif

!-----------------------------------------------------------------------
!   evaluación de los factores k_lm, funciones de la inclinación y de la
!   excentricidad
!-----------------------------------------------------------------------
    print *,'Evaluate k_lm factors...'
    call evalklm()
    print *,'... ready!'
    print *,'Generating data to evaluate obliquity ...'
    call generar_datos_oblicuidad()
    print *,'... ready!'
    print *,'Evaluate eccentricity functions ...'
    call evalge(e,gge)
    print *,'... ready!'

    ! Preliminary computations
    gamalf = exp(gammln(andrade + 1.0))
    gammac = gamalf*cos(0.5*andrade*pi)
    gammas = gamalf*sin(0.5*andrade*pi)
    ai = afit(0.0)
    af = afit(4.5)
    ma = (ai-af)/4.50
    LODi = LODfit(0.0)
    LODf = LODfit(4.5)
    mlod = (lodi-lodf)/4.5
    
    ! Calibration constat for convective flux
    Qcal = Qconv0/(Tref*sqrt(u(Tref)))

    ! Calibration constant for plate velocity
    ucal = u0*(Tref/Qconv0)**2

    Qrad0 = Ur0*Qconv0

    Tpot = Tpot0
    y(1) = Tpot0

    call derivs(t,y,dydx)
!              1       2       3       4       5       6       7       8
10  format(2x,f4.2,1x,f7.2,1x,f6.2,1x,f6.2,1x,f6.2,4x,f4.2,1x,f5.3,1x,f6.3)

    do while (t.le.tf)

        if (ltide) then
            Ur = (Qrad(t) + Qtidal(t,Tpot))/Qconv(Tpot)
        else
            Ur = Qrad(t)/Qconv(Tpot)
        end if
        
        vel = ucal*(Qconv(Tpot)/Tpot)**2
!                        1  2     3              4         5                 6  7
        write(idunit,10) t,Tpot,Qconv(Tpot)*TW,Qrad(t)*TW,Qtidal(t,Tpot)*TW,vel,Ur, &
!                               8           
                         DltetaL_ratio(Tpot)

        call odeint(y,t,t+dt,tol,dt,0.0,derivs,bsstep)
!        call rk4(y,dydx,t,dt,y,derivs) 

        Tpot = y(1)
        
        t = t + dt

    end do

end program korenaga_thermev