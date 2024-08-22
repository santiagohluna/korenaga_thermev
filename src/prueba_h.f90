program prueba_h

    use korenaga_thermev_subs

    implicit none

    real :: Tpot
    real, parameter :: DTi = 1.0
    real, parameter :: Ti0 = 1200.0
    real, parameter :: Tif = Ti0 + 600.0

    open(unit=20,file='../out/prueba_h.out',status='unknown')

    call get_hdTC_data()

    Qcal = Qconv0/(Tref*sqrt(u(Tref))) ! Calibration constat for convective flux

    Tpot = Ti0

    do while (Tpot.le.Tif)

        write(20,*) Tpot,hdTC(Tpot),Qconv(Tpot)*1.e-12

        Tpot = Tpot + DTi
        
    end do
    
end program prueba_h
