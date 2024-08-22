program prueba_delta_etaL

    use korenaga_thermev_subs

    implicit none
    
    real :: Tpot

    Tpot = 1350.0

    print *,'log(DetaL) = ',DltetaL_ratio(Tpot),printout

end program prueba_delta_etaL