program prueba_dinamico

    use korenaga_thermev_subs

    implicit none

    integer :: demid,k
    real :: t,a,n,lod,thp,eps
    real, parameter :: dt = 0.01
    character (len=1) :: chidfit
    character (len=100) :: fmt

    ai = afit(0.0) 
    af = afit(4.5)
    ma = (ai-af)/4.50
    LODi = LODfit(0.0)
    LODf = LODfit(4.5)
    mlod = (lodi-lodf)/4.5

    print *,'Insert DEM ID:'
    read *,demid

    do k=1,3
        ldem(k) = demid == k
    end do

    print *,'Generating data to evaluate obliquity ...'
    call generar_datos_oblicuidad()
    print *,'... ready!'

    write(chidfit,'(i1)') demid

    print *,'Creating out file ...'
    open(unit=50,file='../out/prueba_dinamico_DEMID_'//trim(chidfit)//'.out',status='unknown')
    print *,'... ready!'

    t = ti

    fmt = '(1x,f5.3,1x,f5.3,1x,f5.3,1x,f6.3)'

    do while (t.le.tf)

        call modelo_dinamico(t,a,n,lod,thp)
        
        eps = oblicuidad(t)
        
        write(50,fmt) t,a/a0,lod/LOD0,tidal_bulge(t)

        t = t + dt

    end do
    
end program prueba_dinamico