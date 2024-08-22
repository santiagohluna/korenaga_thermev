module korenaga_thermev_subs

    use, intrinsic:: iso_fortran_env, only: stdin=>input_unit

    implicit none
    
    ! Model parameters
    ! ===================

    real, parameter :: pi = 4.0*atan(1.0)
    
    real, parameter :: day = 86400.0    ! Current day lenght in seconds
    real, parameter :: yr = 365.25*day  ! Current year lenght in seconds
    real, parameter :: Ga = 1.e9*yr     ! Gigayear
    real, parameter :: au = 149597870700.0 ! Astronomical unit

    real, parameter ::   CGU = 6.67408e-11
    real, parameter ::    Me = 5.9722e24
    real, parameter ::    Ml = 7.342e22
    real, parameter :: GMeMl = CGU*(Me+Ml)
    real, parameter ::  mred = Ml*(Me/(Me+Ml))

    real, parameter :: km = 1000.0              ! Kilometer
    real, parameter :: Rc = 3480.e3             ! Core's radius
    real, parameter :: Re = 6370.e3             ! Earth's radius
    real, parameter :: Dm = Re - Rc             ! Mantle thickness
    real, parameter :: Ve = 4.0*pi*(Re**3)/3.0  ! Earth's volume
    real, parameter :: rhot = Me/Ve             ! Mean Earth's density

    real, parameter ::   a0 = 60.142611e0*Re
    real, parameter :: LOD0 = 23.934468e0
    real, parameter :: eps0 = 23.2545*pi/180.0
    real, parameter ::   p0 = 50.467718e0*pi/(180.e0*3600.e0*yr)
    real :: e
    
    real, parameter :: Ct = 7.e27       ! Earth's heat capacity
    real, parameter :: cm = 1100.d0     ! specific heat of mantle
    real, parameter :: alpha = 2.4e-5   ! Mantle thermal expansivity
    real, parameter :: rhom = 3300.0    ! Mantle density
    real, parameter :: g = 9.8          ! Surface gravity
    real, parameter :: phidis = 0.6     ! disgregation point

    real, parameter :: LsD = 2.0        ! Convection cell aspect ratio
    real, parameter :: Rcpb = 200.e3    ! Radius of curvature for plate bending
    real, parameter :: etaL = 1.e23     ! Lithospheric viscosity
    real, parameter :: kapm = 0.86e-6     ! Mantle diffusivity
    
    real, parameter :: Eact = 300.e3    ! Mantle's viscosity activation energy
    real, parameter :: Rgas = 8.32      ! Gas constant
    real, parameter :: etaMref = 1.e21  ! Mantle's viscosity reference value
    real, parameter :: Tref = 1350.0    ! Reference temperature for current mantle viscosity
    real, parameter :: Rac = 1000.0     ! Critical Rayleigh number
    real, parameter :: Toff = 273.0     ! Temperature offset to convert from Celsius to Kelvin or viceversa

    real, parameter :: Cpe = 1.0/sqrt(pi)
    real, parameter :: CMvd = 3.0*(LsD + 2.5)
    real, parameter :: CSvd = 2.5

    real, parameter, dimension(4) :: h0 = (/ 0.372, 0.0164, 0.430, 0.181 /)
    real, parameter, dimension(4) :: lam = (/ 0.155, 0.985, 0.0495, 0.555 /)

    real, parameter :: Qconv0 = 36.e12  ! Current convective thermal flux
    real :: Qrad0
    real, parameter :: u0 = 4.0         ! Current plate velocity
    real, parameter :: Tsup = 0.0

    real :: Qcal                          ! Calibration constat for convective flux

    integer, parameter :: mpol = 3
    integer :: nmaxh,nmaxo
    real, allocatable, dimension(:) :: xa,ya
    real, dimension(1:60) :: xobl,yobl

    logical, dimension(8) :: lreo
    logical, dimension(3) :: ldem
    logical :: ltide
    logical :: lhmax

    integer :: lmax,qmax

    real, parameter :: ti = 0.0
    real, parameter :: tf = 4.5

    integer :: nterms
    real, allocatable :: asuma(:)

    integer, parameter :: lmaxp = 3
    integer, parameter :: qmaxp = 10
    real :: ffi(lmaxp,0:lmaxp,0:lmaxp),gge(lmaxp,0:lmaxp,-qmaxp:qmaxp),alm(lmaxp,0:lmaxp)

!   --------------------------------------------------------------------
!   parametros de los modelos dinamicos
!   --------------------------------------------------------------------
!   ---------
!   demid = 1
!   ---------
    real, parameter, dimension(4) :: ca = (/ -0.0130526, 0.0565061, -0.116111, 1.0 /)
    real, parameter, dimension(4) :: clod = (/ -0.0203081, 0.103267, -0.23361,1.0 /)
!   --------------------------------------------------------------------
!   demid = 2
!   ---------
    real, parameter :: dadt0 = 3.82e-2*ga/yr  ! Bills & Ray (1999)
    real, parameter :: dLODdt0 = 6.653        ! Williams & Boggs (2016)
!   --------------------------------------------------------------------
!   demid = 3
!   ---------
    real :: ai,af,LODi,LODf,ma,mLOD
    real, parameter :: dtconst = 0.5d0 ! Intervalo de tiempo en el que a y LOD permanecen constantes al principio de la evolución dinámica.
    real, parameter :: fcm = 4.5/(4.5-dtconst)
!   --------------------------------------------------------------------
!   Valor actual de la velocidad angular de rotación terrestre
!   --------------------------------------------------------------------
    real, parameter :: thp0 = 2.0*pi/(LOD0*3600.0)
!   --------------------------------------------------------------------

!   --------------------------------------------------------------------
!   parametros reologicos
!   --------------------------------------------------------------------
    real, parameter ::  rigdz = 1e11
    real, parameter ::   flex = 1.0/rigdz
    real, parameter ::  ksubb = Re/(CGU*Me*rhot)
    real, parameter ::  kflex = 0.2
    real, parameter ::   keta = 0.02
    real, parameter ::  andrade = 0.2
    real, parameter :: zandr0 = 1.0
    real :: gamalf,gammac,gammas
    real, parameter ::     k2 = 1.0
    real, parameter ::     qf = 1.0
    real, parameter :: deltat = 0.0
    real, parameter :: epsmay = 0.0
    
    real, parameter :: mu = 0.03

    real, parameter :: h_crit = 1.e-5

    real :: printout

contains

    ! Subroutines and functions
    ! =========================

    subroutine derivs(t, Tpot, dTpotdt)

        ! This subroutine computes the time derivative of the potential temperature of the Earth's interior. It corresponds to Eq. (1) of Korenaga (2006).

        implicit none

        real, intent(in) :: t
        real, intent(in), dimension(:) :: Tpot
        real, intent(out),dimension(:) :: dTpotdt 

        if (ltide) then
            dTpotdt(1) = (Qconv(Tpot(1)) - Qrad(t) - Qtidal(t,Tpot(1)))*Ga/Ct
        else
            dTpotdt(1) = (Qconv(Tpot(1)) - Qrad(t))*Ga/Ct
        end if
        
    end subroutine derivs

    function Qrad(t) result(retval)

        ! This function evaluates the time rate at which radiogenic heat is generated as a funcion of time.
        
        implicit none

        real, intent(in) :: t
        real :: retval,suma
        integer :: k

        suma = 0.0

        do k = 1,4
            suma = suma + h0(k)*exp(lam(k)*t)
        end do
        
        retval = Qrad0*suma
        
    end function 

    function Qconv(Tpot) result(retval)

        ! This function evaluates the time rate at which convective heat escapes from Earth as a function of the internal temperature. It corresponds to Eq. (19) of Korenaga (2006).

        implicit none

        real, intent(in) :: Tpot
        real :: retval
        
        retval = Qcal*Tpot*sqrt(u(Tpot))
        
    end function 

    function Qtidal(t,Tpot)

        implicit none 

        real, intent(in) :: t,Tpot
        real :: Qtidal,ftvf,Tubl,dubl,etavg

        dubl = hdTC(Tpot)*1.e3

        Tubl = Tpot*exp(g*alpha*dubl/cm)

        ftvf = (Rphi(Tubl,dubl)/Re)**3-(Rc/Re)**3

        etavg = etaM(Tpot)
        
        Qtidal = ftvf*pm(t,etavg)
        
    end function Qtidal

    function pm(t,eta)
!   --------------------------------------------------------------------
!   esta función calcula el calor generado por interacción de mareas
!   usando la expresión derivada por efroimsky y makarov (2014)
!   --------------------------------------------------------------------
        implicit none
!       --------------------------------------------------------------------
        integer :: j,l,m,p,q
        real,intent(in) :: t,eta
        real :: rsa,wlmpq,xlmpq,sumapm,rsal
        real :: pm,kr,ki,a,n,lod,thp,i
!       --------------------------------------------------------------------
        call modelo_dinamico(t,a,n,lod,thp)
        rsa = Re/a
!       --------------------------------------------------------------------
!       Evaluación de las funciones de la inclinación
!       --------------------------------------------------------------------
        i = oblicuidad(t)
        call evalfi(i,ffi)
!       --------------------------------------------------------------------
!       calculo de la tasa de produccion de calor por mareas
!       --------------------------------------------------------------------
        j = 1
        do l=2,lmax
            rsal = rsa**(2*l+1)
            do m = 0,l
                do p = 0,l
                    do q = -qmax,qmax
                        wlmpq = real(l-2*p+q)*n - real(m)*thp
                        xlmpq = abs(wlmpq)
                        call reologia(l,wlmpq,eta,kr,ki)
                        asuma(j) = rsal*alm(l,m)*ffi(l,m,p)*gge(l,p,q)*wlmpq*ki*sign(1.0,wlmpq)
                        j = j + 1
                    end do
                end do
            end do
        end do
        sumapm = sumar(j-1,asuma)
        pm = ((CGU*Ml)*(Ml/a))*sumapm
!       --------------------------------------------------------------------
    end function pm

    function Rphi(Tubl,dubl)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: Tubl,dubl
        real :: Rphi,Rini
!       ----------------------------------------------------------------
        Rini = Rast(Tubl,dubl)
        if (meltfdism(Rini,Tubl,dubl)*meltfdism(Re-dubl,Tubl,dubl).lt.0d0) then 
            Rphi = rtbis(meltfdism,Tubl,dubl,Rini,Re-dubl,1.e-3)
        else 
            Rphi = Rlit(Tubl,dubl)
        end if
!   --------------------------------------------------------------------        
    end function Rphi

    function Rast(Tubl,dubl)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: Tubl,dubl
        real :: rast
!       ----------------------------------------------------------------
!       determinacion de la distancia entre el centro de la tierra y el
!       limite inferior de la astenosfera
!       ----------------------------------------------------------------
        if(Tubl.gt.tsol(Re-dubl)) then
            Rast = rtbis(dTmTsol,Tubl,dubl,Rc,Re-dubl,1.e-3)
        else 
            Rast = Re - dubl
        end if
!   --------------------------------------------------------------------
    end function Rast

    function Rlit(Tubl,dubl)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: Tubl,dubl
        real :: Rlit
!       ----------------------------------------------------------------
!       Determinación de la distancia radial entre el centro de la 
!       tierra yla base de la litósfera o límite superior de la astenósfera
!       ----------------------------------------------------------------
        if(Tubl.gt.tsol(Re-dubl)) then
            Rlit = rtbis(dTcndTsol,Tubl,dubl,Re-dubl,Re,1.e-3)
        else 
            Rlit = Re - dubl
        end if
!   --------------------------------------------------------------------
    end function Rlit

    function tsol(r)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: r
        real :: tsol
!       ----------------------------------------------------------------
        if (pdr(r).le.20.0) then
            tsol = 1661.2*exp(log(pdr(r)/1.336 + 1.0)/7.437)
        else
            tsol = 2081.8*exp(log(pdr(r)/101.69 + 1.0)/1.226)
        end if
!   --------------------------------------------------------------------
    end function tsol

    function tliq(r)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: r
        real :: tliq
!       ----------------------------------------------------------------
        if (pdr(r).le.20.d0) then
            tliq = 1982.1*exp(log(pdr(r)/6.594 + 1.0)/5.374)
        else
            tliq = 2006.8*exp(log(pdr(r)/34.65 + 1.0)/1.844)
        end if
!   --------------------------------------------------------------------
    end function tliq

    function pdr(r)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: r
        real :: pdr,rkm
        real, parameter, dimension(3) :: coefm = (/ 4.35256e-6, -0.0901277, 397.012 /)
        real, parameter, dimension(3) :: coefc = (/ -1.52163e-5, -0.0144178, 367.767 /)
!       ----------------------------------------------------------------
        rkm = r*1.e-3
        if (r.le.Rc) then
            pdr = (coefc(1)*rkm + coefc(2))*rkm + coefc(3)
        else
            pdr = (coefm(1)*rkm + coefm(2))*rkm + coefm(3)
        end if
!   --------------------------------------------------------------------
    end function pdr

    function meltfdism(r,Tubl,dubl)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: r,Tubl,dubl
        real :: meltfdism
!       ----------------------------------------------------------------
        meltfdism = fmeltm(Tubl,dubl,r) - phidis
!   --------------------------------------------------------------------
    end function meltfdism

    function fmeltm(r,Tubl,dum)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: r,dum,Tubl
        real :: fmeltm,tdr
!       ----------------------------------------------------------------
        tdr = Tm(Tubl,dum,r)
        if((tdr-tsol(r))*(tliq(r)-tdr).gt.0.0) then 
            fmeltm = (tdr - tsol(r))/(tliq(r) - tsol(r))
        else if (tdr.ge.tliq(r)) then
            fmeltm = 1.0
        else 
            fmeltm = 0.0
        end if
!       ----------------------------------------------------------------
    end function fmeltm

    function tm(Tubl,dum,r)
!   --------------------------------------------------------------------
!   perfil adiabático de temperatura del manto
!   --------------------------------------------------------------------
    implicit none
!   --------------------------------------------------------------------
    real, intent(in) :: r,dum,Tubl
    real :: tm
!   --------------------------------------------------------------------
    tm = Tubl*(1.0 + alpha*g*(Re-dum-r)/cm)
!   --------------------------------------------------------------------
    end function tm

    function tcondubl(Tubl,dum,r)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) ::Tubl,dum,r
        real :: tcondubl,Rin
!       ----------------------------------------------------------------
        Rin = Re - dum
        tcondubl = (Re*Tsup*(r-Rin) + Rin*Tubl*(Re-r))/(dum*r)
!   --------------------------------------------------------------------
    end function tcondubl

    function dTmTsol(r,Tubl,dubl)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: r,Tubl,dubl
        real :: dTmTsol
!       ----------------------------------------------------------------
        dTmTsol = tm(Tubl,dubl,r)-tsol(r)
!   --------------------------------------------------------------------
    end function dTmTsol

    function dTcndTsol(r,Tubl,dubl)
!   --------------------------------------------------------------------        
        implicit none
!       ----------------------------------------------------------------        
        real, intent(in) :: r,Tubl,dubl
        real :: dTcndTsol
!       ----------------------------------------------------------------
        dTcndTsol = tcondubl(Tubl,dubl,r)-tsol(r)
!   --------------------------------------------------------------------        
    end function dTcndTsol

    function u(Tpot)

        implicit none

        real :: u,num,den
        real, intent(in) :: Tpot

        num = Cpe*alpha*rhom*g*Tpot*Dm*hdTC(Tpot)*km
        den = CMvd*etaM(Tpot) + CSvd*etaL*(hdTC(Tpot)*km/Rcpb)**3

        u = num/den

    end function u

    function etaM(Tpot)

        implicit none

        real, intent(in) :: Tpot
        real :: etaM

        etaM = etaMref*exp((Eact/Rgas)*(1.0/(Tpot + Toff) - 1.0/(Tref + Toff)))
        
    end function etaM

    function hdTC(Tpot) result(hm)

        implicit none

        real, intent(in) :: Tpot
        real :: hm,y,dy
        integer :: j,k

        if (lhmax) then

            call hunt(xa,nmaxh-1,Tpot,j)
        
            k = min(max(j-(mpol-1)/2,1),nmaxh+1-mpol)
            
            call polint(xa(k),ya(k),mpol,Tpot,y,dy)
    
            hm = y

        else

            hm = 120.0 + 0.275*(Tpot - 1400.0)

        end if
        
    end function hdTC

    subroutine get_hdTC_data()

        implicit none

        integer :: k,feof
        integer, parameter :: kmax = 100
        real,dimension(kmax) :: Tpot,hh
        real :: x,y
        logical :: ex

        ! Open file

        inquire(file='../data/hdTC.dat',exist=ex)

        if (ex) then
            print *,'Data file found!'
            open(unit=50,file='../data/hdTC.dat',status='old')
        else
            print *,'data file does not exist'
            stop
        end if

        ! Initialize variables
        k = 1
        feof = 0

        print *, 'Processing data file hdTC.dat'

        do while(feof == 0)
            read(50,*,iostat=feof) x,y
            if(feof > 0) then
                print *,'Check data file hdTC.dat'
                stop
            else if(feof < 0) then
                print *, 'Ready!'
            else
                Tpot(k) = x
                hh(k) = y
                k = k + 1
            end if
        end do

        nmaxh = k - 1

        allocate(xa(nmaxh))
        allocate(ya(nmaxh))

        do k=1,nmaxh-1
            xa(k) = Tpot(k)
            ya(k) = hh(k)
        end do

        close(50)
        
    end subroutine get_hdTC_data

    function DltetaL_ratio(Tpot) result(res)

        ! Plate tectonics propension from Korenaga (2010).

        implicit none

        real, intent(in) :: Tpot
        real :: res,theta,gam,DltetaL,DetaLcrit

        theta = Eact*(Tpot-Tsup)/(Rgas*Tpot**2)

        gam = mu/(alpha*(Tpot-Tsup))

        printout = gam

        DltetaL = exp(0.327*theta*(gam**0.647))

        DetaLcrit = 0.25*sqrt(Ra(Tpot))

        res = log10(DltetaL/DetaLcrit)
        
    end function DltetaL_ratio

    function tidal_bulge(t) result(res)

        ! Tectonic activity measure from Zanazzi & Triaud (2019).

        implicit none 

        real, intent(in) :: t
        real :: res,eff_rig,a,n,lod,thp,h_tidal

        eff_rig = 0.5*19.0*rigdz*ksubb

        call modelo_dinamico(t,a,n,lod,thp)

        h_tidal = (Re/a)**3/(1.0 + eff_rig)

        res = log10(h_crit/h_tidal)
        
    end function tidal_bulge

    function h(Tpot)

        implicit none

        real :: h,Nu,theta,gam,DltetaL
        real, intent(in) :: Tpot

        theta = Eact*(Tpot-Tsup)/(Rgas*Tpot**2)

        gam = mu/(alpha*(Tpot-Tsup))

        DltetaL = exp(0.327*theta*gam**0.647)

        Nu = 2.0*(Ra(Tpot)/(Rac*DltetaL))**(1.0/3.0)

        h = Dm/Nu

    end function h

    function Ra(Tpot)

        implicit none

        real, intent(in) :: Tpot
        real :: Ra
    
        Ra = (alpha*rhom*g*(Tpot-Tsup)*Dm**3)/(kapm*etaM(Tpot))
        
    end function Ra

    ! function mu(Tpot)
    !     real, intent(in) :: Tpot
    !     real :: mu
    !     real, parameter, dimension(2) :: fitcof = (/ 9.77205e-5, -0.0790665 /)

    !     mu = fitcof(1)*Tpot + fitcof(2)
        
    ! end function mu

!=======================================================================
    function afit(t)
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: t
        real :: afit
!       ----------------------------------------------------------------
        afit = (t*(t*(ca(1)*t + ca(2)) + ca(3)) + ca(4))*a0
!   --------------------------------------------------------------------
    end function afit
!=======================================================================
    function LODfit(t)
!   --------------------------------------------------------------------        
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: t
        real :: LODfit
!       ----------------------------------------------------------------
        LODfit = (t*(t*(clod(1)*t + clod(2)) + clod(3)) + clod(4))*LOD0
!   --------------------------------------------------------------------    
    end function LODfit
!=======================================================================
    function depsda(t, x)
    !   --------------------------------------------------------------------
    !   Esta subrutina calcula la derivada temporal de la oblicuidad
    !   terrestre utilizando la expresión dada por la ec. (A.6) del trabajo
    !   de Farhat et al. (2022).
    !   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        real, intent(in) :: t,x
        real :: depsda
        real :: a,lod,n,thp,Cthp,frac1,frac2
        real, parameter :: Cthp0 = 0.3306947357075918999972*Me*Re**2
        real, parameter :: kf2 = 0.93
!       ----------------------------------------------------------------    
        call modelo_dinamico(t,a,n,lod,thp)
!       ----------------------------------------------------------------
        Cthp = Cthp0 + 2.0*kf2*(Re**5)*(thp**2 - thp0**2)/(9.0*CGU)
!       ----------------------------------------------------------------
        frac1 = 0.25*mred*n*a/(Cthp*thp)
        frac2 = (thp*cos(x) - 2.0*n)/(thp*cos(x) - n)
!       ----------------------------------------------------------------
        depsda = frac1*sin(x)*frac2
    !   --------------------------------------------------------------------
    end function depsda
!=======================================================================
    subroutine generar_datos_oblicuidad()

        implicit none

        integer :: k
        real :: t,a1,a2,eps,da,lod1,lod2,thp,n
        real :: h1,h2,h3,h4,da2,da6
        real, parameter :: dt = 0.1e0
        real, parameter :: dt2 = 0.5e0*dt

        k = 1
        t = 0.0
        eps = eps0
    
        do while (t.le.5.0)
            call modelo_dinamico(t,a1,n,lod1,thp)
            call modelo_dinamico(t+dt,a2,n,lod2,thp)
            da = a2 - a1
            da2 = 0.5*da
            da6 = da/6.0
!           -----------------------------------------------------------------------
!           Integración de la ec. (A.6) de Farhat et al. (2020) usando el método
!           de Runge-Kutta de orden 4.
!           -----------------------------------------------------------------------
            h1 = depsda(t,eps)
            h2 = depsda(t+dt2,eps+h1*da2)
            h3 = depsda(t+dt2,eps+h2*da2)
            h4 = depsda(t+dt,eps+h3*da)
            eps = eps + da6*(h1 + 2.0*(h2 + h3) + h4)
!           -----------------------------------------------------------------------
            xobl(k) = t
            yobl(k) = eps
!           -----------------------------------------------------------------------
            t = t + dt
            k = k + 1
        end do

        nmaxo = k - 1
    
    end subroutine generar_datos_oblicuidad
!=======================================================================
    function oblicuidad(x)

        real, intent(in) :: x
        real :: y,dy,oblicuidad
        integer :: j,k
    
        call hunt(xobl,nmaxo,x,j)
        k = min(max(j-(mpol-1)/2,1),nmaxo+1-mpol)
        call polint(xobl(k),yobl(k),mpol,x,y,dy)
        
        oblicuidad = y

    end function oblicuidad
!=======================================================================
    subroutine modelo_dinamico(t,a,n,lod,thp)
!   --------------------------------------------------------------------
!   calculo de los valores de a y de lod
!   --------------------------------------------------------------------
        implicit none
!       --------------------------------------------------------------------
        real, intent(in) :: t
        real, intent(out) :: a,n,lod,thp
!       --------------------------------------------------------------------
        if (ldem(1)) then
!       --------------------------------------------------------------------
                a = afit(t)
            lod = LODfit(t)
!       ---------------------------------------------------------------------
        else if (ldem(2)) then
!       --------------------------------------------------------------------
            a = ma*t + a0
            lod = mlod*t + lod0
!       --------------------------------------------------------------------
        else if (ldem(3)) then
!       --------------------------------------------------------------------
            if ((4.5e0-t)*(t-4.5e0+dtconst).ge.0.e0) then
                a = ai
                LOD = LODi
            else 
                a = fcm*ma*t + af
                LOD = fcm*mLOD*t + LODf
            end if
!       --------------------------------------------------------------------
        else
            print *,'error en el identificador del modelo dinamico'
            return
        end if
!       --------------------------------------------------------------------
            n = sqrt(GMeMl/a)/a
        thp = 2.e0*pi/(lod*3600.e0)
!   --------------------------------------------------------------------
    end subroutine modelo_dinamico
!=======================================================================
    subroutine reologia(l,w,eta,kr,ki)
!   --------------------------------------------------------------------
    implicit none
!   --------------------------------------------------------------------
    real, intent(in) :: w,eta
    real, intent(out) :: kr,ki
    real :: x,zeta,bsubl,facl,denkl,numki,numkr,xpr,xpi,lag
    integer,intent(in) :: l
!   --------------------------------------------------------------------
    x = abs(w)
    zeta = zandr0*(100.0*exp(-x/0.2) + 1.0)
    bsubl = real(2*l**2+4*l+3)*ksubb/real(l)
    facl = 1.5/real(l-1)
    if (lreo(1)) then
!    puramente elastico
        kr = facl*flex/(flex + bsubl)
        ki = 0.e0
        return
    else if (lreo(2)) then
!      puramente viscoso
        denkl = (bsubl*x*eta)**2 + 1.e0
        numkr = 1.e0
        numki = -bsubl*x*eta
    else if (lreo(3)) then
!      maxwell
        denkl = ((eta*x)**2)*(flex + bsubl)**2 + 1.e0
        numkr = (flex + bsubl)*flex*(eta*x)**2 + 1.e0
        numki = -bsubl*eta*x
    else if (lreo(4)) then
!      burgers
        xpr = flex*eta*x*(1.e0 + kflex/((kflex*keta*flex*eta*x)**2 + 1.e0))
        xpi = - (1.e0 + keta*(kflex*flex*eta*x)**2/((kflex*keta*flex*eta*x)**2 + 1.e0))         
        numkr = (xpr + bsubl*eta*x)*xpr + xpi**2
        numki = bsubl*eta*x*xpi
        denkl = (xpr + bsubl*eta*x)**2 + xpi**2             
    else if (lreo(5)) then
!     andrade
        xpr = flex*eta*x + ((flex*eta*x)**(1.e0-andrade))*gammac/(zeta**andrade)
        xpi = - (1.e0 + ((flex*eta*x)**(1.e0-andrade))*gammas/(zeta**andrade))
        numkr = (xpr + bsubl*eta*x)*xpr + xpi**2
        numki = bsubl*eta*x*xpi
        denkl = (xpr + bsubl*eta*x)**2 + xpi**2
    else if (lreo(6)) then
!      kaula (1964)
        kr = 0.e0
        ki = k2/qf
        return
    else if (lreo(7)) then
!      singer-mignard
        kr = 0.e0
        ki = k2*x*deltat
        return
    else if (lreo(8)) then
!      efroimsky & lainey (2007)
        lag = epsmay*(epsmay*x)**(-alpha-1.0e0)
        kr = 0.e0
        ki = k2*lag*x
        return
    else
        print *,'error en el identificador de la reologia'
        call abort
    end if
    kr =  facl*numkr/denkl
    ki = -facl*numki/denkl
!   --------------------------------------------------------------------
    return
!   --------------------------------------------------------------------
    end subroutine reologia
!=======================================================================
    subroutine evalklm()
!   --------------------------------------------------------------------
        implicit none
!     ------------------------------------------------------------------
!     evaluacion de las funciones k_lm
!     ------------------------------------------------------------------
!     l=2
        alm(2,0) = 1.0e0
        alm(2,1) = 0.333333333333333e0
        alm(2,2) = 0.0833333333333333e0
!     l=3
        alm(3,0) = 1.0e0
        alm(3,1) = 0.166666666666667e0
        alm(3,2) = 0.0166666666666667e0
        alm(3,3) = 0.00277777777777778e0
!   --------------------------------------------------------------------
    end subroutine evalklm
!=======================================================================
    subroutine evalfi(inc,f2i)
!   --------------------------------------------------------------------
        implicit none
!   --------------------------------------------------------------------
        real, intent(in) :: inc
        real :: ci
        real, intent(out) :: f2i(lmaxp,0:lmaxp,0:lmaxp)
!   --------------------------------------------------------------------
        ci = cos(inc)
!     l=2
        f2i(2,0,0) = 0.140625e0*(-1.0e0 + ci**2)**2
        f2i(2,0,1) = 0.5625e0*(-0.333333333333333e0 + ci**2)**2
        f2i(2,0,2) = 0.140625e0*(-1.0e0 + ci**2)**2
        f2i(2,1,0) = -0.5625e0*(-1.0e0 + ci)*(1.0e0 + ci)**3
        f2i(2,1,1) = 2.25e0*ci**2*(1.0e0 - ci**2)
        f2i(2,1,2) = -0.5625e0*(-1.0e0 + ci)**3*(1.0e0 + ci)
        f2i(2,2,0) = 0.5625e0*(1.0e0 + ci)**4
        f2i(2,2,1) = 2.25e0*(-1.0e0 + ci)**2*(1.0e0 + ci)**2
        f2i(2,2,2) = 0.5625e0*(-1.0e0 + ci)**4
!     l=3
        f2i(3,0,0) = -0.09765625e0*(-1.0e0 + ci**2)**3
        f2i(3,0,1) = 0.87890625e0*(1.0e0 - ci**2)*(-0.2e0 + ci**2)**2
        f2i(3,0,2) = 0.87890625e0*(1.0e0 - ci**2)*(-0.2e0 + ci**2)**2
        f2i(3,0,3) = -0.09765625e0*(-1.0e0 + ci**2)**3
        f2i(3,1,0) = 3.515625e0*(1.0e0 + ci)**4*(0.5e0 - ci + 0.5e0*ci**2) &
        **2/(-1.0e0 + ci)**2
        f2i(3,1,1) = 7.91015625e0*(1.0e0 + ci)**2*(0.0666666666666667e0 + &
        0.666666666666667e0*ci - ci**2)**2
        f2i(3,1,2) = 7.91015625e0*(-1.0e0 + ci)**2*(-0.0666666666666667e0 &
        + 0.666666666666667e0*ci + ci**2)**2
        f2i(3,1,3) = 3.515625e0*(-1.0e0 + ci)**4*(0.5e0 + ci + &
        0.5e0*ci**2)**2/(1.0e0 + ci)**2
        f2i(3,2,0) = -3.515625e0*(-1.0e0 + ci)*(1.0e0 + ci)**5
        f2i(3,2,1) = 31.640625e0*(1.0e0 - ci)*(-0.333333333333333e0 + ci) &
        **2*(1.0e0 + ci)**3
        f2i(3,2,2) = -31.640625e0*(-1.0e0 + ci)**3*(0.333333333333333e0 + &
        ci)**2*(1.0e0 + ci)
        f2i(3,2,3) = -3.515625e0*(-1.0e0 + ci)**5*(1.0e0 + ci)
        f2i(3,3,0) = 3.515625e0*(1.0e0 + ci)**6
        f2i(3,3,1) = 31.640625e0*(-1.0e0 + ci)**2*(1.0e0 + ci)**4
        f2i(3,3,2) = 31.640625e0*(-1.0e0 + ci)**4*(1.0e0 + ci)**2
        f2i(3,3,3) = 3.515625e0*(-1.0e0 + ci)**6
!   --------------------------------------------------------------------
    end subroutine evalfi
!=======================================================================
    subroutine evalge(exc,g2e)
!   --------------------------------------------------------------------
    implicit none
!   --------------------------------------------------------------------
    real, intent(in) :: exc
    real :: e2,e3,e4,e5,e6,e7,e8,e9,e10
    real, intent(out) :: g2e(lmaxp,0:lmaxp,-qmaxp:qmaxp)
!   --------------------------------------------------------------------
        e2 = exc*exc
        e3 = e2*exc
        e4 = e2*e2
        e5 = e4*exc
        e6 = e4*e2
        e7 = e6*exc
        e8 = e6*e2
        e9 = e8*exc
        e10 = e8*e2
!   --------------------------------------------------------------------
!     l=2
!     ---
        g2e(2,0,-10) = 0.e0
        g2e(2,0,-9) = 0.e0
        g2e(2,0,-8) = 0.e0
        g2e(2,0,-7) = 0.e0
        g2e(2,0,-6) = 0.e0
        g2e(2,0,-5) = 0.0040045166015625e0*e10
        g2e(2,0,-4) = e8*(0.00173611111111111e0 + 0.00243055555555556e0*e2)
        g2e(2,0,-3) = e6*(0.000434027777777778e0 + e2*( &
        0.000596788194444444e0 + 0.000629679361979166e0*e2))
        g2e(2,0,-2) = 0.e0
        g2e(2,0,-1) = e2*(0.25e0 + e2*(-0.0625e0 + e2*( &
        0.0169270833333333e0 + e2*(0.00613064236111114e0 + &
        0.0053690592447917e0*e2))))
        g2e(2,0,0) = 1.0e0 + e2*(-5.0e0 + e2*(7.875e0 + e2*( &
        -4.30555555555556e0 + e2*(1.25043402777778e0 - 0.181302083333333e0*e2))))
        g2e(2,0,1) = e2*(12.25e0 + e2*(-53.8125e0 + e2*( &
        85.83984375e0 + e2*(-64.76318359375e0 + 28.4081359863281e0*e2))))
        g2e(2,0,2) = e4*(72.25e0 + e2*(-325.833333333333e0 + e2*( &
        580.215277777778e0 - 547.1625e0*e2)))
        g2e(2,0,3) = e6*(309.906684027778e0 + e2*(-1491.08208550347e0 &
        + 2986.78270975749e0*e2))
        g2e(2,0,4) = e8*(1109.72265625e0 - 5757.64921875e0*e2)
        g2e(2,0,5) = 3536.12958502876e0*e10
        g2e(2,0,6) = 0.e0
        g2e(2,0,7) = 0.e0
        g2e(2,0,8) = 0.e0
        g2e(2,0,9) = 0.e0
        g2e(2,0,10) = 0.e0
        g2e(2,1,-10) = 0.e0
        g2e(2,1,-9) = 0.e0
        g2e(2,1,-8) = 0.e0
        g2e(2,1,-7) = 0.e0
        g2e(2,1,-6) = 0.e0
        g2e(2,1,-5) = 47.9664459228516e0*e10
        g2e(2,1,-4) = e8*(23.16015625e0 + 7.76015625e0*e2)
        g2e(2,1,-3) = e6*(10.97265625e0 + e2*(10.17041015625e0 + 18.3712188720703e0*e2))
        g2e(2,1,-2) = e4*(5.0625e0 + e2*(7.875e0 + e2*(12.9765625e0 + 18.7921875e0*e2)))
        g2e(2,1,-1) = e2*(2.25e0 + e2*(5.0625e0 + e2*(8.96484375e0 + &
        e2*(13.86865234375e0 + 19.7799133300781e0*e2))))
        g2e(2,1,0) = 1.0e0 + e2*(3.0e0 + e2*(6.0e0 + e2*(10.0e0 + exc &
        **2*(15.0e0 + 21.0e0*e2))))
        g2e(2,1,1) = e2*(2.25e0 + e2*(5.0625e0 + e2*(8.96484375e0 + &
        e2*(13.86865234375e0 + 19.7799133300781e0*e2))))
        g2e(2,1,2) = e4*(5.0625e0 + e2*(7.875e0 + e2*(12.9765625e0 + 18.7921875e0*e2)))
        g2e(2,1,3) = e6*(10.97265625e0 + e2*(10.17041015625e0 + 18.3712188720703e0*e2))
        g2e(2,1,4) = e8*(23.16015625e0 + 7.76015625e0*e2)
        g2e(2,1,5) = 47.9664459228516e0*e10
        g2e(2,1,6) = 0.e0
        g2e(2,1,7) = 0.e0
        g2e(2,1,8) = 0.e0
        g2e(2,1,9) = 0.e0
        g2e(2,1,10) = 0.e0
        g2e(2,2,-10) = 0.e0
        g2e(2,2,-9) = 0.e0
        g2e(2,2,-8) = 0.e0
        g2e(2,2,-7) = 0.e0
        g2e(2,2,-6) = 0.e0
        g2e(2,2,-5) = 3536.12958502876e0*e10
        g2e(2,2,-4) = e8*(1109.72265625e0 - 5757.64921875e0*e2)
        g2e(2,2,-3) = e6*(309.906684027778e0 + e2*(-1491.08208550347e0 &
        + 2986.78270975749e0*e2))
        g2e(2,2,-2) = e4*(72.25e0 + e2*(-325.833333333333e0 + e2*( &
        580.215277777778e0 - 547.1625e0*e2)))
        g2e(2,2,-1) = e2*(12.25e0 + e2*(-53.8125e0 + e2*( &
        85.83984375e0 + e2*(-64.76318359375e0 + 28.4081359863281e0*e2))))
        g2e(2,2,0) = 1.0e0 + e2*(-5.0e0 + e2*(7.875e0 + e2*( &
        -4.30555555555556e0 + e2*(1.25043402777778e0 - &
        0.181302083333333e0*e2))))
        g2e(2,2,1) = e2*(0.25e0 + e2*(-0.0625e0 + e2*( &
        0.0169270833333333e0 + e2*(0.00613064236111114e0 + 0.0053690592447917e0*e2))))
        g2e(2,2,2) = 0.e0
        g2e(2,2,3) = e6*(0.000434027777777778e0 + e2*( &
        0.000596788194444444e0 + 0.000629679361979166e0*e2))
        g2e(2,2,4) = e8*(0.00173611111111111e0 + 0.00243055555555556e0*exc**2)
        g2e(2,2,5) = 0.0040045166015625e0*e10
        g2e(2,2,6) = 0.e0
        g2e(2,2,7) = 0.e0
        g2e(2,2,8) = 0.e0
        g2e(2,2,9) = 0.e0
        g2e(2,2,10) = 0.e0
!     ---
!     l=3
!     ---
        g2e(3,0,-10) = 0.e0
        g2e(3,0,-9) = 0.e0
        g2e(3,0,-8) = 0.e0
        g2e(3,0,-7) = 0.e0
        g2e(3,0,-6) = 0.e0
        g2e(3,0,-5) = 6.94444444444445e-5*e10
        g2e(3,0,-4) = e8*(6.7816840277778e-6 + 1.35633680555556e-5*e2)
        g2e(3,0,-3) = 0.e0
        g2e(3,0,-2) = e4*(0.015625e0 + e2*(0.00520833333333334e0 + exc** &
        2*(0.00490993923611111e0 + 0.00393880208333333e0*e2)))
        g2e(3,0,-1) = e2*(1.0e0 + e2*(-2.5e0 + e2*( &
        1.85416666666667e0 + e2*(-0.524305555555555e0 + 0.136197916666667e0*e2))))
        g2e(3,0,0) = 1.0e0 + e2*(-12.0e0 + e2*(49.21875e0 + e2*( &
        -83.21875e0 + e2*(68.0408935546875e0 - 31.233955078125e0*e2))))
        g2e(3,0,1) = e2*(25.0e0 + e2*(-220.0e0 + e2*( &
        736.916666666667e0 + e2*(-1221.72222222222e0 + 1147.2109375e0*exc**2))))
        g2e(3,0,2) = e4*(252.015625e0 + e2*(-2027.36979166667e0 + e2 &
        *(6597.1491156684e0 - 11511.4770507813e0*e2)))
        g2e(3,0,3) = e6*(1660.5625e0 + e2*(-13126.59375e0 + 43691.82890625e0*e2))
        g2e(3,0,4) = e8*(8504.77816433377e0 - 68154.558710395e0*e2)
        g2e(3,0,5) = 36828.8084027778e0*e10
        g2e(3,0,6) = 0.e0
        g2e(3,0,7) = 0.e0
        g2e(3,0,8) = 0.e0
        g2e(3,0,9) = 0.e0
        g2e(3,0,10) = 0.e0
        g2e(3,1,-10) = 0.e0
        g2e(3,1,-9) = 0.e0
        g2e(3,1,-8) = 0.e0
        g2e(3,1,-7) = 0.e0
        g2e(3,1,-6) = 0.e0
        g2e(3,1,-5) = 14.0312673611111e0*e10
        g2e(3,1,-4) = e8*(7.18072509765625e0 + 23.6063720703125e0*e2)
        g2e(3,1,-3) = e6*(3.67361111111111e0 + e2*(14.2152777777778e0 + 36.3644097222222e0*e2))
        g2e(3,1,-2) = e4*(1.890625e0 + e2*(8.421875e0 + e2*( &
        23.4019368489583e0 + 51.6582790798611e0*e2)))
        g2e(3,1,-1) = e2*(1.0e0 + e2*(5.0e0 + e2*(15.0e0 + e2*(35.0e0 + 70.0e0*e2))))
        g2e(3,1,0) = 1.0e0 + e2*(4.0e0 + e2*(11.46875e0 + e2*( &
        26.4756944444444e0 + e2*(53.2151557074653e0 + 96.8403244357639e0*e2))))
        g2e(3,1,1) = e2*(9.0e0 + e2*(16.5e0 + e2*(38.1875e0 + e2*( &
        71.4791666666667e0 + 124.026996527778e0*e2))))
        g2e(3,1,2) = e4*(43.890625e0 + e2*(32.296875e0 + e2*( &
        97.048095703125e0 + 149.09169921875e0*e2)))
        g2e(3,1,3) = e6*(164.694444444444e0 + e2*(-13.3680555555555e0 + 254.317795138889e0*e2))
        g2e(3,1,4) = e8*(532.960510253906e0 - 416.388549804688e0*e2)
        g2e(3,1,5) = 1567.17015625e0*e10
        g2e(3,1,6) = 0.e0
        g2e(3,1,7) = 0.e0
        g2e(3,1,8) = 0.e0
        g2e(3,1,9) = 0.e0
        g2e(3,1,10) = 0.e0
        g2e(3,2,-10) = 0.e0
        g2e(3,2,-9) = 0.e0
        g2e(3,2,-8) = 0.e0
        g2e(3,2,-7) = 0.e0
        g2e(3,2,-6) = 0.e0
        g2e(3,2,-5) = 1567.17015625e0*e10
        g2e(3,2,-4) = e8*(532.960510253906e0 - 416.388549804688e0*e2)
        g2e(3,2,-3) = e6*(164.694444444444e0 + e2*(-13.3680555555555e0 + 254.317795138889e0*e2))
        g2e(3,2,-2) = e4*(43.890625e0 + e2*(32.296875e0 + e2*( &
        97.048095703125e0 + 149.09169921875e0*e2)))
        g2e(3,2,-1) = e2*(9.0e0 + e2*(16.5e0 + e2*(38.1875e0 + e2* &
        (71.4791666666667e0 + 124.026996527778e0*e2))))
        g2e(3,2,0) = 1.0e0 + e2*(4.0e0 + e2*(11.46875e0 + e2*( &
        26.4756944444444e0 + e2*(53.2151557074653e0 + 96.8403244357639e0*e2))))
        g2e(3,2,1) = e2*(1.0e0 + e2*(5.0e0 + e2*(15.0e0 + e2*( 35.0e0 + 70.0e0*e2))))
        g2e(3,2,2) = e4*(1.890625e0 + e2*(8.421875e0 + e2*( &
        23.4019368489583e0 + 51.6582790798611e0*e2)))
        g2e(3,2,3) = e6*(3.67361111111111e0 + e2*(14.2152777777778e0 + 36.3644097222222e0*e2))
        g2e(3,2,4) = e8*(7.18072509765625e0 + 23.6063720703125e0*e2)
        g2e(3,2,5) = 14.0312673611111e0*e10
        g2e(3,2,6) = 0.e0
        g2e(3,2,7) = 0.e0
        g2e(3,2,8) = 0.e0
        g2e(3,2,9) = 0.e0
        g2e(3,2,10) = 0.e0
        g2e(3,3,-10) = 0.e0
        g2e(3,3,-9) = 0.e0
        g2e(3,3,-8) = 0.e0
        g2e(3,3,-7) = 0.e0
        g2e(3,3,-6) = 0.e0
        g2e(3,3,-5) = 36828.8084027778e0*e10
        g2e(3,3,-4) = e8*(8504.77816433377e0 - 68154.558710395e0*e2)
        g2e(3,3,-3) = e6*(1660.5625e0 + e2*(-13126.59375e0 + 43691.82890625e0*e2))
        g2e(3,3,-2) = e4*(252.015625e0 + e2*(-2027.36979166667e0 + exc** &
        2*(6597.1491156684e0 - 11511.4770507813e0*e2)))
        g2e(3,3,-1) = e2*(25.0e0 + e2*(-220.0e0 + e2*( &
        736.916666666667e0 + e2*(-1221.72222222222e0 + 1147.2109375e0*exc**2))))
        g2e(3,3,0) = 1.0e0 + e2*(-12.0e0 + e2*(49.21875e0 + e2*( &
        -83.21875e0 + e2*(68.0408935546875e0 - 31.233955078125e0*e2))))
        g2e(3,3,1) = e2*(1.0e0 + e2*(-2.5e0 + e2*(1.85416666666667e0 &
        + e2*(-0.524305555555555e0 + 0.136197916666666e0*e2))))
        g2e(3,3,2) = e4*(0.015625e0 + e2*(0.00520833333333334e0 + e2 &
        *(0.00490993923611111e0 + 0.00393880208333333e0*e2)))
        g2e(3,3,3) = 0.e0
        g2e(3,3,4) = e8*(6.7816840277778e-6 + 1.35633680555556e-5*e2)
        g2e(3,3,5) = 6.94444444444445e-5*e10
        g2e(3,3,6) = 0.e0
        g2e(3,3,7) = 0.e0
        g2e(3,3,8) = 0.e0
        g2e(3,3,9) = 0.e0
        g2e(3,3,10) = 0.e0
!   --------------------------------------------------------------------
    end subroutine evalge
!=======================================================================
    function sumar(n,arr)
!   --------------------------------------------------------------------
    implicit none
!   --------------------------------------------------------------------
    integer, intent(in) :: n
    integer :: k
    real :: suma,c,t
    real, dimension(n) :: arr
    real :: sumar
!   --------------------------------------------------------------------
        suma = arr(1)
        c = 0.e0
        do k=2,n
        t = suma + arr(k)
        if(abs(suma).ge.abs(arr(k))) then
        c = c + ((suma - t) + arr(k))
        else
        c = c + ((arr(k) - t) + suma)
        end if
        suma = t
        enddo
        sumar = suma + c
!   --------------------------------------------------------------------
    end function sumar
!=======================================================================

    SUBROUTINE hunt(xx,n,x,jlo)
        integer, intent(out) :: jlo
        INTEGER, intent(in) :: n
        REAL, intent(in) :: x,xx(n)
        INTEGER :: inc,jhi,jm
        LOGICAL :: ascnd
        ascnd=xx(n).gt.xx(1)
        if(jlo.le.0.or.jlo.gt.n)then
            jlo=0
            jhi=n+1
            goto 3
        endif
        inc=1
        if(x.ge.xx(jlo).eqv.ascnd)then
1           jhi=jlo+inc
        if(jhi.gt.n)then
            jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
        goto 1
        endif
    else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
            jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
        goto 2
        endif
    endif
3     if(jhi-jlo.eq.1)return
    jm=(jhi+jlo)/2
    if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
    else
        jhi=jm
    endif
    goto 3
    END SUBROUTINE hunt

    ! SUBROUTINE hunt(xx,x,jlo) 
        
    !     IMPLICIT NONE 
        
    !     INTEGER ,INTENT(INOUT) :: jlo 
    !     REAL,INTENT(IN) :: x 
    !     REAL,DIMENSION(:),INTENT(IN) :: xx 
        
    !     !Given an array xx(1:N), and given a value x, returns a value jlo such that x is between xx(jlo) and xx(jlo+1). xx must be monotonic, either increasing or decreasing. jlo=0 or jlo = N is returned to indicate that x is out of range. jlo on input is taken as the initial guess for jlo on output. 
        
    !     INTEGER :: n,inc,jhi,jm 
    !     LOGICAL :: ascnd 

    !     n=size(xx) 
    !     ascnd=(xx(n)>=xx(1))
        
    !     ! True if ascending order of table, false otherwise. 
        
    !     if(jlo<=0.or.jlo>n) then 
    !         ! Input guess not useful. Go immediately to bisection.
    !         jlo=0 
    !         jhi=n+1 
    !     else 
    !         inc = 1 
    !         ! Set the hunting increment. 
    !         if(x>=xx(jlo).eqv.ascnd) then
    !             ! Hunt up: 
    !             do 
    !                 jhi=jlo+inc 
    !                 if(jhi>n) then
    !                     ! Done hunting, since off end of table.
    !                     jhi=n+1 
    !                     exit 
    !                 else 
    !                     if(x<xx(jhi).eqv.ascnd) exit 
    !                     jlo=jhi 
    !                     ! Not done hunting, so double the increment and try again.
    !                     inc=inc+inc 
    !                 endif 
    !             end do 
    !         else ! Huntdown: 
    !             jhi=jlo 
    !             do 
    !                 jlo=jhi-inc 
    !                 if(jlo<1) then !Donehunting,sinceoendoftable. 
    !                     jlo=0 
    !                     exit 
    !                 else 
    !                     if(x>=xx(jlo).eqv.ascnd) exit 
    !                     jhi=jlo ! Not done hunting, so double the increment and try again. 
    !                     inc=inc+inc
    !                 endif 
    !             end do
    !         end if 
    !     end if 
    !     ! Done hunting, value bracketed. 
    !     do ! Hunt is done, so begin the final bisection phase: 
    !         if(jhi-jlo<=1) then 
    !             if(x==xx(n)) jlo=n-1 
    !             if(x==xx(1)) jlo=1 
    !             exit 
    !         else 
    !             jm=(jhi+jlo)/2 
    !             if(x>=xx(jm).eqv.ascnd) then 
    !                 jlo=jm 
    !             else 
    !                 jhi=jm 
    !             end if 
    !         end if 
    !     end do 
        
    ! END SUBROUTINE hunt

    subroutine polint(xp,yp,n,x,y,dy)
    integer :: n,k,m,ns
    real, intent(in) :: x,xp(n),yp(n)
    real, intent(out) :: y,dy
    integer, parameter :: kmax=10
    real :: dlen,dif,dift,ho,hp,w,c(kmax),d(kmax)
    ns=1
    dif=abs(x-xp(1))
    do k=1,n
        dift=abs(x-xp(k))
        if (dift.lt.dif) then
        ns=k
        dif=dift
        endif
        c(k)=yp(k)
        d(k)=yp(k)
    end do
    y=yp(ns)
    ns=ns-1
    do m=1,n-1
        do k=1,n-m
        ho=xp(k)-x
        hp=xp(k+m)-x
        w=c(k+1)-d(k)
        dlen=ho-hp
        if(dlen.eq.0.e0) print *,'failure in polint'
        dlen=w/dlen
        d(k)=hp*dlen
        c(k)=ho*dlen
        end do
        if (2*ns.lt.n-m)then
        dy=c(ns+1)
        else
        dy=d(ns)
        ns=ns-1
        endif
        y=y+dy
    end do
    return
    end subroutine polint

    subroutine set_up_integration(idunit,ur0,Tpot0,dt)

        implicit none

        integer, intent(in) :: idunit 
        real, intent(out) :: ur0,Tpot0,dt
        integer :: demid
        character (len=3) :: hmaxid
        character (len=4) :: chyr,chmon,chday,chhour,chmins,chsecs
        character (len=14) :: exe_ID
        character (len=50) :: fileout
        character (len=10) :: reo,chtide
        character (len=10) :: date
        character (len=8) :: hour
        character (len=100) :: fmt
    
        call read_file_in(idunit,Ur0,Tpot0,dt,demid,hmaxid)

        call timestamp(chyr,chmon,chday,chhour,chmins,chsecs)

        exe_ID = trim(chyr)//trim(chmon)//trim(chday)//trim(chhour)//trim(chmins)//trim(chsecs)

        fileout = trim('K06_'//exe_ID//'.out')

        date = trim(chday)//'/'//trim(chmon)//'/'//trim(chyr)
        hour = trim(chhour)//':'//trim(chmins)//':'//trim(chsecs)

        if (ltide) then 
            chtide = 'yes'
        else 
            chtide = 'no '
        end if

        if (lreo(1)) then
            reo = 'elast'
        else if (lreo(2)) then
            reo = 'visc'
        else if (lreo(3)) then
            reo = 'max'
        else if (lreo(4)) then
            reo = 'bur'
        else if (lreo(5)) then
            reo = 'and'
        else if (lreo(6)) then
            reo = 'qconst'
        else if (lreo(7)) then
            reo = 'dtconst'
        else if (lreo(8)) then
            reo = 'e&l2007'
        end if

        open(unit=idunit,file='../data/exe_DB.dat',status='old',position='append',action='write')

!                   1      2     3       4      5     6      7      8     9     10
        fmt = '(2x,a10,1x,a8,1x,f7.2,3x,f4.2,4x,a3,5x,i1,7x,f5.3,1x,a7,1x,a3,4x,a22)'
!                          1    2     3    4    5            6   7  8    9       10 
        write(idunit,fmt) date,hour,Tpot0,Ur0,trim(chtide),demid,e,reo,hmaxid,trim(fileout)

        close(idunit)

        call create_out_file(idunit,fileout)

    end subroutine set_up_integration

    subroutine read_file_in(idunit,Ur0,Tpot0,dt,demid,hmaxid)

        implicit none 

        integer, intent(in) :: idunit
        real, intent(out) :: Ur0,Tpot0,dt
        integer, intent(out) :: demid
        character (len=3), intent(out) :: hmaxid
        character (len=50) algo
        integer :: tidefl,idreo,k

        open(unit=idunit,file='../in/korenaga06_thermev.in')

        read(idunit,*) algo
        read(idunit,*) Tpot0
        read(idunit,*) algo
        read(idunit,*) Ur0
        read(idunit,*) algo
        read(idunit,*) tidefl
        read(idunit,*) algo
        read(idunit,*) e
        read(idunit,*) algo
        read(idunit,*) demid
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) idreo
        read(idunit,*) algo
        read(idunit,*) lmax,qmax
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) algo
        read(idunit,*) hmaxid
        read(idunit,*) algo
        read(idunit,*) dt

        ltide = tidefl.eq.1

        do k=1,3
            ldem(k) = demid == k
        end do

        do k = 1,8
            lreo(k) = idreo == k
        end do

        lhmax = hmaxid == 'K06'
    
        close(idunit)
        
    end subroutine read_file_in

    subroutine timestamp(chyr,chmon,chday,chhour,chmins,chsecs)
!   --------------------------------------------------------------------
!   This subroutine prints current year, month, day, hour, minutes and 
!   seconds as characters in order to be included as part of the output
!   file.
!   --------------------------------------------------------------------
        implicit none
!       ----------------------------------------------------------------
        integer :: tvals(8),year,month,days,hour,mins,secs
        character (len=4), intent(out) :: chyr,chmon,chday,chhour, &
                                          chmins,chsecs
!       ---------------------------------------------------------------
        call date_and_time(values=tvals)
         year = tvals(1)
        month = tvals(2)
         days = tvals(3)
         hour = tvals(5)
         mins = tvals(6)
         secs = tvals(7)
!       ---------------------------------------------------------------
        write(chyr,'(i4)') year
!       ---------------------------------------------------------------
        if (month.lt.10) then
            write(chmon,'(i1)') month
            chmon = trim('0'//chmon)
        else
            write(chmon,'(i2)') month
        end if
!       ---------------------------------------------------------------
        if (days.lt.10) then
            write(chday,'(i1)') days
            chday = trim('0'//chday)
        else
            write(chday,'(i2)') days
        end if
!       ---------------------------------------------------------------
        if (hour.lt.10) then
            write(chhour,'(i1)') hour
            chhour = trim('0'//chhour)
        else
            write(chhour,'(i2)') hour
        end if
!       ---------------------------------------------------------------
        if (mins.lt.10) then
            write(chmins,'(i1)') mins
            chmins = trim('0'//chmins)
        else
        write(chmins,'(i2)') mins
        end if
!       ---------------------------------------------------------------
        if (secs.lt.10) then
            write(chsecs,'(i1)') secs
            chsecs = trim('0'//chsecs)
        else
            write(chsecs,'(i2)') secs
        end if
!   --------------------------------------------------------------------
    end subroutine timestamp

    subroutine create_out_file(idunit,fileout)
!   -------------------------------------------------------------------
        implicit none
!       ---------------------------------------------------------------
        integer, intent(in) :: idunit
        character (len=50), intent(in) :: fileout
!       ---------------------------------------------------------------
!       This subroutine creates the output file
!       ---------------------------------------------------------------

        open(unit=idunit,file='../out/'//trim(fileout),status='unknown')
!                    1     2     3     4     5     6     7      8
100     format(a1,1x,a4,1x,a7,1x,a6,1x,a6,1x,a9,1x,a4,1x,a5,1x,a14)
!                                1       2         3         4          5        6     7
        write(idunit,100) '#','  t  ',' T_i ',' Q_conv ',' Q_rad ',' Q_tidal ',' u ',' Ur ', &
!                                   8
                             ' Deta_L ratio '
!                              1             2             3             4 
        write(idunit,100) '#',repeat('=',4),repeat('=',7),repeat('=',6),repeat('=',6), &
!                              5             6              7             8
                            repeat('=',9),repeat('=',4),repeat('=',5),repeat('=',14)
!   --------------------------------------------------------------------

    end subroutine create_out_file

    function gammln(xx)
        implicit none
        real :: gammln
        real, intent(in) :: xx
        integer :: j
        real :: ser,stp,tmp,x,y,cof(6)
        save cof,stp
        cof(1) = 76.18009172947146
        cof(2) = -86.50532032941677
        cof(3) = 24.01409824083091
        cof(4) = -1.231739572450155
        cof(5) = .1208650973866179e-2
        cof(6) = -.5395239384953e-5
        stp = 2.5066282746310005
        x=xx
        y=x
        tmp=x+5.5
        tmp=(x+0.5)*log(tmp)-tmp
        ser=1.000000000190015
        do j=1,6
          y=y+1.0
          ser=ser+cof(j)/y
        end do
        gammln=tmp+log(stp*ser/x)
        return
    end function gammln

    FUNCTION rtbis(func,par1,par2,x1,x2,xacc)
    
        implicit none

        INTEGER :: JMAX
        real :: x1,x2,xacc,par1,par2
        REAL :: rtbis,func
        PARAMETER (JMAX=40)
        INTEGER :: j
        REAL :: dx,f,fmid,xmid
        fmid=func(x2,par1,par2)
        f=func(x1,par1,par2)
        if(f*fmid.ge.0.) then 
            print *, 'root must be bracketed in rtbis'
            call abort()
        end if
        if(f.lt.0.)then
            rtbis=x1
            dx=x2-x1
        else
            rtbis=x2
            dx=x1-x2
        endif
        do 11 j=1,JMAX
            dx=dx*.5
            xmid=rtbis+dx
            fmid=func(xmid,par1,par2)
            if(fmid.le.0.)rtbis=xmid
            if(abs(dx).lt.xacc .or. fmid.eq.0.) return
    11    continue
        print *, 'too many bisections in rtbis'
    END

end module korenaga_thermev_subs