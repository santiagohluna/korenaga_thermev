module integrators

    implicit none
    
    contains

    SUBROUTINE rk4(y,dydx,x,h,yout,derivs) 
            
        USE nrtype; USE nrutil, ONLY : assert_eq 
        
        IMPLICIT NONE 
        
        REAL(SP),DIMENSION(:),INTENT(IN) :: y,dydx 
        REAL(SP),INTENT(IN) :: x,h 
        REAL(SP),DIMENSION(:),INTENT(OUT) :: yout
        
        INTERFACE 
            SUBROUTINE derivs(x,y,dydx) 
                USE nrtype 
                IMPLICIT NONE 
                REAL(SP),INTENT(IN) :: x 
                REAL(SP),DIMENSION(:),INTENT(IN) :: y 
                REAL(SP),DIMENSION(:),INTENT(OUT) :: dydx 
            END SUBROUTINE derivs 
        END INTERFACE 
        
        ! Given values for the N variables y and their derivatives dydx known at x,use the fourth order Runge-Kutta method to advance the solution over an interval h and return the incremented variables as yout, which need not be a distinct array from y. y, dydx and yout are all of length N. The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x. 
        
        INTEGER(I4B) :: ndum 
        REAL(SP) :: h6,hh,xh 
        REAL(SP),DIMENSION(size(y)) :: dym,dyt,yt 
        
        ndum=assert_eq(size(y),size(dydx),size(yout),'rk4') 
        hh=h*0.5_sp 
        h6=h/6.0_sp 
        xh=x+hh 
        yt=y+hh*dydx ! First step. 
        
        call derivs(xh,yt,dyt) ! Second step. 
        yt=y+hh*dyt 
    
        call derivs(xh,yt,dym) ! Third step. 
        yt=y+h*dym 
        dym=dyt+dym 
    
        call derivs(x+h,yt,dyt) ! Fourth step. 
        
        yout=y+h6*(dydx+dyt+2.0_sp*dym) ! Accumulate increments with proper weights. 
    
    END SUBROUTINE rk4

    SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,intsch) 
        USE nrtype; USE nrutil, ONLY:nrerror,reallocate 
        USE ode_path 
        
        IMPLICIT NONE 
        REAL(SP),DIMENSION(:),INTENT(INOUT) :: ystart 
        REAL(SP),INTENT(IN) :: x1,x2,eps,h1,hmin 
        INTERFACE 
            SUBROUTINE derivs(x,y,dydx) 
                USE nrtype 
                IMPLICIT NONE 
                REAL(SP),INTENT(IN) :: x 
                REAL(SP),DIMENSION(:),INTENT(IN) :: y 
                REAL(SP),DIMENSION(:),INTENT(OUT) :: dydx 
            END SUBROUTINE derivs 
            SUBROUTINE intsch(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs) 
                USE nrtype 
                IMPLICIT NONE 
                REAL(SP),DIMENSION(:),INTENT(INOUT) :: y
                REAL(SP),DIMENSION(:),INTENT(IN) :: dydx,yscal 
                REAL(SP),INTENT(INOUT) :: x 
                REAL(SP),INTENT(IN) :: htry,eps 
                REAL(SP),INTENT(OUT) :: hdid,hnext 
                INTERFACE 
                    SUBROUTINE derivs(x,y,dydx) 
                        USE nrtype 
                        IMPLICIT NONE 
                        REAL(SP),INTENT(IN) :: x 
                        REAL(SP),DIMENSION(:),INTENT(IN) :: y 
                        REAL(SP),DIMENSION(:),INTENT(OUT) :: dydx 
                    END SUBROUTINE derivs 
                END INTERFACE 
            END SUBROUTINE intsch 
        END INTERFACE 
            
        REAL(SP),PARAMETER :: TINY=1.0e-30_sp
        INTEGER(I4B),PARAMETER :: MAXSTP=10000 
        
        ! Runge-Kutta driver with adaptive step size control. Integrate the array of starting values ystart from x1 to x2 with accuracy eps, storing intermediate results in the module variables in ode_path. h1 should be set as a guessed first step size, hmin as the minimum  allowed step size (can be zero). On output ystart is replaced by values at the end of the integration interval. derivs is the user-supplied subroutine for calculating the right-hand side derivative, while intsch is the name of the stepper routine to be used. 
        
        INTEGER(I4B) :: nstp 
        REAL(SP) :: h,hdid,hnext,x,xsav 
        REAL(SP),DIMENSION(size(ystart)) :: dydx,y,yscal 
        x=x1 
        h=sign(h1,x2-x1) 
        nok=0 
        nbad=0 
        kount=0 
        y(:)=ystart(:) 
        nullify(xp,yp) 
        ! Pointers nullified here, but memory not deallocated. If odeint is called multiple times,calling program should deallocate xp and yp between calls. 
        if(save_steps) then 
            xsav=x-2.0_sp*dxsav 
            allocate(xp(256)) 
            allocate(yp(size(ystart),size(xp))) 
        endif 
        do nstp=1,MAXSTP ! Take at most MAXSTP steps. 
            call derivs(x,y,dydx) 
            yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY  ! Scaling used to monitor accuracy. This general purpose choice can be modified if need be. 
            if(save_steps.and.(abs(x-xsav)>abs(dxsav)))& ! Store intermediate results. 
            call save_a_step 
            if((x+h-x2)*(x+h-x1)>0.0) h=x2-x ! If step size can overshoot, decrease. 
            call intsch(y,dydx,x,h,eps,yscal,hdid,hnext,derivs) 
            if(hdid==h) then 
                nok=nok+1 
            else 
                nbad=nbad+1 
            end if 
            if((x-x2)*(x2-x1)>=0.0) then ! Are we done? 
                ystart(:)=y(:) 
                if(save_steps) call save_a_step ! Save final step. 
                RETURN ! Normal exit. 
            end if 
            if(abs(hnext)<hmin)& 
                call nrerror('step size smaller than minimum in odeint') 
                h=hnext 
        end do 
        call nrerror('too many steps in odeint')

        CONTAINS 
        SUBROUTINE save_a_step 
            kount=kount+1 
            if(kount>size(xp)) then 
                xp=>reallocate(xp,2*size(xp)) 
                yp=>reallocate(yp,size(yp,1),size(xp)) 
            end if 
            xp(kount)=x 
            yp(:,kount)=y(:) 
            xsav=x 
        END SUBROUTINE save_a_step 
    END SUBROUTINE odeint

    SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs) 
        USE nrtype; USE nrutil, ONLY:assert_eq,swap 
        IMPLICIT NONE 
        INTEGER(I4B),INTENT(IN) :: nstep 
        REAL(SP),INTENT(IN) :: xs,htot 
        REAL(SP),DIMENSION(:),INTENT(IN) :: y,dydx 
        REAL(SP),DIMENSION(:),INTENT(OUT) :: yout 
        INTERFACE 
            SUBROUTINE derivs(x,y,dydx) 
                USE nrtype 
                IMPLICIT NONE 
                REAL(SP),INTENT(IN) :: x 
                REAL(SP),DIMENSION(:),INTENT(IN) :: y 
                REAL(SP),DIMENSION(:),INTENT(OUT) :: dydx 
            END SUBROUTINE derivs 
        END INTERFACE 
        ! Modified midpoint step. Dependent variable vector y and its derivative vector dydx are input at xs. Also input is htot, the total step to be taken, and nstep, the number of substeps to be used. The output is returned as yout, which need not be a distinct array from y; if it is distinct, however, then y and dydx are returned undamaged. y, dydx, and yout must all have the same length. 
        INTEGER(I4B) :: n,ndum 
        REAL(SP) :: h,h2,x 
        REAL(SP),DIMENSION(size(y)) :: ym,yn 
        ndum=assert_eq(size(y),size(dydx),size(yout),'mmid') 
        h=htot/nstep ! Step size this trip. 
        ym=y 
        yn=y+h*dydx ! Firststep. 
        x=xs+h 
        call derivs(x,yn,yout) ! Will use yout for temporary storage of derivatives. 
        h2=2.0_sp*h 
        do n=2,nstep ! General step.
            call swap(ym,yn) 
            yn=yn+h2*yout
            x=x+h 
            call derivs(x,yn,yout) 
        end do 
        yout=0.5_sp*(ym+yn+h*yout) ! Last step.
    END SUBROUTINE mmid

    SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs) 
        USE nrtype; USE nrutil, ONLY:arth,assert_eq,cumsum,iminloc,nrerror,& 
        outerdiff,outerprod,upper_triangle 
        USE nr, ONLY:mmid,pzextr 
        IMPLICIT NONE 
        REAL(SP),DIMENSION(:),INTENT(INOUT) :: y 
        REAL(SP),DIMENSION(:),INTENT(IN) :: dydx,yscal 
        REAL(SP),INTENT(INOUT) :: x 
        REAL(SP),INTENT(IN) :: htry,eps 
        REAL(SP),INTENT(OUT) :: hdid,hnext 
        INTERFACE 
            SUBROUTINE derivs(x,y,dydx) 
                USE nrtype 
                IMPLICIT NONE 
                REAL(SP),INTENT(IN) :: x 
                REAL(SP),DIMENSION(:),INTENT(IN) :: y 
                REAL(SP),DIMENSION(:),INTENT(OUT) :: dydx 
            END SUBROUTINE derivs 
        END INTERFACE 
        INTEGER(I4B),PARAMETER :: IMAX=9,KMAXX=IMAX-1 
        REAL(SP),PARAMETER :: SAFE1=0.25_sp,SAFE2=0.7_sp,REDMAX=1.0e-5_sp,& 
        REDMIN=0.7_sp,TINY=1.0e-30_sp,SCALMX=0.1_sp 
        ! Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust stepsize. Input are the dependent variable vector y and its derivative dydx at the starting  value of the independent variable x. Also input are the step size to be attempted htry, the required accuracy eps, and the vector yscal against which the error is scaled. On output, y and x are replaced by their new values, hdid is the step size that was actually accomplished, and hnext is the estimated next step size. derivs is the user-supplied subroutine that computes the right-hand-side derivatives. y, dydx, and yscal must all have the same length. Be sure to set htry on successive steps to the value of hnext returned from the previous step, as is the case if the routine is called by odeint. Parameters: KMAXX is the maximum row number used in the extrapolation; IMAX is the next rownumber; SAFE1 and SAFE2 are safety factors; REDMAX is the maximum factor used when a step size is reduced, REDMIN the minimum; TINY prevents division by zero; 1/SCALMX is the maximum factor by which a step size can be increased. 
        INTEGER(I4B) :: k,km,ndum 
        INTEGER(I4B),DIMENSION(IMAX) :: nseq=(/2,4,6,8,10,12,14,16,18/) 
        INTEGER(I4B),SAVE :: kopt,kmax 
        REAL(SP),DIMENSION(KMAXX,KMAXX),SAVE :: alf 
        REAL(SP),DIMENSION(KMAXX) :: err 
        REAL(SP),DIMENSION(IMAX),SAVE :: a 
        REAL(SP),SAVE :: epsold=-1.0_sp,xnew 
        REAL(SP) :: eps1,errmax,fact,h,red,scale,wrkmin,xest 
        REAL(SP),DIMENSION(size(y)) :: yerr,ysav,yseq 
        LOGICAL(LGT) :: reduct 
        LOGICAL(LGT),SAVE :: first=.true. 
        ndum=assert_eq(size(y),size(dydx),size(yscal),'bsstep') 
        if(eps/=epsold) then ! A new tolerance, so reinitialize. 
            hnext=-1.0e29_sp !"Impossible" values. 
            xnew=-1.0e29_sp 
            eps1=SAFE1*eps 
            a(:)=cumsum(nseq,1) 
            ! Compute α(k,q): 
            where(upper_triangle(KMAXX,KMAXX)) alf=eps1**& 
            (outerdiff(a(2:),a(2:))/outerprod(arth(&
            3.0_sp,2.0_sp,KMAXX),(a(2:)-a(1)+1.0_sp))) 
            epsold=eps 
            do kopt=2,KMAXX-1 ! Determine optimal row number for convergence. 
                if(a(kopt+1)>a(kopt)*alf(kopt-1,kopt)) exit 
            end do 
            kmax=kopt 
        end if 
        h=htry 
        ysav(:)=y(:) ! Save the starting values. 
        if(h/=hnext.or.x/=xnew) then ! A new step size or a new integration: Reestablish the order window. 
            first=.true. 
            kopt=kmax 
        end if 
        reduct=.false. 
        main_loop: do 
            do k=1,kmax ! Evaluate the sequence of modified midpoint integrations. 
                xnew=x+h 
                if(xnew==x) call nrerror('stepsizeunderflowinbsstep') 
                call mmid(ysav,dydx,x,h,nseq(k),yseq,derivs) 
                xest=(h/nseq(k))**2 ! Squared, since error series is even. 
                call pzextr(k,xest,yseq,y,yerr) ! Perform extrapolation. 
                if(k/=1) then ! Compute normalized error estimate (k). 
                    errmax=maxval(abs(yerr(:)/yscal(:))) 
                    errmax=max(TINY,errmax)/eps ! Scale error relative to tolerance. 
                    km=k-1 
                    err(km)=(errmax/SAFE1)**(1.0_sp/(2*km+1)) 
                end if 
                if(k/=1.and.(k>=kopt-1.or.first)) then ! In order window. 
                    if(errmax<1.0) exit main_loop ! Converged. 
                    if(k==kmax.or.k==kopt+1) then ! Check for possible step size reduction. 
                        red=SAFE2/err(km) 
                        exit 
                    else if(k==kopt) then 
                        if(alf(kopt-1,kopt)<err(km))then 
                            red=1.0_sp/err(km) 
                            exit 
                        end if 
                    else if(kopt==kmax) then 
                        if(alf(km,kmax-1)<err(km)) then 
                            red=alf(km,kmax-1)*SAFE2/err(km) 
                            exit 
                        end if 
                    else if(alf(km,kopt)<err(km)) then 
                        red=alf(km,kopt-1)/err(km) 
                        exit 
                    end if 
                end if 
            end do 
            red=max(min(red,REDMIN),REDMAX) ! Reduce step size by at least REDMIN and atmost REDMAX. 
            h=h*red 
            reduct=.true. 
        end do main_loop ! Try again. 
        x=xnew ! Successful step taken. 
        hdid=h 
        first=.false. 
        kopt=1+iminloc(a(2:km+1)*max(err(1:km),SCALMX)) ! Compute optimal row for convergence and corresponding stepsize. 
        scale=max(err(kopt-1),SCALMX) 
        wrkmin=scale*a(kopt) 
        hnext=h/scale 
        if(kopt>=k.and.kopt/=kmax.and..not.reduct) then ! Check for possible order increase, but not if step size was just reduced. 
            fact=max(scale/alf(kopt-1,kopt),SCALMX) 
            if(a(kopt+1)*fact<=wrkmin) then 
                hnext=h/fact
                kopt=kopt+1 
            end if 
        end if 
    END SUBROUTINE bsstep

    SUBROUTINE pzextr(iest,xest,yest,yz,dy) 
        USE nrtype; USE nrutil, ONLY:assert_eq,nrerror 
        IMPLICIT NONE 
        INTEGER(I4B),INTENT(IN) :: iest 
        REAL(SP),INTENT(IN) :: xest 
        REAL(SP),DIMENSION(:),INTENT(IN) :: yest 
        REAL(SP),DIMENSION(:),INTENT(OUT) :: yz,dy 
        ! Use polynomial extrapolation to evaluate N functions at x=0 by fitting a polynomial to a sequence of estimates with progressively smaller values x=xest, and corresponding  function vectors yest. This call is number iest in the sequence of calls. Extrapolated function values are outputas yz, and their estimated error is output as dy. yest, yz, and dy are arrays of length N. 
        INTEGER(I4B),PARAMETER :: IEST_MAX=16 
        INTEGER(I4B) :: j,nv 
        INTEGER(I4B),SAVE::nvold=-1 
        REAL(SP) :: delta,f1,f2 
        REAL(SP),DIMENSION(size(yz)) :: d,tmp,q 
        REAL(SP),DIMENSION(IEST_MAX),SAVE::x 
        REAL(SP),DIMENSION(:,:),ALLOCATABLE,SAVE :: qcol 
        nv=assert_eq(size(yz),size(yest),size(dy),'pzextr') 
        if(iest>IEST_MAX) &
        call nrerror('pzextr: probable misuse, too much extrapolation') 
        if(nv/=nvold) then ! Set up internal storage. 
            if(allocated(qcol)) deallocate(qcol) 
            allocate(qcol(nv,IEST_MAX)) 
            nvold=nv 
        end if 
        x(iest)=xest ! Save current independent variable. 
        dy(:)=yest(:) 
        yz(:)=yest(:) 
        if(iest==1) then ! Store first estimate in first column. 
            qcol(:,1)=yest(:) 
        else 
            d(:)=yest(:) 
            do j=1,iest-1 
                delta=1.0_sp/(x(iest-j)-xest) 
                f1=xest*delta 
                f2=x(iest-j)*delta 
                q(:)=qcol(:,j) ! Propagate tableau 1 diagonal more.
                qcol(:,j)=dy(:) 
                tmp(:)=d(:)-q(:) 
                dy(:)=f1*tmp(:) 
                d(:)=f2*tmp(:) 
                yz(:)=yz(:)+dy(:) 
            end do 
            qcol(:,iest)=dy(:) 
        end if 
    END SUBROUTINE pzextr
    
end module integrators