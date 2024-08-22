MODULE ode_path 
    ! On output nok and nbad are the number of good and bad (but retried and fixed) steps taken. If save steps is set to true in the calling program, then intermediate values are stored in xp and yp at intervals greater than dxsav. kount is the total number of saveds teps. 
    USE nrtype 
    INTEGER(I4B) :: nok,nbad,kount 
    LOGICAL(LGT),SAVE :: save_steps=.false. 
    REAL(SP) :: dxsav 
    REAL(SP),DIMENSION(:),POINTER :: xp 
    REAL(SP),DIMENSION(:,:),POINTER :: yp 
END MODULE ode_path