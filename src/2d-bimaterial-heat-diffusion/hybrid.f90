PROGRAM fde_2D_heattransfer
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	
	!Variable declaration
	!-------------------------------------------------------------------------------------------------------------------------------------- 
	REAL(KIND = 8), PARAMETER		          :: PI = 4.*ATAN(1.0d0)						! pi to double precision
	INTEGER, PARAMETER                        :: N = 200                  				    ! grid size
	INTEGER, PARAMETER                        :: M = 200                  				    ! grid size
	INTEGER, PARAMETER                        :: Z = 100                 				    ! grid size
	REAL(KIND = 8), PARAMETER		          :: alpha1 = 0.001_8          				    ! m^2/s) thermal diffusivity of Material(M1)
	REAL(KIND = 8), PARAMETER		          :: alpha2 = 0.007_8          				    ! m^2/s) thermal diffusivity of Material(M2)
	REAL(KIND = 8), PARAMETER 		          :: xl = 0.0_8								    ! lower limit of x  
	REAL(KIND = 8), PARAMETER 		          :: xh = 4.0_8  							    ! metres) upper limit of x
	REAL(KIND = 8), PARAMETER 		          :: yl = 0.0_8								    ! lower limit of y  
	REAL(KIND = 8), PARAMETER 		          :: yh = 4.0_8  								! metres) upper limit of y
	REAL(KIND = 8)					          :: dx,Lx,dy,Ly,dt,t,dx2,dy2,dt_dx2,dt_dy2     ! Constant values
	REAL(KIND = 8), PARAMETER		          :: tmax = 2.0_8                 			    ! (seconds) End time
		
	INTEGER(KIND = 4) 				          :: pid		                                ! Proc ID
	INTEGER(KIND = 4) 				          :: np	    								    ! Num Proc
	INTEGER(KIND = 4) 				          :: ierr		    							! Error status  
	REAL(KIND = 8) 					          :: walltime  								    ! Wall time	
	INTEGER, PARAMETER                        :: root = 0                                   ! Root processor ID
									          
	!Variables local to processors
	REAL(KIND = 8)					          :: lambdax,lambday                		    ! Variables for Stability Criteria
	REAL(KIND = 8), ALLOCATABLE		          :: x(:)                   					! array of grid location values in x direction
	REAL(KIND = 8), ALLOCATABLE		          :: y(:)                   					! array of grid location values in y direction
	REAL(KIND = 8), ALLOCATABLE	              :: f(:),Tpf(:,:),Tp0(:,:)                     ! Matrix of Temperature
	REAL(KIND = 8), ALLOCATABLE	              :: Ap(:,:)               	                    ! Matrix for constant values
	INTEGER                                   :: ml,mh,kl,kh                                ! Partition 	
	
	REAL(KIND = 8)						  	  :: xval(1:N+1)                                ! array of grid location values in x direction
	REAL(KIND = 8)						  	  :: yval(1:M+1)                                ! array of grid location values in y direction
	REAL(KIND = 8)						  	  :: fval(1:M+1)                                ! array of grid location values in y direction
	REAL(KIND = 8)                    	  	  :: Tf(1:M+1,1:N+1)                            ! Matrix of Temperature	
	INTEGER, PARAMETER                    	  :: iunit = 30            	                    ! Variable to save data in file
	INTEGER, PARAMETER                        :: iwrite = 50                 				! Variable to save data in file
	INTEGER                                   :: iter             				            ! Variable to save data in file
	CHARACTER(80)                         	  :: filename                                   ! variable for file name
	
	
	! Initialise MPI   
	CALL MPI_INIT(ierr)                                                         			! Initialize MPI Program
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr) 							    			    ! Get Size of Communicator (Number of procs)s
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)    							    			! Get Rank of Communicator (proc id)
	
	!--------------------------------------------------------------------------------------------------------------------------------------
		
	! Set all the values for the domain Variables
	CALL GetVars(pid,root,N,M,xl,xh,yl,yh,dx,Lx,dy,Ly,Z,tmax,dt,t,dx2,dy2,dt_dx2,dt_dy2)
	
	! Partition domain
	CALL Partition(pid,np,N,M,Lx,Ly,dx,dy,ml,mh,kl,kh)
	
	ALLOCATE(x(kl:kh))
	ALLOCATE(y(1:M+1))
	ALLOCATE(Ap(1:M+1,kl:kh))
	
	
	! Ensure if stability criteria is met
	CALL Stability_Criteria(pid,np,root,M,Lx,dt,dx2,dy2,ml,mh,kl,kh,x,Ap,alpha1,alpha2,lambdax,lambday)	
	
	! Store Wall Time (Real World Time)
	IF (pid == root) THEN
		walltime = MPI_Wtime( )
  	END IF
	
	
	ALLOCATE(f(kl:kh))	
	ALLOCATE(Tp0(1:M+1,kl:kh))
	ALLOCATE(Tpf(1:M+1,kl:kh))	
		
	! Initialise Temperature values
	CALL Initialise(pid,np,pi,N,M,ml,mh,kl,kh,dx,dy,x,y,Lx,Ly,f,Tp0,Tpf,Tf)	
	
	! Gather data
	CALL MPI_GATHER(x(ml:mh),mh-ml+1,MPI_DOUBLE_PRECISION,xval(2),mh-ml+1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
	
	CALL MPI_GATHER(f(ml:mh),mh-ml+1,MPI_DOUBLE_PRECISION,fval(2),mh-ml+1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
	
	CALL MPI_GATHER(Tpf(1:M+1,ml:mh),(M+1)*(mh-ml+1),MPI_DOUBLE_PRECISION,Tf(1,2),&
	(M+1)*(mh-ml+1),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
		
	! Write Initial Temperature Data
	IF (pid.EQ.0) THEN
		! Write Temperature Data
		iter =1
		xval(1) = 0
		xval(N+1) = Lx
		yval = y
		CALL Write_PlotFile(pid,root,N,M,t,xval,yval,Tf,filename,iunit,iter)
	END IF
	
	t = t+dt ! Initial stepping through time
	
	DO WHILE (t.LE.tmax)
	
		CALL HeatEqn(M,Tp0,ml,mh,kl,kh,dt_dx2,dt_dy2,Ap,Tpf)
		
		CALL MPIswap(pid,np,M,N,ml,mh,kl,kh,Tpf)
		
		! Gather data
		CALL MPI_GATHER(Tpf(1:M+1,ml:mh),(M+1)*(mh-ml+1),MPI_DOUBLE_PRECISION,Tf(1,2),(M+1)*(mh-ml+1),&
		MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)	
		
		! Write Temperature Data
		IF (pid.EQ.0) THEN
			iter = iter +1
			IF(MOD(iter,iwrite).eq.0) THEN			
				CALL Write_PlotFile(pid,root,N,M,t,xval,yval,Tf,filename,iunit,iter)
			END IF
		END IF
		Tp0 = Tpf
		t = t + dt
	END DO
	
	!! Run time
	!IF (pid == root) THEN
	!	walltime = MPI_Wtime ( ) - walltime
	!	WRITE ( *, '(a)' ) ' '
	!	WRITE ( *, '(a,f12.5)' ) '  Run time seconds : ', walltime
	!END IF
	
	CALL MPI_FINALIZE ( ierr )
	
END PROGRAM fde_2D_heattransfer

!-------------------------------------------------------------------------------------------------
SUBROUTINE GetVars( pid,root,N,M,xl,xh,yl,yh,dx,Lx,dy,Ly,Z,tmax,dt,t,dx2,dy2,dt_dx2,dt_dy2)

	IMPLICIT NONE
	INCLUDE 'mpif.h'
	INTEGER, INTENT(IN)                       :: pid,root,N,M,Z
	REAL(KIND = 8),	   INTENT(IN)             :: xl,xh,yl,yh,tmax
	REAL(KIND = 8),	   INTENT(OUT)            :: dx,Lx,dy,Ly,dt,t,dx2,dy2,dt_dx2,dt_dy2
	REAL(KIND = 8)                            :: bcast_data(1:10)
	INTEGER(KIND = 4) 				          :: ierr
	
	!$omp parallel if(t<1)
	IF (pid.EQ.root) then
		Lx = xh-xl 					! Domain length in x direction
		dx = Lx/N					! Grid spacing in x direction
		Ly = yh-yl 					! Domain length in y direction
		dy = Ly/M					! Grid spacing in y direction
		dt = tmax/Z                 ! calculation of time interval
		t = 0                       ! initial value for time	
		dx2 = dx**2
		dy2 = dy**2
		dt_dx2 = dt/(dx**2)
		dt_dy2 = dt/(dy**2)	
		bcast_data = [Lx,Ly,dx,dy,t,dt,dx2,dy2,dt_dx2,dt_dy2]		
	END IF
	!$omp end parallel
	
	! Broadcast the domain and time grid size values from root to all processors
	CALL MPI_BCAST(bcast_data,10,MPI_REAL8,root,MPI_COMM_WORLD,ierr)
	
	!Unpack the broadcasted data
	Lx 	   = bcast_data(1)
	Ly 	   = bcast_data(2)
	dx 	   = bcast_data(3)
	dy 	   = bcast_data(4)
	t  	   = bcast_data(5)
	dt 	   = bcast_data(6)
	dx2    = bcast_data(7)
	dy2    = bcast_data(8)
	dt_dx2 = bcast_data(9)
	dt_dy2 = bcast_data(10)	

END SUBROUTINE GetVars 
!-------------------------------------------------------------------------------------------------
! Divide the domain between number of processors
SUBROUTINE Partition(pid,np,N,M,Lx,Ly,dx,dy,ml,mh,kl,kh)

	IMPLICIT NONE
	INTEGER(KIND = 4),INTENT(IN)              ::  pid,np,N,M 
	REAL(KIND = 8), INTENT(IN)                ::  Lx,Ly,dx,dy	
	INTEGER(KIND = 4),INTENT(OUT)             ::  ml,mh,kl,kh		
	INTEGER  nlocal,overld,i
	
	! Calculate Grids in partition
	!$omp end parallel
	IF (MOD((N+1)-2,np).ne.0) THEN
		overld = MOD((N+1)-2,np)                 
		nlocal = ((N+1)-2-overld)/np			
	ELSE
		nlocal  = (N-2)/np		
	END IF  
	!$omp end parallel
	
	ml      = pid*nlocal + 1 + 1  ! First point in partition
	mh      = ml + nlocal - 1     ! Last point in partition	
	
	IF (pid.EQ.(np-1)) then
		kl = N
		kh = N+1		
	ELSE 
		kl = ml-1 
	    kh = mh+1 
	END IF
	
END SUBROUTINE Partition
!!-------------------------------------------------------------------------------------------------	
SUBROUTINE Stability_Criteria(pid,np,root,M,Lx,dt,dx2,dy2,ml,mh,kl,kh,x,Ap,alpha1,alpha2,lambdax,lambday)

	IMPLICIT NONE
	INCLUDE 'mpif.h'
	INTEGER(KIND = 4),INTENT(IN)              ::  root,pid,np,M 
	REAL(KIND = 8), INTENT(IN)                ::  Lx,dt,dx2,dy2,alpha1,alpha2
	INTEGER(KIND = 4),INTENT(IN)              ::  ml,mh,kl,kh
	REAL(KIND = 8),INTENT(INOUT)              ::  Ap(1:M+1,kl:kh)	
	REAL(KIND = 8),INTENT(OUT)                ::  lambdax,lambday	
	REAL(KIND = 8),INTENT(IN)                 ::  x(kl:kh)	
	INTEGER                                   ::  i,j,ierr
	REAL(KIND = 8)                            ::  stable,allones		
	REAL(KIND = 8)                            ::  stability(1:np+1)
	
	stable = 1
	allones = 1
	
	! Setting diffusion value based location of grid point in the domain
	!$omp parallel if(Lx>0)
	!$omp parallel do default(none) shared(Lx,kl,kh,alpha1,alpha2,M)
	!$omp private(i,j,Ap,x,y,f)
	DO j= kl,kh
		DO i=1,M+1
		    IF (x(j).LT.Lx/2.0_8) THEN	
				Ap(i,j)= alpha1
			ELSE IF (x(j).GT.Lx/2.0_8) THEN	
				Ap(i,j)= alpha2
			ELSE
				Ap(i,j) = (alpha1 +alpha2)/2	
			END IF
		END DO
	END DO
	!$omp end parallel do
	!$omp end parallel
	
	!$omp parallel if(lambdax>0.and.lambday>0)
	!$omp parallel do default (none) shared (lambdax, lambday,dt,dx2,dy2,Stability) &
	!$omp private (i,j,Ap)
	! Calculating stability criteria (lambda) at each grid point
	DO j=kl,kh
		DO i=1,M+1
			lambdax = (Ap(i,j) * dt)/(dx2)
			lambday = (Ap(i,j) * dt)/(dy2)
			IF ((lambdax.GE.0.5).AND.(lambday.GE.0.5)) THEN
				stable = 0
			ENDIF
		END DO	
	END DO
	!$omp end parallel do
	!$omp end parallel
	! Gather stability values to root processor
	CALL MPI_GATHER(stable,1,MPI_DOUBLE_PRECISION,stability,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
	
	! Verify if stability criteria is met
	IF (pid.EQ.root) THEN
		DO i = 1,np
			IF (stability(i).NE.1) THEN
				allones = 0			
	        END IF
		END DO
		IF (allones.EQ.0) THEN
			WRITE(*,*) "The program is going to exit since stability criteria is not met"
			CALL EXIT(0)
		ELSE
			WRITE(*,*) "Stability criteria is met"
		END IF
	END IF	       
END SUBROUTINE Stability_Criteria
!!-------------------------------------------------------------------------------------------------	
SUBROUTINE Initialise(pid,np,pi,N,M,ml,mh,kl,kh,dx,dy,x,y,Lx,Ly,f,Tp0,Tpf,Tf)         

	IMPLICIT none	
	INTEGER, INTENT(IN)                       :: pid,np
	INTEGER, INTENT(IN)                       :: N,M,ml,mh,kl,kh
	REAL(KIND = 8), INTENT(IN)                :: pi
	REAL(KIND = 8), INTENT(INOUT)             :: x(kl:kh),y(1:M+1)
	REAL(KIND = 8), INTENT(IN)                :: Lx,Ly,dx,dy
	REAL(KIND = 8), INTENT(OUT) 			  :: Tp0(1:M+1,kl:kh),Tpf(1:M+1,kl:kh),Tf(1:M+1,1:N+1)	  
	REAL(KIND = 8), INTENT(INOUT) 			  :: f(kl:kh)
	INTEGER ::  i,j			
	
	! Calculate grid spacing in x-direction
	
	DO i = kl,kh
		x(i) = (i-1)*dx
	END DO
	! Calculate grid spacing in y-direction
	DO i = 1,M+1
		y(i) = (i-1)*dy
	END DO
	
	! Calculate initial temperatures at boundaries
	!$omp parallel do default(none) shared(Lx,Ly,pi,kl,kh)
	!$omp private(i,x,y,f)
	DO i = kl,kh
		f(i) = 	100*(x(i)/Lx) + ((sin(pi * (y(i)/Ly)))**2)
	END DO	
	!$omp end parallel do
	
	!$omp parallel do default(none) shared(kl,kh,M)
	!$omp private(i,j,f,Tp0,Tpf,f)
	DO j=kl,kh			
		DO i=1,M+1							
			IF ((i.EQ.0).OR.(i.EQ.M+1)) THEN
				! Top and Bottom boundary conditions
				Tp0(i,j) = f(j)
				Tf(i,j) = f(j)
			ELSE IF (j.EQ.N+1) THEN
				! Right boundary condition
				Tp0(i,j) = 100
				Tf(i,j) = 100
			ELSE
				! Left boundary condition
				Tp0(i,j) = 0
				Tf(i,j) = 0
			END IF
			Tpf(i,j) = Tp0(i,j)
		END DO
	END DO	
	!$omp end parallel do
END SUBROUTINE Initialise	
!-------------------------------------------------------------------------------------------------
SUBROUTINE HeatEqn(M,Tp0,ml,mh,kl,kh,dt_dx2,dt_dy2,Ap,Tpf)  
	IMPLICIT NONE
	INTEGER :: M,ml,mh,kl,kh
	REAL(KIND = 8), INTENT(IN)                :: dt_dx2,dt_dy2
	REAL(KIND = 8), INTENT(IN)                :: Ap(1:M+1,kl:kh)
	REAL(KIND = 8), INTENT(INOUT)             :: Tp0(1:M+1,kl:kh)
	REAL(KIND = 8), INTENT(INOUT)             :: Tpf(1:M+1,kl:kh)	
	INTEGER                                   :: i,j
	
	!$omp parallel do default (none) shared (Tpf,dt_dy2,dt_dx2) &
	!$omp private (i,j,Ap,Tp0)
	DO j = ml, mh		
		DO i= 2, M					
			Tpf(i,j) = ((Ap(i,j+1)+ Ap(i,j)/2)*dt_dx2*(Tp0(i,j+1)-Tp0(i,j)))-&
					   ((Ap(i,j)+ Ap(i,j-1)/2)*dt_dx2*(Tp0(i,j)-Tp0(i,j-1)))+&
					   (dt_dy2*((Ap(i+1,j)*Tp0(i+1,j))-(2*Ap(i,j)*Tp0(i,j))+(Ap(i-1,j)*Tp0(i-1,j))))+ Tp0(i,j)				
		END DO					
	END DO		
	!$omp end parallel do
	Tp0(2:M,ml:mh) = Tpf(2:M,ml:mh)	
	
END SUBROUTINE HeatEqn
!-------------------------------------------------------------------------------------------------
SUBROUTINE MPIswap(pid,np,M,N,ml,mh,kl,kh,Tpf)
 
	IMPLICIT none
	INCLUDE 'mpif.h'  
	INTEGER(KIND = 4), INTENT(IN)             :: pid,np
	INTEGER(KIND = 4), INTENT(IN)             :: M,N,ml,mh,kl,kh
	REAL(KIND = 8),    INTENT(INOUT)          :: Tpf(1:M+1,kl:kh)  
	INTEGER, ALLOCATABLE                      :: reqs(:), stats(:,:)
	INTEGER, PARAMETER 						  :: TagL = 1            ! Tag Send LEFT
	INTEGER, PARAMETER 						  :: TagR = 2            ! Tag Send Right 
	INTEGER                                   :: i,j,q,ierr

	IF (pid.EQ.0) THEN
		ALLOCATE(reqs(2))
		ALLOCATE(stats(MPI_STATUS_SIZE,2))
    ELSEIF (pid.EQ.(np-1))	THEN
		ALLOCATE(reqs(2))
		ALLOCATE(stats(MPI_STATUS_SIZE,2))
	ELSE
		ALLOCATE(reqs(4))
		ALLOCATE(stats(MPI_STATUS_SIZE,4))
	END IF
		
	q = 0
	
	IF ( pid .GT. 0) THEN 		! Send to left		
		q = q + 1

		CALL MPI_ISEND ( Tpf(:,ml), M+1, MPI_DOUBLE_PRECISION, pid-1, TagL, &
		MPI_COMM_WORLD, reqs(q), ierr )		
	END IF
	
	IF ( pid .LT. np - 1 ) THEN	! Receive from right 
		q = q+1
				
		CALL MPI_IRECV ( Tpf(:,kh), M+1,  MPI_DOUBLE_PRECISION, pid+1, TagL, &
		MPI_COMM_WORLD, reqs(q), ierr )		
	END IF
		
	IF ( pid .LT. np - 1 ) THEN	! Send to right
		q = q + 1
		
		CALL MPI_ISEND ( Tpf(:,mh), M+1, MPI_DOUBLE_PRECISION, pid+1, TagR, &
		MPI_COMM_WORLD,reqs(q), ierr)		
	END IF	
	
	IF ( pid .GT.0 ) THEN		! Receive from left
		q = q + 1
		
		CALL MPI_IRECV ( Tpf(:,kl), M+1, MPI_DOUBLE_PRECISION, pid-1, TagR, &
		MPI_COMM_WORLD,reqs(q), ierr )		
	END IF
	
	! Waiting for all the requests to complete
	CALL MPI_WAITALL ( q , reqs , stats , ierr )
	DEALLOCATE(reqs)
	DEALLOCATE(stats)
	
END SUBROUTINE MPIswap 
!-------------------------------------------------------------------------------------------------
SUBROUTINE Write_PlotFile(pid,root,N,M,t,xval,yval,Tf,filename,iunit,iter)
 
	IMPLICIT none	  
	INTEGER(KIND = 4), INTENT(IN)             :: pid,root,iunit
	INTEGER(KIND = 4), INTENT(IN)             :: iter
	INTEGER(KIND = 4), INTENT(IN)             :: N,M
	REAL(KIND = 8),    INTENT(IN)             :: t
	REAL(KIND = 8),    INTENT(INOUT)          :: xval(1:N+1),yval(1:M+1)
	REAL(KIND = 8),    INTENT(INOUT)          :: Tf(1:M+1,1:N+1)
	CHARACTER(80),     INTENT(INOUT)          :: filename  	
	INTEGER                                   :: i,j,ierr		
		
	! Set up a file to write Temperature values
	WRITE(filename,'(a,i4.4,a)') 'TecPlot2D_Temperature_',iter,'.tec'	
	OPEN ( UNIT = iunit, FILE = filename, FORM = 'formatted', ACCESS = 'sequential', STATUS = 'replace', IOSTAT = ierr )	
	IF ( ierr /= 0 ) THEN
		WRITE ( *, '(a)' ) '  Error opening file : tecplot_2D '
		STOP
	END IF
	WRITE ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
	WRITE ( iunit, '(a)' ) 'Variables=' // trim ( '"Time","X","Y","T"' )
	WRITE ( iunit, '(a,i6,a,i6,a,a)' ) 'Zone I=', M+1, ', J=', N+1, ', F=POINT'
	DO j=1,N+1;DO i=1,M+1; WRITE(iunit,'(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6)') t, xval(j), yval(i), Tf(i,j);END DO;END DO	
	
END SUBROUTINE Write_PlotFile 