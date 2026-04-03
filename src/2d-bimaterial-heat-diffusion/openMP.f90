PROGRAM fde_2D_heattransfer
	IMPLICIT NONE
	
	!Variable declaration
	REAL(KIND = 8), PARAMETER		:: PI = 4.*ATAN(1.0d0)						! pi to double precision
	INTEGER, PARAMETER              :: N = 200                  				! grid size
	INTEGER, PARAMETER              :: M = 200                  				! grid size
	INTEGER, PARAMETER              :: Z = 100                 					! grid size
	REAL(KIND = 8), PARAMETER		:: alpha1 = 0.001_8          				! m^2/s) thermal diffusivity of Material(M1)
	REAL(KIND = 8), PARAMETER		:: alpha2 = 0.007_8          				! m^2/s) thermal diffusivity of Material(M2)
	REAL(KIND = 8), PARAMETER 		:: xl = 0.0_8								! lower limit of x  
	REAL(KIND = 8), PARAMETER 		:: xh = 4.0_8  								! metres) upper limit of x
	REAL(KIND = 8), PARAMETER 		:: yl = 0.0_8								! lower limit of y  
	REAL(KIND = 8), PARAMETER 		:: yh = 4.0_8  								! metres) upper limit of y
	REAL(KIND = 8)					:: dx,Lx,dy,Ly,dt,t,dx2,dy2,dt_dx2,dt_dy2   ! Constant values
	REAL(KIND = 8)					:: lambdax,lambday   						! Variables for Stability Criteria		
	REAL(KIND = 8), PARAMETER		:: tmax = 2.0_8                 			! (seconds) End time
	REAL(KIND = 8)					:: x(1:N+1)                   				! array of grid location values in x direction
	REAL(KIND = 8)					:: y(1:M+1)                   				! array of grid location values in y direction
	REAL(KIND=8)                    :: Tf(1:M+1,1:N+1),T0(1:M+1,1:N+1)          ! Matrix of Temperature
	REAL(KIND=8)                    :: f(1:N+1),A(1:M+1,1:N+1)               	! Matrix for constant values
	INTEGER                         :: i,j,stability            				! Variables for iteration and stability
	INTEGER                         :: iunit,ierr,iter             				! Variable to save data in file
	INTEGER, PARAMETER              :: iwrite = 50                 				! Variable to save data in file
	CHARACTER(80)                   :: filename                 				! variable for file name
	REAL(KIND = 8)					:: time1,time2              				! Variables to measure run-time
	!-------------------------------------------------------------------------------------------------------- 
	Lx = xh-xl 					! Domain length in x direction
	dx = Lx/N					! Grid spacing in x direction
	Ly = yh-yl 					! Domain length in y direction
	dy = Ly/M					! Grid spacing in y direction
	dt = tmax/Z                 ! calculation of time interval
	t = 0                       ! initial value for time	
	iter = 1
		
	! Calculating grid locations in x and y directions
	DO j=1,N+1
		x(j) = (j-1)*dx	        	
	END DO
	
	DO i=1,N+1
		y(i) = (i-1)*dy 
	END DO
	
	DO j=1,N+1
		f(j) = 	100*(x(j)/Lx) + ((sin(pi * (y(j)/Ly)))**2)
	END DO
	
	! Calculate the constant values	
	!$omp parallel do default (none) shared (A,x,Lx) &
	!$omp private (i,j)	
	DO j= 1,N+1
		DO i=1,M+1
		    IF (x(j).LT.Lx/2.0_8) THEN	
				A(i,j)= alpha1
			ELSE IF (x(j).GT.Lx/2.0_8) THEN	
				A(i,j)= alpha2
			ELSE
				A(i,j) = (alpha1 +alpha2)/2	
			END IF
		END DO
	END DO
	!$omp end parallel do
	dx2 = dx**2
	dy2 = dy**2
	dt_dx2 = dt/(dx**2)
	dt_dy2 = dt/(dy**2)	
	!-------------------------------------------------------------------------------------------------------- 	
	! Ensure if stability criteria is met
	!$omp parallel do default (none) shared (lambdax, lambday,dt,dx2,dy2,Stability) &
	!$omp private (i,j,A)
	DO j=1,N+1
		DO i=1,N+1
			lambdax = (A(i,j) * dt)/(dx2)
			lambday = (A(i,j) * dt)/(dy2)
			IF ((lambdax.LT.0.5).AND.(lambday.LT.0.5)) THEN
				stability=1
			ELSE				
				stability=0
			ENDIF
		END DO	
	END DO
	!$omp end parallel do
	
	IF (stability.EQ.1) THEN
		WRITE(*,*) "Stability criteria is met"
	ELSE
		WRITE(*,*) "The program is going to exit since stability criteria is not met"
		CALL EXIT(0)
	END IF
	!-------------------------------------------------------------------------------------------------------- 
		
	! Start Timer	
	CALL CPU_TIME(time1)	
	
	! Initial and Boundary COnditions
	DO j=1,N+1			
		DO i=1,M+1							
			IF ((i.EQ.0).OR.(i.EQ.M+1)) THEN
				T0(i,j) = f(j)			
			ELSE IF (j.EQ.N+1) THEN
				T0(i,j) = 100
			ELSE
				T0(i,j) = 0
			END IF
			Tf(i,j) = T0(i,j)
		END DO
	END DO
	
	! Write Temperature values
	WRITE(filename,'(a,i4.4,a)') 'TecPlot2D_Temperature_',iter,'.tec'	
	OPEN ( UNIT = iunit, FILE = filename, FORM = 'formatted', ACCESS = 'sequential', STATUS = 'replace', IOSTAT = ierr )	
	IF ( ierr /= 0 ) THEN
		WRITE ( *, '(a)' ) '  Error opening file : tecplot_2D '
		STOP
	END IF
	WRITE ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
	WRITE ( iunit, '(a)' ) 'Variables=' // trim ( '"Time","X","Y","T"' )
	WRITE ( iunit, '(a,i6,a,i6,a,a)' ) 'Zone I=', M+1, ', J=', N+1, ', F=POINT'
	DO j=1,N+1;DO i=1,M+1; WRITE(iunit,'(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6)')  t, x(j), y(i), Tf(i,j) ;END DO;END DO
	CLOSE(UNIT = iunit)	 ! Close file
	
	t = t+dt ! Initial stepping through time
				
	DO WHILE ((t.LT.tmax).OR.(t.EQ.tmax))	
		!$omp parallel do default (none) shared (Tf,dt_dy2,dt_dx2) &
		!$omp private (i,j,A,T0)
		DO j=2,N		
			DO i=2,M					
				Tf(i,j) = ((A(i,j+1)+ A(i,j)/2)*dt_dx2*(T0(i,j+1)-T0(i,j)))-&
						  ((A(i,j)+ A(i,j-1)/2)*dt_dx2*(T0(i,j)-T0(i,j-1)))+&
						  (dt_dy2*((A(i+1,j)*T0(i+1,j))-(2*A(i,j)*T0(i,j))+(A(i-1,j)*T0(i-1,j))))+ T0(i,j)				
			END DO					
		END DO
		!$omp end parallel do
		T0 = Tf	
		t = t+dt
		iter = iter+1
		
		! Write Temperature to a file
		IF(MOD(iter,iwrite).eq.0) THEN
			WRITE(filename,'(a,i4.4,a)') 'TecPlot2D_Temperature_',iter,'.tec'	
			OPEN ( UNIT = iunit, FILE = filename, FORM = 'formatted', ACCESS = 'sequential', STATUS = 'replace', IOSTAT = ierr )	
			IF ( ierr /= 0 ) THEN
				WRITE ( *, '(a)' ) '  Error opening file : tecplot_2D '
				STOP
			END IF
			WRITE ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
			WRITE ( iunit, '(a)' ) 'Variables=' // trim ( '"Time","X","Y","T"' )
			WRITE ( iunit, '(a,i6,a,i6,a,a)' ) 'Zone I=', M+1, ', J=', N+1, ', F=POINT'
			DO j=1,N+1;DO i=1,M+1; WRITE(iunit,'(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6)')  t, x(j), y(i), Tf(i,j) ;END DO;END DO
			CLOSE(UNIT = iunit)	 ! Close file
		END IF	
		
	END DO
	
	! Finish Timer
	CALL CPU_TIME(time2)
	WRITE(*,'(a)') "# Simulation Finnished "
	WRITE(*,'(a15,f14.10)') "# Total WTime = ",  time2 - time1		
	
END PROGRAM fde_2D_heattransfer

