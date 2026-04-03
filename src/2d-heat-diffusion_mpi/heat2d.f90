!Name = Shubhrasmita Mandal (Smita)
!********************************************************************************************************************************
!Using MPI
!The program finds the temperature of a 2D heat diffusion problem.
!Heat transfer is  steady and the boundary coditions are as follows:
!f(x) = sin^2(2πx/Lx) at upper boundary and f(x) = -sin^2(2πx/Lx) at lower boundary
!Temperature at all other boundary is zero. FLuid inside the boundaries has spatially varying velocities as follows:
!U(x, y) = sin(πx/Lx) cos(πy/Ly) + sin(2πx/Lx) cos(2πy/Ly) ; V (x, y) = − cos(πx/Lx) sin(πy/Ly) − cos(2πx/Lx) sin(2πy/Ly)
!Heat diffusion equation is given as: alpha*U*(dT/dx)+alpha*V(dT/dy) -(d^2T/dx^2) -(d^2T/dy^2) -beta*T=0
!Using first and second order centered finite difference discretisation in space.
!*********************************************************************************************************************************
!Steps for compilation:
!To compile the program into executable, following command can be used:
!mpif90 smita.f90 -o smita.exe
!Further, to run the executable, use the following command:
!mpirun -np 6 smita.exe
!*********************************************************************************************************************************

PROGRAM mpi_2d_ht
	IMPLICIT NONE	
	INCLUDE 'mpif.h'
	
	!Variable declaration
    REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979323846 	
	INTEGER :: kx,ky,kz,il,ih,jl,jh,i,j,k
	INTEGER :: px,py,pil,pjl,pih,pjh
	INTEGER :: nprocs,pid,ierr,load_p,pid_start,pid_stop,left,right,tag,iter,q
	REAL (KIND=8) :: dx,dy,meshlx,meshly,rmax,rcurrent
	REAL (KIND=8) :: alpha,beta,time1,time2	
	REAL(KIND=8), ALLOCATABLE :: aw(:,:),ae(:,:),ap(:,:),an(:,:),as(:,:),u(:,:),v(:,:)
	REAL(KIND=8), ALLOCATABLE :: x(:),y(:),t(:,:),b(:,:),f(:),xval(:),yval(:),Temp(:,:),rsdl(:,:)
	REAL(KIND=8), ALLOCATABLE :: r(:,:),kt(:,:),ROW_TYPE(:),l_gr(:),g_gr(:)
	INTEGER, ALLOCATABLE :: reqs(:), stats(:,:),rcounts1(:),rcounts2(:),disp1(:),disp2(:)
	
	CALL MPI_INIT(ierr) ! Initialize MPI Program
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr) ! Get Size of Communicator (Number of procs)s
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)    ! Get Rank of Communicator (proc id)
		
	!-----------------------------------------
	! USER INPUT
	!-----------------------------------------
	kx = 100; ky = 100; kz = 100		! MESHSIZE	
	meshlx = 1.0; meshly = 1.0 	        ! DOMAIN SIZE
	rmax = 0.000001_8					! MAX RESIDUAL 	  
    alpha = 80.0_8
	beta = 0.5_8
	tag = 101	
	
	!-----------------------------------------
	! Grid Size
	!-----------------------------------------
	dx = meshlx/REAL(kx-1)
	dy = meshly/REAL(ky-1)
	
	! Interior Indices
	il = 2; jl = 2;
	ih = kx-1; jh = ky-1	
	
	!-------------------------------------------------
	! Divide the domain between the processes
	!-------------------------------------------------
	load_p = kx/nprocs
	IF (pid.EQ.0) THEN
		pid_start = (pid * load_p) + 1
		pid_stop = (pid_start + load_p)
	ELSE IF (pid.EQ.(nprocs-1)) THEN
		pid_start = (pid * load_p)
		pid_stop = kx
	ELSE
		pid_start = (pid * load_p)
		pid_stop = (pid_start + load_p) + 1		
	END IF	
	
	!---------------------------------------------------
	! No. of nodes in x and y direction for each process
	!---------------------------------------------------
	px = (pid_stop - pid_start) + 1
	py = ky
	
	!-----------------------------------------
	! Interior indices for each process
	!-----------------------------------------
	IF (pid.EQ.0) THEN
		pil = il;pjl = jl
		pih = pid_stop;pjh = jh	
	ELSE IF(pid.EQ.nprocs-1) THEN
		pil = pid_start;pjl = jl
		pih = ih;pjh = jh
	ELSE
		pil = pid_start;pjl = jl
		pih = pid_stop;pjh = jh
	ENDIF
		
	!--------------------------------------------
	!Initialize variables for Neighboring process
	!--------------------------------------------
	IF (pid.EQ.0) THEN
		left = MPI_PROC_NULL
		right = pid + 1
		ALLOCATE(reqs(2))
		ALLOCATE(stats(MPI_STATUS_SIZE,2))
    ELSEIF (pid.EQ.(nprocs-1))	THEN
		right = MPI_PROC_NULL
		left = pid - 1	
		ALLOCATE(reqs(2))
		ALLOCATE(stats(MPI_STATUS_SIZE,2))
	ELSE
		left = pid - 1
		right = pid + 1	
		ALLOCATE(reqs(4))
		ALLOCATE(stats(MPI_STATUS_SIZE,4))
	END IF
	
	!-----------------------------------------
	! Allocate memory
	!-----------------------------------------
	ALLOCATE(aw(ky,kx))
	ALLOCATE(ae(ky,kx))
	ALLOCATE(an(ky,kx))
	ALLOCATE(as(ky,kx))
	ALLOCATE(ap(ky,kx))
	ALLOCATE(t(ky,kx))	           
	ALLOCATE(b(ky,kx))
	ALLOCATE(u(ky,kx))
	ALLOCATE(v(ky,kx))
	ALLOCATE(r(ky,kx))
	ALLOCATE(kt(ky,kx))
	ALLOCATE(f(kx))
	ALLOCATE(x(kx))
	ALLOCATE(y(ky))	
	ALLOCATE(ROW_TYPE(kx))
    ALLOCATE(l_gr(ky*kx))
    ALLOCATE(g_gr(ky*kx))
	ALLOCATE(rcounts1(nprocs))
	ALLOCATE(rcounts2(nprocs))
	ALLOCATE(disp1(nprocs))
	ALLOCATE(disp2(nprocs))
	
	
	IF (pid.EQ.0) THEN
		ALLOCATE(Temp(ky,kx))
		ALLOCATE(rsdl(ky,kx))
		ALLOCATE(xval(kx))
		ALLOCATE(yval(ky))
	ENDIF	
	
	!--------------------------------
	! Initialize Solution Parameters
	!--------------------------------	
	! x coordinate
	DO i = pid_start,pid_stop
		x(i) = (i-1)*dx	
		! Upper/Lower Boundary conditions
        f(i) = (sin((2 * pi * x(i))/meshlx)) ** 2		
	END DO
	
	! y coordinate
	DO j = 1,py
		y(j) = (j-1)*dy											
	END DO	
	
	!Calculate coefficient matrices
	DO j = pid_start,pid_stop
		DO i = 1,ky
			IF (i.EQ.1) THEN 
				t(i,j) = -f(j)
				kt(i,j) = -f(j)
				b(i,j) = -f(j)
			ELSEIF (i.EQ.py) THEN
				t(i,j) = f(j)
				kt(i,j) = f(j)
				b(i,j) = f(j)
			ELSE
				t(i,j) = 0.0
				kt(i,j) = 0.0
				b(i,j) = 0.0				
			ENDIF
		END DO
	END DO	

	CALL GetCoefficients(pid_start,pid_stop,kx,ky,meshlx,meshly,x,y,dx,dy,alpha,beta,u,v,aw,ae,an,as,ap)		
	!-----------------------------------------------------------------
	!-----------------------------------------------------------------
	
	IF (pid.EQ.0) THEN	
		WRITE(*,'(a)') "# Initializing Simulation - Jacobian Solver for Steady-State Heat "		
		WRITE(*,'(a20,i14.1)') "# Nx = ", kx
		WRITE(*,'(a20,i14.1)') "# Ny = ", ky
		WRITE(*,'(a20,f14.1)') "# Lx = ", meshlx
		WRITE(*,'(a20,f14.1)') "# Ly = ", meshly
		WRITE(*,'(a20,E20.10)') "# Max Residual = ", rmax	
		WRITE(*,'(a)') "# Solution Initialized "
	END IF
	
	!------------------------------------------------------------------
	! Solver Start: Jacobian Method
	!------------------------------------------------------------------
	
	! Start Timer	
	CALL CPU_TIME(time1)
	
	iter = 0
	
	! Get Initial Residual
	CALL GetResidual(kx,ky,pil,pih,pjl,pjh,aw,ae,an,as,ap,b,t,r)	
	
	! Gather residuals to Root
	IF (pid.EQ.0) THEN
		DO j = 1,kx
			DO i = 1,ky				
				rsdl(i,j) = r(i,j)
			END DO
		END DO
	END IF	
	k = 0
	DO j = pid_start+2,pid_stop
		DO i = 1,ky
			k = k+1
			l_gr(k) = r(i,j)			
		ENDDO
	ENDDO
	
	IF(pid.EQ.0) THEN
	
		DO i = 0,nprocs-1
			IF(i.EQ.0) THEN
				rcounts1(i+1) = ((load_p + 1)-2)*ky
				rcounts2(i+1) = ((load_p + 1)-2)
				disp1(i+1) = 0
				disp2(i+1) = 0
			ELSEIF ((i.GT.0).AND.(i.LT.nprocs-1)) THEN
				rcounts1(i+1) = ((load_p + 2)-2)*ky
				rcounts2(i+1) = ((load_p + 2)-2)
				disp1(i+1) = ((load_p*i) - 1)*ky
				disp2(i+1) = (load_p*i) - 1
			ELSE
				rcounts1(i+1) = (((kx - (load_p*i)) + 1)-2)*ky
				rcounts2(i+1) = (((kx - (load_p*i)) + 1)-2)
				disp1(i+1) = ((load_p*i) -1)*ky
				disp2(i+1) = (load_p*i) -1
			ENDIF
			
		ENDDO
	ENDIF	
	CALL MPI_GATHERV(l_gr(1:k),k,MPI_REAL8,g_gr(1:ky*kx),rcounts1,disp1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	
	IF (pid.EQ.0) THEN
		k=1
		DO j=pid_start+2,kx
			DO i=1,ky
				rsdl(i,j) = g_gr(k)
				k = k+1
			ENDDO
		ENDDO
	ENDIF	
	
	! Calculate Domain averaged residual for stopping criterion	
	IF (pid.EQ.0) THEN			
		rcurrent = SUM(SUM(ABS(rsdl(jl:jh,il:ih)),1),1) / ((kx-2)*(ky-2))	
		WRITE(*,*) "Residual at iteration ",iter,":",rcurrent		
	END IF	
	
	! Broascast average residual value for the domain to all processes		
	CALL MPI_BCAST(rcurrent,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	
	DO WHILE ((iter.LT.kz).AND.(rcurrent.GT.rmax))
		
		!------------------------------------------
		! Communication between processes 
		!------------------------------------------
		DO i = pjl,pjh			
			q = 0
			! Communication of non-contiguous data			
			CALL MPI_TYPE_VECTOR(kx,1,ky, MPI_REAL8 ,ROW_TYPE ,ierr )
			CALL MPI_TYPE_COMMIT(ROW_TYPE ,ierr)
			CALL MPI_BARRIER (MPI_COMM_WORLD,ierr)			
			DO j = pil,pih						
				IF ((j.EQ.(pil + 1)).AND.(left.NE.MPI_PROC_NULL)) THEN
					q = q + 1						
					CALL MPI_ISEND(t(i,j),1,ROW_TYPE,left,tag,MPI_COMM_WORLD,reqs(q),ierr)
										
				ELSEIF ((j.EQ.(pih - 1)).AND.(right.NE.MPI_PROC_NULL)) THEN
					q = q + 1						
					CALL MPI_ISEND(t(i,j),1,ROW_TYPE,right,tag,MPI_COMM_WORLD,reqs(q),ierr)					
					
				ELSEIF ((j.EQ.pil).AND.(left.NE.MPI_PROC_NULL)) THEN
					q = q + 1						
					CALL MPI_IRECV(t(i,j),1,ROW_TYPE,left,tag,MPI_COMM_WORLD,reqs(q),ierr)
										
				ELSEIF ((j.EQ.pih).AND.(right.NE.MPI_PROC_NULL)) THEN
					q = q + 1
					CALL MPI_IRECV(t(i,j),1,ROW_TYPE,right,tag,MPI_COMM_WORLD,reqs(q),ierr)						
				END IF				
				
			END DO
			! Waiting for all the requests to complete
			CALL MPI_WAITALL ( q , reqs , stats , ierr )
			CALL MPI_TYPE_FREE(ROW_TYPE )				
		END DO
		
		! Calculate temperature at each node for kth iteration			
		CALL GetTemperature(kx,ky,pil,pih,pjl,pjh,aw,ae,an,as,ap,b,t,kt)
		
		DO j = il,ih
			DO i = jl,jh
				t(i,j) = kt(i,j)
			END DO
		END DO
				
		!-------------------------------------------------
		! Calculate residual at each node
        !-------------------------------------------------		
		CALL GetResidual(kx,ky,pil,pih,pjl,pjh,aw,ae,an,as,ap,b,kt,r)
		
		! Gather residuals and temperatures to Root
		IF (pid.EQ.0) THEN
			DO j = 1,kx
				DO i = 1,ky
					Temp(i,j) = kt(i,j)
					rsdl(i,j) = r(i,j)
				END DO
			END DO
		END IF	
			
		k = 0
		
		DO j = pid_start+2,pid_stop
			
			DO i = 1,ky
				k = k+1
				l_gr(k) = kt(i,j)				
			ENDDO
		ENDDO		
		CALL MPI_GATHERV(l_gr(1:k),k,MPI_REAL8,g_gr(1:ky*kx),rcounts1,disp1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
		
		IF (pid.EQ.0) THEN
			k=1
			DO j=pid_start+2,kx
				DO i=1,ky
					Temp(i,j) = g_gr(k)
					k = k+1
				ENDDO
			ENDDO
		ENDIF
		
		k = 0
		DO j = pid_start+2,pid_stop
			DO i = 1,ky
				k = k+1
				l_gr(k) = r(i,j)				
			ENDDO
		ENDDO		
		CALL MPI_GATHERV(l_gr(1:k),k,MPI_REAL8,g_gr(1:ky*kx),rcounts1,disp1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
		
		IF (pid.EQ.0) THEN
			k=1
			DO j=pid_start+2,kx
				DO i=1,ky
					rsdl(i,j) = g_gr(k)
					k = k+1
				ENDDO
			ENDDO
		ENDIF	
						
		IF (pid.EQ.0) THEN			
			! Calculate Domain averaged residual for stopping criterion
			rcurrent = SUM(SUM(ABS(rsdl(jl:jh,il:ih)),1),1) / ((kx-2)*(ky-2))
			WRITE(*,*) "Residual at iteration ",iter+1,":",rcurrent
		END IF		
		
		! Broascast average residual value for the domain to all processes				
		CALL MPI_BCAST(rcurrent,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)		
		iter = iter + 1
	END DO	
		
	!------------------------
	! Finish Timer
	!------------------------
	CALL CPU_TIME(time2)
	
	!------------------------------------
	! Gather x and y node values to Root
	!------------------------------------
	IF (pid.EQ.0) THEN	
		WRITE(*,'(a)') "# Simulation Finished "
	    WRITE(*,'(a15,f14.10)') "# Total WTime = ",  time2 - time1		
		
		DO i = 1,kx
			xval(i) = x(i)
		END DO
			
		DO j = 1,ky
			yval(j) = y(j)
			
		END DO		
				
	END IF	
	
	! Gather x and y node values to Root
	k = (pid_stop - (pid_start+2)) + 1	    
	CALL MPI_GATHERV(x(pid_start+2:pid_stop),k,MPI_REAL8,xval(il+1),rcounts2,disp2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)	
	
	!--------------------
	! Plot Data
	!--------------------
	IF (pid.EQ.0) THEN		
		CALL tecplot_2D(kx,ky,xval,yval,Temp,rcurrent)
	END IF	
	
	CALL MPI_FINALIZE(ierr)
	
END PROGRAM mpi_2d_ht

!-------------------------------------------------
! Subroutine to calculate coefficient matrices
!-------------------------------------------------
SUBROUTINE GetCoefficients(pid_start,pid_stop,kx,ky,meshlx,meshly,x,y,dx,dy,alpha,beta,u,v,aw,ae,an,as,ap)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: kx,ky,pid_start,pid_stop
	REAL(KIND = 8), INTENT(IN) :: meshlx,meshly,dx,dy,alpha,beta
	REAL(KIND = 8), DIMENSION(kx), INTENT(IN) :: x
	REAL(KIND = 8), DIMENSION(ky), INTENT(IN) :: y
	REAL(KIND = 8), DIMENSION(ky,kx), INTENT(OUT) :: u,v,aw,ae,an,as,ap
	REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979323846 
	INTEGER :: i,j	
	
	DO j = pid_start,pid_stop
		DO i = 1,ky
			u(i,j) = (sin((pi * x(j))/meshlx) * cos((pi * y(i))/meshly)) + (sin((2 * pi * x(j))/meshlx) * cos((2 * pi * y(i))/meshly))
			v(i,j) = -(cos((pi * x(j))/meshlx) * sin((pi * y(i))/meshly)) - (cos((2 * pi * x(j))/meshlx) * sin((2 * pi * y(i))/meshly))
			aw(i,j) = (-alpha * u(i,j) * dx * (dy**2)) - (2 * (dy**2))								
			ae(i,j) = (alpha * u(i,j) * dx * (dy**2)) - (2 * (dy**2))	
			an(i,j) = (alpha * v(i,j) * dy * (dx**2)) - (2 * (dx**2))									
			as(i,j) = (-alpha * v(i,j) * dy * (dx**2)) - (2 * (dx**2))						
			ap(i,j) = (4 * (dy**2)) + (4 * (dx**2)) - (2 * beta * (dx**2) * (dy**2))
		END DO
	END DO
END SUBROUTINE GetCoefficients

!-------------------------------------------------
! Subroutine to calculate temperature for domain
!-------------------------------------------------
SUBROUTINE GetTemperature(kx,ky,pil,pih,pjl,pjh,aw,ae,an,as,ap,b,t,kt)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: kx,ky,pil,pih,pjl,pjh
	REAL(KIND = 8), DIMENSION(ky,kx), INTENT(IN) :: aw,ae,an,as,ap,b,t
	REAL(KIND = 8), DIMENSION(ky,kx), INTENT(INOUT) :: kt	
	INTEGER :: i,j	
	
	DO j = pil,pih
		DO i = pjl,pjh 					
			kt(i,j) = (-b(i,j) - (aw(i,j)*t(i,j-1)) - (ae(i,j)*t(i,j+1)) - (as(i,j)*t(i-1,j)) - (an(i,j)*t(i+1,j)))/ap(i,j)				
		END DO
	END DO
END SUBROUTINE GetTemperature
		
!-------------------------------------------------
! Subroutine to get Residual Matrix
!-------------------------------------------------
SUBROUTINE GetResidual(kx,ky,pil,pih,pjl,pjh,aw,ae,an,as,ap,b,kt,r)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: pil,pih,pjl,pjh,kx,ky
	REAL(KIND = 8), DIMENSION(ky,kx), INTENT(IN) :: aw,ae,an,as,ap,b,kt
	REAL(KIND = 8), DIMENSION(ky,kx), INTENT(OUT) :: r
	
	INTEGER :: i,j			
	r(:,:) = 0.0
	
	DO j = pil,pih
		DO i = pjl,pjh			
			r(i,j) = -b(i,j) - (aw(i,j)*kt(i,j-1)) - (ae(i,j)*kt(i,j+1)) - (as(i,j)*kt(i-1,j)) - (an(i,j)*kt(i+1,j)) - (ap(i,j)*kt(i,j))			
		END DO
	END DO		
END SUBROUTINE GetResidual

!----------------------------------------------
! Write the data to an output file and screen
!----------------------------------------------
SUBROUTINE tecplot_2D (kx,ky,xval,yval,Temp,rcurrent)	
	IMPLICIT NONE	
	INTEGER, INTENT(IN) :: kx,ky
	REAL(KIND = 8), DIMENSION(kx), INTENT(IN) :: xval
	REAL(KIND = 8), DIMENSION(ky), INTENT(IN) :: yval
	REAL(KIND = 8), DIMENSION(ky,kx), INTENT(IN) :: Temp 
	REAL(KIND = 8), INTENT(IN) :: rcurrent	
	INTEGER :: iunit,i,j,ierr
	REAL(KIND = 8) :: maxt			
	CHARACTER(80), PARAMETER ::  file_name = 'TecPlot2D.tec'
	
	maxt = 0.0_8
	
	OPEN ( UNIT = iunit, FILE = file_name, FORM = 'formatted', ACCESS = 'sequential', STATUS = 'replace', IOSTAT = ierr )
	
	IF ( ierr /= 0 ) THEN
		WRITE ( *, '(a)' ) '  Error opening file : tecplot_2D '
		STOP
	END IF
	
	WRITE ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
	WRITE ( iunit, '(a)' ) 'Variables=' // trim ( '"X","Y","T"' )
	
	WRITE ( iunit, '(a)' ) ' '
	WRITE ( iunit, '(a,i6,a,i6,a,a)' ) 'Zone I=', ky, ', J=', kx, ', F=POINT'
	
	DO i = 1, ky
		DO j = 1, kx
			WRITE ( iunit, '(2f10.3,g15.6)' ) xval(j), yval(i), Temp(i,j)
			IF (Temp(i,j).GT.maxt) THEN
				maxt = Temp(i,j)
			END IF
		END DO
	END DO	
	WRITE(*,*) "Solution Matrix:"
	WRITE(*,*) "Residual at solution domain: ", rcurrent
	WRITE(*,*) "Maximum Temperature at solution domain: ", maxt
	
	CLOSE ( UNIT = iunit )		
	
	END SUBROUTINE tecplot_2D
!----------------------------------------------------------------------------- 		
	
	
	
	