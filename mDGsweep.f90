! --------------------------------------------------------------------
! Modal DG sweep update for DG/FV hybrid 1D Advection
! By: Devin Light
! --------------------------------------------------------------------

SUBROUTINE mDGsweep(rhoq,rhop,u,uedge,dxel,nelem,N,wghts,nodes,DG_C,DG_LUC,DG_L,DG_DL,IPIV,dt,dodghybrid,doposlimit)
	USE mDGmod
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- External functions
	REAL(KIND=DOUBLE), EXTERNAL :: B ! RHS function for evolution ODE for kth expansion coefficent

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem,N ! Number of elements, highest degree of Legendre polynomial being used
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt ! Element width, Finite Volume sub-cell width, time step
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes ! Quadrature weights and node locations
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_C, DG_LUC,DG_L,DG_DL! C matrix, used to xfer between grids, LU decomp of C
	INTEGER, DIMENSION(0:N), INTENT(IN):: IPIV ! Pivot array for RHS when using DG_LUC
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(IN) :: u ! Velocities at quadrature locations within each element
	REAL(KIND=DOUBLE), DIMENSION(1:nelem), INTENT(IN) :: uedge ! Edge velocities at RHS of each element
	LOGICAL, INTENT(IN) :: dodghybrid,doposlimit
	! --- Outputs
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(INOUT) :: rhoq,rhop ! Soln as sub-cell averages within each element at FV cell centers

	! --- Local Variables
	INTEGER :: i,j,k
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqBAR,rpBAR,foorq ! Reshaped cell averages
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1) :: A,A1,A2,R,R1,R2 ! Reshaped DG coefficent matricies
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: utild ! Reshaped velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: uedgetild ! Periodically extended edge velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flxrq, flxrp! Array of fluxes F(j,j+1) (flx(0) is left flux at left edge of domain)
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: fcfrq,fcfrp ! Flux correction factors for positivity limiting

	! -- DGESV parameters
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: FOO_y
	REAL(KIND=DOUBLE), DIMENSION(1:nelem) :: hold
	INTEGER :: ierr
	
	REAL(KIND=4) :: ta,tb
	REAL(KIND=4), DIMENSION(2) :: time
	
	! ########################################
    ! A(k,j) gives a_k(t) in the jth element for rhoq 
	! R(k,j) gives a_k(t) in the jth element for rho
    ! ########################################

    ! Reform incoming values to be more convienent
    DO j=1,nelem
        rqBAR(:,j) = rhoq(1+(N+1)*(j-1) : (N+1)*j)
		rpBAR(:,j) = rhop(1+(N+1)*(j-1) : (N+1)*j)
        utild(:,j) = u(1+(N+1)*(j-1) : (N+1)*j)
    END DO
		uedgetild(1:nelem) = uedge(1:nelem)

    ! For hybrid, values incoming assumed to be cell averages, project onto basis, giving coefficents a_k(t) 
	! which corresp. to series that, when averaged over the evenly spaced subcells in each element, result in the incoming values.
	IF(dodghybrid) THEN

		! Construct RHS using IPIV array
		DO i=0,N
			hold = rqBAR(i,1:nelem)
			rqBAR(i,1:nelem) = rqBAR(IPIV(i)-1,1:nelem)
			rqBAR(IPIV(i)-1,1:nelem) = hold

			hold = rpBAR(i,:)
			rpBAR(i,:) = rpBAR(IPIV(i)-1,:)
			rpBAR(IPIV(i)-1,:) = hold
		ENDDO

		DO j=1,nelem

			FOO_y = 0D0
			! Solve Ly=RHS for y
			FOO_y(0) = rqBAR(0,j)

			DO k=1,N
				FOO_y(k) = rqBAR(k,j) - SUM(DG_LUC(k,0:k-1)*FOO_y(0:k-1))
			ENDDO

			! Solve Ux=y for x
			A(N,j) = (1D0/DG_LUC(N,N))*FOO_y(N)
			DO k=N-1,0,-1
				A(k,j) = (1D0/DG_LUC(k,k))*(FOO_y(k) - SUM(DG_LUC(k,k+1:N)*A(k+1:N,j)))
			ENDDO
		ENDDO

		DO j=1,nelem
			FOO_y = 0D0
			! Solve Ly=RHS for y
			FOO_y(0) = rpBAR(0,j)
			DO k=1,N
				FOO_y(k) = rpBAR(k,j) - SUM(DG_LUC(k,0:k-1)*FOO_y(0:k-1))
			ENDDO
			! Solve Ux=y for x
			R(N,j) = (1D0/DG_LUC(N,N))*FOO_y(N)
			DO k=N-1,0,-1
				R(k,j) = (1D0/DG_LUC(k,k))*(FOO_y(k) - SUM(DG_LUC(k,k+1:)*R(k+1:,j)))
			ENDDO
		ENDDO

	ELSE
	! Otherwise, just use reshaped values
!	A(:,1:nelem) = rqBAR(:,1:nelem)
		STOP 'ERROR in mDGsweep.f90: cannot call modal DG without cell averages'
	END IF

	IF(doposlimit) THEN
		DO j=1,nelem
!			A(0,j) = MAX(A(0,j),0D0)
!			write(*,*) 2D0*A(0,j) - (2D0/(N+1))*SUM(rqBAR(:,j))
			IF(A(0,j) .lt. 0D0) THEN
				write(*,*) '!!! ~~ WARNING:: MASS AFTER PROJECTION IS NEGATIVE! j=',j
				write(*,*) 'DG MASS:',2D0*A(0,j)
				write(*,*) 'AV MASS:',(2D0/(DBLE(N+1)))*SUM(rqBAR(:,j))
			ENDIF
		ENDDO		
	ENDIF
	! -- Enforce periodicity
	A(:,0) = A(:,nelem)
	A(:,nelem+1) = A(:,1)
	R(:,0) = R(:,nelem)
	R(:,nelem+1) = R(:,1)
	uedgetild(0) = uedgetild(nelem)

    ! #######################
    ! Time step using SSPRK3
    ! #######################

!	DO j=0,nelem
!		flxrq(j) = NUMFLUX(A,j,j+1,uedgetild(j),N,nelem)
!		flxrp(j) = NUMFLUX(R,j,j+1,uedgetild(j),N,nelem)
!	ENDDO

	CALL NUMFLUX(A,uedgetild,N,nelem,flxrq)
	CALL NUMFLUX(R,uedgetild,N,nelem,flxrp)
	
	! Determine flux correction factors
	fcfrq = 1d0 ! By default, do no corrections
	fcfrp = 1d0
	IF(doposlimit) THEN
		CALL FLUXCOR(A,A,flxrq,DG_C,dxel,dt,nelem,N,1,fcfrq)
		CALL FLUXCOR(R,R,flxrp,DG_C,dxel,dt,nelem,N,1,fcfrp)
	END IF

	! -- First step
	DO j=1,nelem	
		DO k=0,N
			A1(k,j) = A(k,j) + dt*B(A,flxrq,utild,dxel,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(:,k))
			R1(k,j) = R(k,j) + dt*B(R,flxrp,utild,dxel,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(:,k))
		END DO
	END DO

	! -- Enforce periodicity
    A1(:,0) = A1(:,nelem)
	A1(:,nelem+1) = A1(:,1)
	R1(:,0) = R1(:,nelem)
	R1(:,nelem+1) = R1(:,1)

	! Compute fluxes
!	DO j=0,nelem
!		flxrq(j) = NUMFLUX(A1,j,j+1,uedgetild(j),N,nelem)
!		flxrp(j) = NUMFLUX(R1,j,j+1,uedgetild(j),N,nelem)
!	END DO

	CALL NUMFLUX(A1,uedgetild,N,nelem,flxrq)
	CALL NUMFLUX(R1,uedgetild,N,nelem,flxrp)


	! Determine flux correction factors
	IF(doposlimit) THEN
		fcfrq = 1D0
		fcfrp = 1D0
		CALL FLUXCOR(A1,A,flxrq,DG_C,dxel,dt,nelem,N,2,fcfrq)
		CALL FLUXCOR(R1,R,flxrp,DG_C,dxel,dt,nelem,N,2,fcfrp)
	END IF

	! -- Second step
	DO j = 1, nelem
		DO k = 0, N
			A2(k,j) = (0.75D0)*A(k,j) + (0.25D0)*(A1(k,j) + dt*B(A1,flxrq,utild,dxel,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(:,k)))
			R2(k,j) = (0.75D0)*R(k,j) + (0.25D0)*(R1(k,j) + dt*B(R1,flxrp,utild,dxel,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(:,k)))
		END DO
	END DO

	! -- Enforce periodicity
	A2(:,0) = A2(:,nelem)
	A2(:,nelem+1) = A2(:,1)
	R2(:,0) = R2(:,nelem)
	R2(:,nelem+1) = R2(:,1)

	! Compute fluxes
!	DO j=0,nelem
!		flxrq(j) = NUMFLUX(A2,j,j+1,uedgetild(j),N,nelem)
!		flxrp(j) = NUMFLUX(R2,j,j+1,uedgetild(j),N,nelem)
!	END DO
	CALL NUMFLUX(A2,uedgetild,N,nelem,flxrq)
	CALL NUMFLUX(R2,uedgetild,N,nelem,flxrp)


	! Determine flux correction factors
	IF(doposlimit) THEN
		fcfrq = 1D0
		fcfrp = 1D0
		CALL FLUXCOR(A2,A,flxrq,DG_C,dxel,dt,nelem,N,3,fcfrq)
		CALL FLUXCOR(R2,R,flxrp,DG_C,dxel,dt,nelem,N,3,fcfrp)
	END IF

	! -- Third step
	DO j = 1,nelem
		DO k = 0, N
			A(k,j) = (1D0/3D0)*A(k,j) + (2D0/3D0)*(A2(k,j) + dt*B(A2,flxrq,utild,dxel,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(:,k)))
			R(k,j) = (1D0/3D0)*R(k,j) + (2D0/3D0)*(R2(k,j) + dt*B(R2,flxrp,utild,dxel,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(:,k)))
		END DO
	END DO

	! #########
	! END SPPRK3 TIMESTEP
	! #########

	! #######
	! BEGIN CELL AVERAGING
	! #######

	IF(dodghybrid) THEN
    ! After time stepping, use DG_C to average the series expansion to cell-averaged values
    ! on evenly spaced grid for finite volumes
	    DO j = 1,nelem
			rqBAR(:,j) = MATMUL(DG_C,A(:,j))
			rpBAR(:,j) = MATMUL(DG_C,R(:,j))
	    END DO
	ELSE
	! Otherwise, just send back DG coefficents
!		DO j=1,nelem
!			rqBAR(:,j) = A(:,j)
!		END DO
		STOP 'ERROR in mDGsweep.f90: cannot call modal DG without cell averages'
	END IF

go to 999
	IF(doposlimit) THEN
		DO j=1,nelem

		IF(A(0,j) .lt. 0D0) THEN
			write(*,*) 'j=',j,'Ex Mj=',2D0*A(0,j)
		END IF

		IF(dodghybrid) THEN
			IF(SUM(rqBAR(:,j)) .lt. 0D0) THEN
				write(*,*) 'j=',j,'Av Mj=',(2D0/(N+1))*SUM(rqBAR(:,j))
			END IF
		END IF


		END DO
	END IF
999 continue
	! #######
	! END CELL AVERAGING
	! #######

	! Mass filling to prevent negative subcell averages within each element
	IF(doposlimit) THEN
		IF(.NOT. dodghybrid) THEN
			write(*,*) 'WARNING: TRYING TO MASS FILL EXPANSION COEFFICENTS -- THIS IS BAD!'
		END IF

		! Write pre- mass redistribution cell values

		CALL MFILL(rqBAR,N,nelem)
		CALL MFILL(rpBAR,N,nelem)

		DO j=1,nelem
		if(minval(rqBAR(:,j)) .lt. 0D0) then
			write(*,*) 'NEGATIVE CELL AVERAGE IN ELEMENT j=',j
			write(*,*) 'MINIMUM IS:',minval(rqBAR(:,j))
		end if
		end do
	END IF

    ! Reform original rp and rq vectors to send back
	! -- Note: There are 2 ghost cells that we don't update
    DO j=1,nelem
        rhoq(1+(N+1)*(j-1) : (N+1)*j) = rqBAR(:,j)
		rhop(1+(N+1)*(j-1) : (N+1)*j) = rpBAR(:,j)
    END DO



END SUBROUTINE mDGsweep

REAL(KIND=KIND(1D0)) FUNCTION B(Ain,flx,utild,dxel,wghts,nodes,k,j,nelem,N,fluxcf,Leg,dLeg) 
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem, N
	INTEGER, INTENT(IN) :: k,j ! Mode number, element number
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes,dLeg
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1), INTENT(IN) :: Ain
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: utild
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: Leg

	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: flx ! Array of element edge fluxes flx(j) = F(j,j+1)
	
	! Flux correction factor. fluxcf(j)=reduction in flux coming through face (j+1/2)
	! fluxcf(0)=reduction factor for flux coming through left boundary of domain.
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: fluxcf
 
	! --- Local variables
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: HOLDER1
	INTEGER :: i

	DO i=0,N
		HOLDER1(i) = wghts(i)*utild(i,j)*dLeg(i)*SUM(Ain(:,j)*Leg(:,i))
	END DO

	B = SUM(HOLDER1)


	IF(k .eq. 0) THEN
		B = B - fluxcf(j)*flx(j) + fluxcf(j-1)*flx(j-1)
	ELSE
		B = B - flx(j) + ((-1D0)**k)*flx(j-1)
	END IF
	
	B = B*((2*k+1)/dxel)

END FUNCTION B

REAL(KIND=KIND(1D0)) FUNCTION NUMFLUX(Ain,jleft,jright,u,N,nelem)
	USE mDGmod
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: N,nelem,jleft,jright
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1), INTENT(IN) :: Ain
	REAL(KIND=DOUBLE), INTENT(IN) :: u ! Velocity at element interface

	! -- Local variables
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: foo
	INTEGER :: k

	foo = 0d0

	! -- Note that it is assumed that left and right hand limits of velocities are the same
	! -- Also note that because the Legendre basis functions are all +1 at xi = 1D0, and +-1 at xi = -1D0, there may be
	!	 some optimization that we can do here..
	IF(u>0D0) THEN
		NUMFLUX = u*SUM(Ain(:,jleft))
	ELSE
		DO k=0,N
			foo(k) = ((-1D0)**k)*Ain(k,jright) ! (-1**k) since legendre polys alternate sign at -1
		END DO
		NUMFLUX = u*SUM(foo)
	END IF

END FUNCTION NUMFLUX

SUBROUTINE NUMFLUX(A,u,N,nelem,flx)
	USE mDGmod
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: N,nelem
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1), INTENT(IN) :: A
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: u

	! -- Outputs	
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flx

	! -- Local variables
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: foo
	INTEGER :: k,j

	DO j=0,nelem
		IF(u(j) .ge. 0D0) then
			flx(j) = u(j)*SUM(A(:,j))
		ELSE
			DO k=0,N
				foo(k) = ((-1D0)**k)*A(k,j+1)
			ENDDO
			flx(j) = u(j)*SUM(foo)
		ENDIF
	ENDDO

END SUBROUTINE NUMFLUX


SUBROUTINE FLUXCOR(Acur,Apre,flx,DG_C,dxel,dt,nelem,N,substep,fluxcf)
	! Computes flux reductions factors to prevent total mass within each element from going negative
	! Outputs fluxcf. fluxcf(j) is the reduction factor for the right face of element j, 
	!  with fluxcf(0) being the factor for the left domain interface
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: N, nelem,substep
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1), INTENT(IN) :: Acur,Apre
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: flx
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_C
	! -- Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(OUT) :: fluxcf
	! -- Local variables
	REAL(KIND=DOUBLE) :: Pj,Qj,eps
	REAL(KIND=DOUBLE), DIMENSION(0:nelem+1) :: R ! Reduction ratio for outward fluxes so that element j has non-negative values (1D0 indicates no limiting needed)
	INTEGER :: j,k

	eps = 1D-16 ! Small parameter used to ensure no division by 0

	DO j=1,nelem
		! Compute maximum allowable flux out of element j based on which substep of ssprk3 we are on
		SELECT CASE(substep)
			CASE(1) 
				Qj = (dxel/dt)*Acur(0,j)
!				Qj = (dxel**2/dt)*SUM(MATMUL(DG_C,Acur(:,j)))
		IF(Qj .lt. 0D0) THEN
			write(*,*) 'Stage 1 Qj < 0! j=',j,'Qj=',Qj
			write(*,*) Acur(0,j)
		END IF

			CASE(2) 
				Qj = (dxel/dt)*(3D0*Apre(0,j) + 1D0*Acur(0,j))
!				Qj = (dxel**2/dt)*( 3D0*SUM(MATMUL(DG_C,Apre(:,j))) + 1D0*SUM(MATMUL(DG_C,Acur(:,j))) )
		IF(Qj .lt. 0D0) THEN
			write(*,*) 'Stage 2 Qj < 0! j=',j,'Qj=',Qj
			write(*,*) Apre(0,j),Acur(0,j)
		END IF

			CASE(3) 
				Qj = (dxel/(2D0*dt))*(1D0*Apre(0,j) + 2D0*Acur(0,j))
!				Qj = (dxel**2/(2D0*dt))*( 1D0*SUM(MATMUL(DG_C,Apre(:,j))) + 2D0*SUM(MATMUL(DG_C,Acur(:,j))) )
		IF(Qj .lt. 0D0) THEN
			write(*,*) 'Stage 3 Qj < 0! j=',j,'Qj=',Qj
			write(*,*) Apre(0,j),Acur(0,j)
		END IF

		END SELECT

		! Compute actual flux out of element j
		Pj = DMAX1(0D0,flx(j)) - DMIN1(0D0,flx(j-1)) + eps

		! Compute reduction ratio
		R(j) = DMIN1(1D0,Qj/Pj)
	END DO
	! Periodicity
	R(0) = R(nelem)
	R(nelem+1) = R(1)

	! Compute flux corection factors
	fluxcf = R(0:nelem)
	DO j=0,nelem
		! If flux at right edge is negative, use limiting ratio in element to the right of current one
		! (since that is where we are pulling mass from)
		IF(flx(j) < 0D0) THEN
			fluxcf(j) = R(j+1)
		END IF
	END DO

END SUBROUTINE FLUXCOR


SUBROUTINE MFILL(rhoq,N,nelem)
	! Subroutine for mass filling within an element to remove negative cell averaged values
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! -- Inputs
	INTEGER, INTENT(IN) :: N,nelem
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(INOUT) :: rhoq
	! -- Local Variables
	INTEGER :: j,k
	REAL(KIND=DOUBLE) :: r,Mp,Mt

	DO j=1,nelem
		Mp = 0D0
		Mt = 0D0

		DO k=0,N
			Mt = Mt + rhoq(k,j)
			rhoq(k,j) = MAX(0D0,rhoq(k,j)) ! Zero out negative masses
			Mp = Mp + rhoq(k,j)
		ENDDO
!		IF(Mp .eq. 0D0) THEN
!			write(*,*) 'OH NO!! BIG TROUBLE!'
!		ENDIF

!		r = MAX(Mt,0D0)/(Mp+1D-14)
		r = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))
		IF(r .gt. 1D0) THEN
			write(*,*) 'WARNING REDUCTION RATIO > 1.0!!'
		ENDIF
		rhoq(:,j) = r*rhoq(:,j) ! Reduce remaining positive masses by reduction factor

	ENDDO

END SUBROUTINE MFILL
