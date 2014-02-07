! --------------------------------------------------------------------
! Modal DG sweep update for DG/FV hybrid 1D Advection
! By: Devin Light
! --------------------------------------------------------------------

SUBROUTINE mDGsweep(rhoq,rhop,u,uedge,dxel,nelem,N,wghts,nodes,DG_C,DG_LUC,DG_L,DG_DL,IPIV,dt, & 
			        dodghybrid,doposlimit,dorhoupdate,time,transient)
	USE mDGmod
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- External functions
	REAL(KIND=DOUBLE), EXTERNAL :: B,tfcn ! RHS function for evolution ODE for kth expansion coefficent

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem,N ! Number of elements, highest degree of Legendre polynomial being used
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt,time ! Element width, Finite Volume sub-cell width, time step, and current time
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes ! Quadrature weights and node locations
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_C, DG_LUC,DG_L,DG_DL! C matrix, used to xfer between grids, LU decomp of C
	INTEGER, DIMENSION(0:N), INTENT(IN):: IPIV ! Pivot array for RHS when using DG_LUC
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(IN) :: u ! Velocities at quadrature locations within each element
	REAL(KIND=DOUBLE), DIMENSION(1:nelem), INTENT(IN) :: uedge ! Edge velocities at RHS of each element

	LOGICAL, INTENT(IN) :: dodghybrid,doposlimit,dorhoupdate,transient
	! --- Outputs
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(INOUT) :: rhoq,rhop ! Soln as sub-cell averages within each element at FV cell centers

	! --- Local Variables
	INTEGER :: i,j,k
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqBAR,rpBAR,foorq ! Reshaped cell averages
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: A,A1,A2,R,R1,R2 ! Reshaped DG coefficent matricies
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: utild,uTmp ! Reshaped velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: uedgetild,uEdgeTmp ! Periodically extended edge velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flxrq, flxrp! Array of fluxes F(j,j+1) (flx(0) is left flux at left edge of domain)
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: fcfrq,fcfrp ! Flux correction factors for positivity limiting
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqQuadVals,rhopQuadVals,qQuadVals
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1) :: rqEdgeVals,rhopEdgeVals,qEdgeVals,ones
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: qOnes

	! -- DGESV parameters
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: FOO_y
	REAL(KIND=DOUBLE), DIMENSION(1:nelem) :: hold
	INTEGER :: ierr
		
	ones = 1D0
	qOnes = 1D0

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
	uedgetild(0) = uedgetild(nelem)

    ! For DG hybrid, values incoming assumed to be cell averages, project onto basis, giving coefficents a_k(t) 
	! which corresp. to series that, when averaged over the evenly spaced subcells in each element, result in the incoming values.

	! Project rho onto basis functions
	! Construct RHS using IPIV array
	DO i=0,N
			hold = rpBAR(i,:)
			rpBAR(i,:) = rpBAR(IPIV(i)-1,:)
			rpBAR(IPIV(i)-1,:) = hold
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

!	! Enforce periodicity
!	R(:,0) = R(:,nelem)
!	R(:,nelem+1) = R(:,1)

	! Project rho*q onto basis functions
	! Construct RHS using IPIV array
	DO i=0,N
		hold = rqBAR(i,1:nelem)
		rqBAR(i,1:nelem) = rqBAR(IPIV(i)-1,1:nelem)
		rqBAR(IPIV(i)-1,1:nelem) = hold
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

!	! -- Enforce periodicity
!	A(:,0) = A(:,nelem)
!	A(:,nelem+1) = A(:,1)


    ! ####################################################
    ! UPDATE rhoQ and rho ; Time step using SSPRK3
    ! ####################################################

	! Update velocities to current time
	uTmp = utild
	uEdgeTmp = uedgetild
	if(transient) then
		uTmp = uTmp*tfcn(time)
		uEdgeTmp = uEdgeTmp*tfcn(time)
	endif

	CALL evalExpansion(R,DG_L,rhopQuadVals,rhopEdgeVals,N,nelem)
	CALL evalExpansion(A,DG_L,rqQuadVals,rqEdgeVals,N,nelem)

!	qQuadVals = rqQuadVals/rhopQuadVals
!	qEdgeVals = rqEdgeVals/rhopEdgeVals

	! Note that first update step for rho uses flux (rho*u) from beginning of timestep (assumed to be absorbed into incoming velocity)
!	CALL NUMFLUX(ones,qEdgeVals,uEdgeTmp,N,nelem,flxrp,flxrq) 
	CALL NUMFLUX(rhopEdgeVals,rqEdgeVals,uEdgeTmp,N,nelem,flxrp,flxrq) 


	! Determine flux correction factors
	fcfrp = 1d0 ! By default, do no corrections
	fcfrq = 1d0 ! By default, do no corrections
	IF(doposlimit) THEN
		CALL FLUXCOR(R,R,flxrp,DG_C,dxel,dt,nelem,N,1,fcfrp)
		CALL FLUXCOR(A,A,flxrq,DG_C,dxel,dt,nelem,N,1,fcfrq)
	END IF

	! -- First step
	DO j=1,nelem	
		DO k=0,N
			A1(k,j) = A(k,j) + (dt/dxel)*B(qQuadVals(:,j),qOnes(:),flxrq,uTmp,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(k,:))
!			R1(k,j) = R(k,j) + (dt/dxel)*B(qOnes(:),qOnes(:),flxrp,uTmp,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(k,:))
R1(k,j) = R(k,j)
		END DO
	END DO

A = A1
R = R1
go to 100
	! -- Enforce periodicity
!	R1(:,0) = R1(:,nelem)
!	R1(:,nelem+1) = R1(:,1)
!   A1(:,0) = A1(:,nelem)
!	A1(:,nelem+1) = A1(:,1)

	! Update velocities to new time
uTmp = utild
uEdgeTmp = uedgetild
if(transient) then
	uTmp = uTmp*tfcn(time+dt)
	uEdgeTmp = uEdgeTmp*tfcn(time+dt)
endif

	CALL evalExpansion(R1,DG_L,rhopQuadVals,rhopEdgeVals,N,nelem)
	CALL evalExpansion(A1,DG_L,rqQuadVals,rqEdgeVals,N,nelem)

	qQuadVals = rqQuadVals/rhopQuadVals
	qEdgeVals = rqEdgeVals/rhopEdgeVals

	! Compute fluxes
	CALL NUMFLUX(rhopEdgeVals,qEdgeVals,uEdgeTmp,N,nelem,flxrp,flxrq) 

	! Determine flux correction factors
	IF(doposlimit) THEN
		fcfrp = 1D0
		fcfrq = 1D0 
		CALL FLUXCOR(R1,R,flxrp,DG_C,dxel,dt,nelem,N,2,fcfrp)
		CALL FLUXCOR(A1,A,flxrq,DG_C,dxel,dt,nelem,N,2,fcfrq)
	END IF

	! -- Second step
	DO j = 1, nelem
		DO k = 0, N
			A2(k,j) = (3D0/4D0)*A(k,j) + (1D0/4D0)*(A1(k,j) + & 
						(dt/dxel)*B(qQuadVals(:,j),rhopQuadVals(:,j),flxrq,uTmp,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(k,:)))
!			R2(k,j) = (3D0/4D0)*R(k,j) + (1D0/4D0)*(R1(k,j) + &
!						(dt/dxel)*B(qOnes(:),rhopQuadVals(:,j),flxrp,uTmp,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(k,:)))
R2(k,j) = R(k,j)
		END DO
	END DO

	! -- Enforce periodicity
!	R2(:,0) = R2(:,nelem)
!	R2(:,nelem+1) = R2(:,1)

	! Update velocities to new time
uTmp = utild
uEdgeTmp = uedgetild
if(transient) then
	uTmp = uTmp*tfcn(time+dt/2d0)
	uEdgeTmp = uEdgeTmp*tfcn(time+dt/2d0)
endif

	CALL evalExpansion(R2,DG_L,rhopQuadVals,rhopEdgeVals,N,nelem)
	CALL evalExpansion(A2,DG_L,rqQuadVals,rqEdgeVals,N,nelem)

	qQuadVals = rqQuadVals/rhopQuadVals
	qEdgeVals = rqEdgeVals/rhopEdgeVals


	! Compute fluxes
	CALL NUMFLUX(rhopEdgeVals,qEdgeVals,uEdgeTmp,N,nelem,flxrp,flxrq) 

	! Determine flux correction factors
	IF(doposlimit) THEN
		fcfrp = 1D0
		fcfrq = 1D0 
		CALL FLUXCOR(R2,R,flxrp,DG_C,dxel,dt,nelem,N,3,fcfrp)
		CALL FLUXCOR(A2,A,flxrq,DG_C,dxel,dt,nelem,N,3,fcfrq)
	END IF

	! -- Third step
	DO j = 1,nelem
		DO k = 0, N
			A(k,j) = A(k,j)/3D0 + 2D0*(A2(k,j) + &
						(dt/dxel)*B(qQuadVals(:,j),rhopQuadVals(:,j),flxrq,uTmp,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(k,:)))/3D0
!			R(k,j) = R(k,j)/3D0 + 2D0*(R2(k,j) + &
!						(dt/dxel)*B(qOnes(:),rhopQuadVals(:,j),flxrp,uTmp,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(k,:)))/3D0
R(k,j) = R(k,j)
		END DO
	END DO

100 continue

	! #########
	! END SPPRK3 TIMESTEP ; 	BEGIN CELL AVERAGING
	! #########

    ! After time stepping, use DG_C to average the series expansion to cell-averaged values
    ! on evenly spaced grid for finite volumes

    DO j = 1,nelem
		rpBAR(:,j) = MATMUL(DG_C,R(:,j))
		rqBAR(:,j) = MATMUL(DG_C,A(:,j))
    END DO

	! #######
	! END CELL AVERAGING
	! #######

	! Mass filling to prevent negative subcell averages within each element
	IF(doposlimit) THEN
		CALL MFILL(rpBAR,N,nelem)
		CALL MFILL(rqBAR,N,nelem)
	END IF

    ! Reform original rp and rq vectors to send back
	! -- Note: There are 2 ghost cells that we don't update

    DO j=1,nelem
		rhop(1+(N+1)*(j-1) : (N+1)*j) = rpBAR(:,j)
		rhoq(1+(N+1)*(j-1) : (N+1)*j) = rqBAR(:,j)
    END DO


!	IF(dorhoupdate) THEN
!	ELSE	
!	! Density at physical time levels is already known ; just set values
!		rhop(:) = jcbn1d(:)
!	ENDIF

END SUBROUTINE mDGsweep

REAL(KIND=KIND(1D0)) FUNCTION B(qQuadVals,rhopQuadVals,flx,utild,wghts,nodes,k,j,nelem,N,fluxcf,Leg,dLeg) 
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem, N
	INTEGER, INTENT(IN) :: k,j ! Mode number, element number
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes,dLeg
!	REAL(KIND=DOUBLE), DIMENSION(0:N,0:nelem+1), INTENT(IN) :: Ain
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: qQuadVals,rhopQuadVals
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: utild
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: Leg

	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: flx ! Array of element edge fluxes flx(j) = F(j,j+1)
	
	! Flux correction factor. fluxcf(j)=reduction in flux coming through face (j+1/2)
	! fluxcf(0)=reduction factor for flux coming through left boundary of domain.
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: fluxcf
 
	! --- Local variables
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: HOLDER
	INTEGER :: i

!	DO i=0,N
!		HOLDER(i) = SUM(Ain(:,j)*Leg(:,i))
!	END DO

!	B = SUM(wghts(:)*utild(:,j)*dLeg(:)*HOLDER)

	B = SUM(wghts(:)*utild(:,j)*dLeg(:)*qQuadVals(:)*rhopQuadVals(:))

	IF(k .eq. 0) THEN
		B = B - fluxcf(j)*flx(j) + fluxcf(j-1)*flx(j-1)
	ELSE
		B = B - flx(j) + ((-1D0)**k)*flx(j-1)
	END IF
	
	B = (2D0*k+1D0)*B

END FUNCTION B

SUBROUTINE NUMFLUX(rhopEdgeVals,qEdgeVals,rhouEdge,N,nelem,flxrp,flxrq) 
!	USE mDGmod
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: N,nelem
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(IN) :: rhopEdgeVals,qEdgeVals
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: rhouEdge

	! -- Outputs	
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flxrp,flxrq

	! -- Local variables
	INTEGER :: j

	DO j=0,nelem
		flxrp(j) = 0.5D0*rhopEdgeVals(1,j)*(rhouEdge(j)+DABS(rhouEdge(j)))+0.5D0*rhopEdgeVals(0,j+1)*(rhouEdge(j)-DABS(rhouEdge(j)))
!		flxrq(j) = 0.5D0*qEdgeVals(1,j)*rhopEdgeVals(1,j)*(rhouEdge(j)+DABS(rhouEdge(j)))+ &
!					0.5D0*qEdgeVals(0,j+1)*rhopEdgeVals(0,j+1)*(rhouEdge(j)-DABS(rhouEdge(j)))
		flxrq(j) = 0.5D0*qEdgeVals(1,j)*(rhouEdge(j)+DABS(rhouEdge(j)))+0.5D0*qEdgeVals(0,j+1)*(rhouEdge(j)-DABS(rhouEdge(j)))
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
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: Acur,Apre
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
!		IF(Qj .lt. 0D0) THEN
!			write(*,*) 'Stage 1 Qj < 0! j=',j,'Qj=',Qj
!			write(*,*) Acur(0,j)
!		END IF

			CASE(2) 
				Qj = (dxel/dt)*(3D0*Apre(0,j) + 1D0*Acur(0,j))
!				Qj = (dxel**2/dt)*( 3D0*SUM(MATMUL(DG_C,Apre(:,j))) + 1D0*SUM(MATMUL(DG_C,Acur(:,j))) )
!		IF(Qj .lt. 0D0) THEN
!			write(*,*) 'Stage 2 Qj < 0! j=',j,'Qj=',Qj
!			write(*,*) Apre(0,j),Acur(0,j)
!		END IF

			CASE(3) 
				Qj = (dxel/(2D0*dt))*(1D0*Apre(0,j) + 2D0*Acur(0,j))
!				Qj = (dxel**2/(2D0*dt))*( 1D0*SUM(MATMUL(DG_C,Apre(:,j))) + 2D0*SUM(MATMUL(DG_C,Acur(:,j))) )
!		IF(Qj .lt. 0D0) THEN
!			write(*,*) 'Stage 3 Qj < 0! j=',j,'Qj=',Qj
!			write(*,*) Apre(0,j),Acur(0,j)
!		END IF

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
!	fluxcf = R(0:nelem)
	DO j=0,nelem
		! If flux at right edge is negative, use limiting ratio in element to the right of current one
		! (since that is where we are pulling mass from)
!		IF(flx(j) < 0D0) THEN
!			fluxcf(j) = R(j+1)
!		END IF
		fluxcf(j) = R(j) - (1D0-INT(SIGN(1D0,flx(j))))/2D0*(R(j)-R(j+1))
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

		r = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))
!		IF(r .gt. 1D0) THEN
!			write(*,*) 'WARNING REDUCTION RATIO > 1.0!!'
!		ENDIF
		rhoq(:,j) = r*rhoq(:,j) ! Reduce remaining positive masses by reduction factor

	ENDDO

END SUBROUTINE MFILL

SUBROUTINE evalExpansion(A,leg_quad,quadVals,edgeVals,N,nelem)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,N
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: leg_quad
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(INOUT) :: A

	! -- Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(OUT) :: quadVals
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(OUT) :: edgeVals
	! -- Local variables
	INTEGER :: i,k,j

	DO j=1,nelem

		DO i=0,N
			quadVals(i,j) = SUM(A(:,j)*leg_quad(:,i))
		ENDDO

		edgeVals(0,j) = SUM(A(:,j)*(/ ((-1D0)**i , i=0,N) /)) ! Expansion value at left edge of element
		edgeVals(1,j) = SUM(A(:,j)) ! Expansion value at right edge of element

	ENDDO

	! Extend edgeVals periodically
	edgeVals(:,0) = edgeVals(:,nelem)
	edgeVals(:,nelem+1) = edgeVals(:,1)

END SUBROUTINE evalExpansion
