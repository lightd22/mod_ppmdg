! --------------------------------------------------------------------
! Modal DG sweep update for DG/FV hybrid 1D Advection
! By: Devin Light
! --------------------------------------------------------------------

SUBROUTINE mDGsweep(rhoq0,rho0,rhoqIn,rhoIn,u,uedge,dxel,nelem,N,wghts,nodes,DG_C,DG_LUC,DG_L,DG_DL,IPIV,dt, & 
			        doposlimit,dorhoupdate,stage)
	USE mDGmod
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- External functions
	REAL(KIND=DOUBLE), EXTERNAL :: B ! RHS function for evolution ODE for kth expansion coefficent

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem,N ! Number of elements, highest degree of Legendre polynomial being used
	INTEGER, INTENT(IN) :: stage
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt! Element width, Finite Volume sub-cell width, time step, 
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes ! Quadrature weights and node locations
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_C, DG_LUC,DG_L,DG_DL! C matrix, used to xfer between grids, LU decomp of C
	INTEGER, DIMENSION(0:N), INTENT(IN):: IPIV ! Pivot array for RHS when using DG_LUC
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(IN) :: u ! Velocities at quadrature locations within each element at 3 time levels (tn, tn+dt,tn+dt/2)
	REAL(KIND=DOUBLE), DIMENSION(1:nelem), INTENT(IN) :: uedge ! Edge velocities at RHS of each element at 3 time levels
	LOGICAL, INTENT(IN) :: doposlimit,dorhoupdate

	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(IN) :: rho0,rhoq0 ! Solution at time level n
	! --- Outputs
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(INOUT) :: rhoqIn,rhoIn ! Output solution at next stage level

	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))) :: tmp
	! --- Local Variables
	INTEGER :: i,j,k
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqBAR0,rqBAR,rhoBAR0,rhoBAR,qBAR ! Reshaped cell averages
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: A0,AIn,AOut,R0,RIn,ROut,Q ! Reshaped DG coefficent matricies
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: utild ! Reshaped velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: uedgetild ! Periodically extended edge velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flxrq, flxrp! Array of fluxes F(j,j+1) (flx(0) is left flux at left edge of domain)
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: fcfrq,fcfrp ! Flux correction factors for positivity limiting
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqQuadVals,rhoQuadVals,qQuadVals,qQuadVals2
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1) :: rqEdgeVals,rhoEdgeVals,qEdgeVals,ones,qEdgeVals2
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: qOnes

	REAL(KIND=DOUBLE), DIMENSION(1:nelem) :: error

	REAL(KIND=DOUBLE) :: M0, M1,Mq

	ones = 1D0
	qOnes = 1D0

	! ########################################
    ! A(k,j) gives a_k(t) in the jth element for rhoq 
	! R(k,j) gives a_k(t) in the jth element for rho
    ! ########################################
		
    ! Reshape incoming values to be more convienent
    DO j=1,nelem
		rhoBAR0(:,j) = rho0(1+(N+1)*(j-1) : (N+1)*j)
		rqBAR0(:,j) = rhoq0(1+(N+1)*(j-1) : (N+1)*j)

        rqBAR(:,j) = rhoqIn(1+(N+1)*(j-1) : (N+1)*j)
		rhoBAR(:,j) = rhoIn(1+(N+1)*(j-1) : (N+1)*j)

		qBAR(:,j) = rqBAR(:,j)/rhoBAR(:,j)

        utild(:,j) = u(1+(N+1)*(j-1) : (N+1)*j)
    END DO

	uedgetild(1:nelem) = uedge(1:nelem)
	uedgetild(0) = uedgetild(nelem)

    ! For DG hybrid, values incoming assumed to be cell averages -> Project onto basis, giving coefficents a_k(t) 
	! which corresp. to series that, when averaged over the evenly spaced subcells in each element, result in the incoming values.

	CALL projectAverages(R0,DG_LUC,IPIV,rhoBAR0,N,nelem)
	CALL projectAverages(A0,DG_LUC,IPIV,rqBAR0,N,nelem)

	CALL projectAverages(AIn,DG_LUC,IPIV,rqBAR,N,nelem)
	CALL projectAverages(RIn,DG_LUC,IPIV,rhoBAR,N,nelem)

	CALL projectAverages(Q,DG_LUC,IPIV,qBAR,N,nelem)

    ! UPDATE rhoQ and rho ;

!	CALL evalExpansion(A,DG_L,rqQuadVals,rqEdgeVals,N,nelem)
	CALL evalExpansion(Q,DG_L,qQuadvals,qEdgevals,N,nelem)
	CALL evalExpansion(R0,DG_L,rhoQuadVals,rhoEdgeVals,N,nelem) ! Maybe should use R0 here? (Ideally shouldn't matter, since fwd steps for rho should yield a constant)

	CALL NUMFLUX(rhoEdgeVals,qEdgeVals,uEdgeTild,N,nelem,flxrp,flxrq) 

	! Determine flux correction factors
	fcfrp = 1d0 ! By default, do no corrections
	fcfrq = 1d0 ! By default, do no corrections
!	IF(doposlimit) THEN
!		CALL FLUXCOR(Rin,R0,flxrp,DG_C,dxel,dt,nelem,N,stage,fcfrp) ! Maybe don't need this? Shouldn't need to limit rho
!		CALL FLUXCOR(Ain,A0,flxrq,DG_C,dxel,dt,nelem,N,stage,fcfrq)
!	END IF

	! Do forward step
	DO j=1,nelem	
		DO k=0,N
			Aout(k,j) = Ain(k,j) + (dt/dxel)*B(qQuadVals(:,j),rhoQuadVals,flxrq,uTild,wghts,nodes,k,j,nelem,N,fcfrq,DG_L,DG_DL(k,:))
		END DO
	END DO

	IF(dorhoupdate) THEN
		DO j=1,nelem	
			DO k=0,N
				Rout(k,j) = Rin(k,j) + (dt/dxel)*B(qOnes(:),rhoQuadVals,flxrp,uTild,wghts,nodes,k,j,nelem,N,fcfrp,DG_L,DG_DL(k,:))	
			END DO
		END DO
	ELSE
		Rout = R0
!		SELECT CASE(stage)
!			CASE(2)
!				DO j=1,nelem
!					Aout(:,j) = 0.75D0*Ain(:,j)+0.25D0*Aout(:,j) 
!				ENDDO
!			CASE(3)
!				DO j=1,nelem
!					Aout(:,j) = (1D0/3D0)*Ain(:,j)+(2D0/3D0)*Aout(:,j) 
!				ENDDO
!			CASE DEFAULT
!				Aout = Aout
!		END SELECT
	ENDIF


    ! After time stepping, use DG_C to average the series expansion to cell-averaged values
    ! on evenly spaced grid for finite volumes
	DO j=1,nelem
		rqBAR(:,j) = MATMUL(DG_C,Aout(:,j))
		rhoBAR(:,j) = MATMUL(DG_C,Rout(:,j))
	ENDDO

	! #######
	! END CELL AVERAGING
	! #######

	! Mass filling to prevent negative subcell averages within each element
	IF(doposlimit) THEN
!		CALL MFILL(rpBAR,N,nelem)
		CALL MFILL(rqBAR,N,nelem)
!		CALL MFILL(qBAR,N,nelem)

!		rqBAR(:,:) = rpBAR(:,:)*qBAR(:,:) ! Necessary to ensure q = (rq)/r
	END IF

    ! Reform original rp and rq vectors to send back
    DO j=1,nelem
		rhoIn(1+(N+1)*(j-1) : (N+1)*j) = rhoBAR(:,j)
		rhoqIn(1+(N+1)*(j-1) : (N+1)*j) = rqBAR(:,j)
    END DO

END SUBROUTINE mDGsweep

REAL(KIND=KIND(1D0)) FUNCTION B(qQuadVals,rhopQuadVals,flx,utild,wghts,nodes,k,j,nelem,N,fluxcf,Leg,dLeg) 
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem, N
	INTEGER, INTENT(IN) :: k,j ! Mode number, element number
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes,dLeg
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

SUBROUTINE NUMFLUX(rhoEdgeVals,qEdgeVals,uEdge,N,nelem,flxrp,flxrq) 
!	USE mDGmod
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: N,nelem
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(IN) :: rhoEdgeVals,qEdgeVals
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: uEdge

	! -- Outputs	
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flxrp,flxrq

	! -- Local variables
	INTEGER :: j

	DO j=0,nelem
		flxrp(j) = 0.5D0*rhoEdgeVals(1,j)*(uEdge(j)+DABS(uEdge(j)))+0.5D0*rhoEdgeVals(0,j+1)*(uEdge(j)-DABS(uEdge(j)))
		flxrq(j) = 0.5D0*rhoEdgeVals(1,j)*qEdgeVals(1,j)*(uEdge(j)+DABS(uEdge(j)))+ & 
				   0.5D0*rhoEdgeVals(0,j+1)*qEdgeVals(0,j+1)*(uEdge(j)-DABS(uEdge(j)))

!		flxrq(j) = 0.5D0*qEdgeVals(1,j)*(rhouEdge(j)+DABS(rhouEdge(j)))+0.5D0*qEdgeVals(0,j+1)*(rhouEdge(j)-DABS(rhouEdge(j)))
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

SUBROUTINE projectAverages(A,DG_LUC,IPIV,avgs,N,nelem)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,N
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_LUC
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: avgs
	INTEGER, DIMENSION(0:N), INTENT(IN) :: IPIV
	! -- Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: A
	! -- Local variables
	INTEGER :: i,j,k
	REAL(KIND=DOUBLE), DIMENSION(1:nelem) :: hold
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: fooBAR
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: FOO_y

	fooBAR = avgs

	DO i=0,N ! Reorder RHS according to IPIV
			hold = fooBAR(i,1:nelem)
			fooBAR(i,:) = fooBAR(IPIV(i)-1,:)
			fooBAR(IPIV(i)-1,1:nelem) = hold
	ENDDO

	DO j=1,nelem
		FOO_y = 0D0
		! Solve Ly=RHS for y
		FOO_y(0) = fooBAR(0,j)
		DO k=1,N
			FOO_y(k) = fooBAR(k,j) - SUM(DG_LUC(k,0:k-1)*FOO_y(0:k-1))
		ENDDO
		! Solve Ux=y for x
		A(N,j) = (1D0/DG_LUC(N,N))*FOO_y(N)
		DO k=N-1,0,-1
			A(k,j) = (1D0/DG_LUC(k,k))*(FOO_y(k) - SUM(DG_LUC(k,k+1:N)*A(k+1:N,j)))
		ENDDO
	ENDDO

END SUBROUTINE projectAverages

SUBROUTINE averageQ(qBAR,A,R,wghts,nodes,N,nelem)
	USE mDGmod
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! Inputs
	INTEGER, INTENT(IN) :: N,nelem
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,nodes
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: A,R
	! Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(OUT) :: qBAR
	! Local Variables
	INTEGER :: i,j,k
	REAL(KIND=DOUBLE) :: dxi,xiCenter
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: rq,rho
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N) :: z,leg

	dxi = 2D0/DBLE(N+1)
	DO i=0,N
		xiCenter = dxi/2D0 + (i-1)*dxi-1D0
		z(i,0:N) = xiCenter+nodes(0:N)*dxi/2D0
	ENDDO

	DO j=1,nelem
		DO i=0,N
			CALL evalLegendre(leg,z(i,:),N,N)
			qBAR(i,j) = 0.5D0*SUM( wghts(:)*SUM(A(:,j)*leg(:,i))/SUM(R(:,j)*leg(:,i)) )
		ENDDO
	ENDDO
END SUBROUTINE averageQ
