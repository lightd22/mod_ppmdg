program execute
 !compileagain
  USE mDGmod ! Dev: Module containing stuff for DG
  USE netCDF

  implicit none

  logical :: doposlimit, domonlimit, dosellimit, dosemilagr, dopcm, &
       dowenosplit, dowenounsplit,dodghybrid,dotimesteptest
  logical :: transient, nofluxns, nofluxew, oddstep, firststep, transient2, &
       transient3
  real(kind=8) :: time, lambdamax, epslambda
  integer :: nmethod, nmethod2, nmaxcfl ! Dev: nmethod2 is used for sweeps in y-direction
  integer :: start_res,norder

  nmaxcfl=2

  dosemilagr = .false.  ! use CFL<=1 routines by default
  nofluxns = .false.
  nofluxew = .false.
  transient = .false.
  transient2 = .false.
  transient3 = .false.

  norder = 4 ! order of reconstructing polynomials for DG method

  write(*,*) '================'
  write(*,*) 'TEST #0: Uniform field in deformation flow '
  write(*,*) '================'
  start_res = 8*(norder+1)
  transient = .true.
!  call test2dweno(100,start_res,start_res,2,3,0.d0,0.d0,20,0.01D0) !1D0/(2D0*4D0-1D0)

  write(*,*) '================'
  write(*,*) 'TEST #1: Constant Diagonal Advection '
  write(*,*) '================'
  transient = .false.
!  call test2dweno(1,start_res,start_res,2,3,0.d0,0.d0,20,0.08d0) !1D0/(2D0*4D0-1D0)
!  call test2dweno(99,start_res,start_res,2,3,0.d0,0.d0,20,0.01D0) !1D0/(2D0*4D0-1D0)

  write(*,*) '================'
  write(*,*) 'TEST #2: Smooth cosbell deformation'
  write(*,*) '================'
  transient = .true.
!  call test2dweno(6,start_res,start_res,2,3,0.d0,0.d0,20,0.08D0)

  write(*,*) '================'
  write(*,*) 'TEST #3: Standard cosbell deformation'
  write(*,*) '================'
  transient = .true.
  call test2dweno(5,start_res,start_res,2,3,0.d0,0.d0,20,0.08d0)

  write(*,*) '================'
  write(*,*) 'TEST #4: Solid body rotation of cylinder'
  write(*,*) '================'
  transient = .false.
!  call test2dweno(101,start_res,start_res,2,3,0.d0,0.d0,20,0.08d0)

  write(*,*) '================'
  write(*,*) 'TEST #5: Solid body rotation of cylinder (scaled to compare to franks output)'
  write(*,*) '================'
  transient = .false.
!  call test2dweno(201,start_res,start_res,2,3,0.d0,0.d0,20,0.05d0)
   

!  nofluxns = .true.
!  nofluxew = .true.
!  transient3 = .true.
!  call test2dweno(14,50,50,2,3,0.d0,0.d0,8,1.d0) ! Each of these use no-flux conditions in both NS and EW, save for later
!  call test2dweno(15,50,50,2,3,0.d0,0.d0,8,1.d0)
!  nofluxns = .false.
!  nofluxew = .false.

contains

  subroutine test2dweno(ntest,nx0,ny0,nscale,nlev,xstr,ystr,noutput,maxcfl)
    implicit none

    integer, intent(in) :: ntest, nx0, ny0, nscale, nlev, noutput
    real (kind=8), intent(in) :: maxcfl, xstr, ystr

    real (kind=8), allocatable, dimension(:,:) :: q,q0,dqdt,u,v,uout,vout, &
         utilde,vtilde,jcbn, rho, rhoprime, rhoq, u2,v2,u2tilde,v2tilde
    real (kind=8), allocatable, dimension(:) :: dx, qi, qiface, x, xf, upos, &
         dy, qj, qjface, y, yf, vpos
    real (kind=8), allocatable, dimension(:,:) :: xlambda, xmonlimit, ylambda, ymonlimit
   
    real(kind=8), dimension(nlev) :: e1, e2, ei
    integer :: npad, tmp_method(27)
    integer :: nx, ny, nstep, i ,j , n, p, ierr, nout, imethod, ncheck, n1, n2, nn

    real (kind=8) :: dt, tfinal, pi, &
         tmp_time, tmp_umax, tmp_vmax, cnvg1, cnvg2, cnvgi, &
         tmp_qmin, tmp_qmax,calculatedMu,dxm,dym
    REAL(KIND=8), DIMENSION(1:2) :: domainCenter, xEdge, yEdge
    character(len=40) :: cdf_out 
    real (kind=8), external :: tfcn
    real(kind=4), dimension(2) :: tstart,tend
    real(kind=4) :: t0,tf
	real(kind=4) :: talpha,tbeta
    character(len=8) :: outdir

	! DG variables
	INTEGER :: nex,ney,quadOrder
	REAL(KIND=8) :: dxel,dyel
    REAL(KIND=8), allocatable, dimension(:) :: DG_xec,DG_yec, DG_nodes,DG_wghts,DG_x,DG_y,DG_FOO
	INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV
    REAL(KIND=8), allocatable, dimension(:,:) :: DGu0,DGv0, DGuedge0,DGvedge0,DG_C,DG_LUC,DG_L,DG_DL


    pi = DACOS(-1D0)

    if(nlev.lt.1) STOP 'nlev should be at least 1 in test2dweno'

    n2 = 2
	go to 300
    tmp_method(1) = 1  ! PPM no limiting
    tmp_method(2) = 2  ! PPM, positive-definite using FCT
    tmp_method(3) = 6  ! PPM, positive selective limiting using FCT
    tmp_method(4) = 37 ! PPM, selective limiting using PMOD
    tmp_method(5) = 38 ! PPM, selective and positivity limiting using PMOD
	tmp_method(6) = 8 ! PPM, positive-defininte using PMOD
	tmp_method(7) = 98 ! DG, avg, no limiting
	tmp_method(8) = 99 ! DG, avg, mass filling (positive-definite)
    tmp_method(9) = 100 ! PPM/DG Hybrid, no limiting
	tmp_method(10) = 101 ! PPM/DG Hybrid, positive-definite limiting using FCT
	tmp_method(11) = 102 ! PPM/DG Hybrid, positive-definite limiting using PMOD

	300 continue
	tmp_method(1) = 98
	tmp_method(2) = 99
	tmp_method(3) = 1
	tmp_method(4) = 2



    do nn = 1,n2
       imethod = tmp_method(nn)

      
       dosellimit = .false.
       domonlimit = .false.
       doposlimit = .false.
       dopcm = .false.  ! use PPM by default
       dowenosplit = .false.  ! use PPM by default
       dowenounsplit = .false.  ! use PPM by default
	   dodghybrid = .false.

       !-------------------
       ! Dev: nmethod2 will be used for y-sweeps. For most methods it will be same as nmethod, but it will be different for DG/PPM
       !-------------------

       select case(imethod)
       case(1)
          write(*,*) 'PPM FCT, No limiting'
          outdir = 'pfctnon/'
          nmethod = 1
          nmethod2 = 1
       case(2)
          write(*,*) 'PPM FCT, Positive'
          doposlimit = .true.
          outdir = 'pfctpos/'
          nmethod = 2
          nmethod2 = 2
       case(3)
          write(*,*) 'PPM FCT, Monotonic'
          domonlimit = .true.
          outdir = 'pfctmon/'
          nmethod = 3
          nmethod2 = 3
       case(4)
          write(*,*) 'PPM FCT, Monotonic, Positive'
          domonlimit = .true.
          doposlimit = .true.
          outdir = 'pfctpmn/'
          nmethod = 4
          nmethod2 = 4
       case(5)
          write(*,*) 'PPM FCT, Selective'
          dosellimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'pfctsel/'
          nmethod = 5
          nmethod2 = 5
       case(6)
          write(*,*) 'PPM FCT, Selective, Positive'
          dosellimit = .true. 
          doposlimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'pfctpse/'
          nmethod = 6
          nmethod2 = 6
       case(7)
          write(*,*) 'PPM PMOD, No limiting'
          outdir = 'ppmdnon/'
          nmethod = 21
          nmethod2 = 21
       case(8)
          write(*,*) 'PPM PMOD, Positive'
          doposlimit = .true.
          outdir = 'ppmdpos/'
          nmethod = 22
          nmethod2 = 22
       case(9)
          write(*,*) 'PPM PMOD, Monotonic'
          domonlimit = .true.
          outdir = 'ppmdmon/'
          nmethod = 23
          nmethod2 = 23
       case(10)
          write(*,*) 'PPM PMOD, Monotonic, Positive'
          domonlimit = .true.
          doposlimit = .true.
          outdir = 'ppmdpmn/'
          nmethod = 24
          nmethod2 = 24
       case(11)
          write(*,*) 'PPM PMOD, Selective'
          dosellimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'ppmdsel/'
          nmethod = 25
          nmethod2 = 25
       case(12)
          write(*,*) 'PPM PMOD, Selective, Positive'
          dosellimit = .true. 
          doposlimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'ppmdpse/'
          nmethod = 26
          nmethod2 = 26
       case(37)
          write(*,*) 'PPM, Selective, Clean'
          dosellimit = .true. 
          doposlimit = .false.
          domonlimit = .false.
          outdir = 'ppmsecl/'
          nmethod = 81
          nmethod2 = 81
       case(38)
          write(*,*) 'PPM, Selective, Clean, Positive'
          dosellimit = .true. 
          doposlimit = .true.
          domonlimit = .false.
          outdir = 'ppmpscl/'
          nmethod = 82
          nmethod2 = 82
		case(98)
		  write(*,*) 'DG, averages, no limiting'
		  write(*,*) 'WARNING: Should only be used with periodic BCs'
		  dodghybrid = .true.
		  doposlimit = .false.
		  outdir = 'dgnolim/'
		  nmethod = 99
		  nmethod2 = 99
          quadOrder = norder!CEILING( (3D0/2D0)*norder )
		  write(*,FMT='(A5,i1,A15,i1)') ' N = ',norder,', Quad Order = ',quadOrder
		case(99)
		  write(*,*) 'PD-DG, averages, element mass redist'
		  write(*,*) 'WARNING: Should only be used with periodic BCs'
		  dodghybrid = .true.
		  doposlimit = .true.
		  outdir = 'dgmfill/'
		  nmethod = 99
		  nmethod2 = 99
		  write(*,*) 'N = ',norder
        case(100)
          write(*,*) 'PPM/DG Hybrid, no limiting'
          write(*,*) 'WARNING: Should only be used with periodic BCs'
          outdir = 'ppmdghy/'
		  dodghybrid = .true.
          nmethod = 99
          nmethod2 = 1
          write(*,*) 'Number of GLL nodes=',norder+1
		case(101)
		  write(*,*) 'PPM/DG Hybrid, FCT, Positive, N=5'
		  write(*,*) 'WARNING: Should only be used with periodic BCs'
		  outdir = 'ppmdgfc/'
		  dodghybrid = .true.
          doposlimit = .true.
		  nmethod = 99
		  nmethod2 = 2
		  write(*,*) 'Number of GLL nodes=',norder+1
		case(102)
		  write(*,*) 'PPM/DG Hybrid, PMOD, Positive N=5'
		  write(*,*) 'WARNING: Should only be used with periodic BCs'
		  outdir = 'ppmdgpm/'
		  dodghybrid = .true.
		  doposlimit = .true.
		  nmethod = 99
		  nmethod2 = 22
		  write(*,*) 'Number of GLL nodes=',norder+1
       end select

       if(dopcm.and.domonlimit) then
          npad = 6
       elseif(dopcm.and.dosellimit.and..not.doposlimit) then
          npad = 3
       elseif(dosellimit.or.doposlimit) then
          npad = 4
       else
          npad = 3
       end if
       write(*,*) 'npad = ', npad

    ! For DG: Compute gaussian nodes in [-1,1] + weights for quadratures
    ! (only needs to be done once)

    allocate(DG_nodes(0:norder), DG_wghts(0:norder),DG_L(0:norder,0:norder),DG_DL(0:norder,0:norder),STAT=ierr)
    DG_nodes = 0D0
    DG_wghts = 0D0
    IF(nmethod .eq. 99 .OR. nmethod2 .eq. 99) THEN
        CALL quad_nodes(norder+1,DG_nodes)
        CALL quad_weights(norder+1,DG_nodes, DG_wghts)
		! Fill in array of Legendre polynomials evaluated at quad nodes + Leg. derivative at quad nodes
		DO i=0,norder
			DO j=0,norder
				DG_L(i,j) = legendre(DG_nodes(j),i)
				DG_DL(i,j) = dlegendre(DG_nodes(j),i)
			ENDDO
		ENDDO
    ENDIF

    xEdge(1) = 0D0
    xEdge(2) = 1D0
    yEdge(:) = xEdge(:)
    domainCenter(1) = SUM(xEdge)/SIZE(xEdge)
    domainCenter(2) = SUM(yEdge)/SIZE(yEdge)

    write(*,*) 'Domain is: [',xEdge(1),',',xEdge(2),'].'
    write(*,*) 'Warning: Not all tests have been implemented for non-unit square domain!'
 
    do p = 1,nlev
       
       t0 = etime(tstart)

       nx = nx0*nscale**(p-1)
       ny = ny0*nscale**(p-1)
                                        !Dev: When using DG, each element in x-direction will have norder+1 quad nodes.
       nex = INT((nx)/(norder+1))  	    !For the matrix C (used in exchanging between DG/FV) to be nonsingular, nx must be 
	   ney = INT((ny)/(norder+1))		!a multiple of norder+1

       dxel = (xEdge(2)-xEdge(1))/nex           		
	   dyel = (yEdge(2)-yEdge(1))/ney

       allocate(q(1-npad:nx+npad,1-npad:ny+npad), &
            q0(1:nx,1:ny),dqdt(1:nx,1:ny),jcbn(1:nx,1:ny), &
            u(0:nx,1:ny),utilde(0:nx,1:ny),dx(1:nx),x(1:nx),xf(0:nx), &
            v(1:nx,0:ny),vtilde(1:nx,0:ny),dy(1:ny),y(1:ny),yf(0:ny), &
            u2(0:nx,1:ny),u2tilde(0:nx,1:ny),v2(1:nx,0:ny),v2tilde(1:nx,0:ny),&
            uout(0:nx,1:ny),vout(1:nx,0:ny), &
            rho(1:nx,1:ny),rhoprime(1:nx,1:ny),rhoq(1:nx,1:ny), &
            xlambda(0:nx,1:ny),xmonlimit(0:nx,1:ny), &
            ylambda(1:nx,0:ny),ymonlimit(1:nx,0:ny), &
            STAT=ierr)
       allocate(DGu0(1:nx,1:ny),DGv0(1:nx,1:ny), DG_x(1:nx),DG_y(1:ny),DGuedge0(1:nex,1:ny),DGvedge0(1:nx,1:ney) &
				,STAT=ierr)
       allocate(DG_C(0:norder,0:norder),DG_LUC(0:norder,0:norder),IPIV(0:norder),DG_FOO(0:norder),&
				DG_xec(1:nex),DG_yec(1:ney),STAT=ierr)

       
       DGu0(:,:) = 0.D0
	   DGv0(:,:) = 0.D0
       
       ! set up x grid
    
       !----------------------
       ! Set up x-grid for PPM
       !----------------------

!        do i = 1,nx
!           dx(i) = 1.d0 + xstr*COS(2.d0*pi*DBLE(i)/DBLE(nx))
!        end do
!        dx = dx/SUM(dx) ! normalize for unit length
        dx = (xEdge(2)-xEdge(1))/nx
        xf(0) = xEdge(1)
        do i = 1,nx
              xf(i) = xf(i-1) + dx(i)
              x(i) = 0.5d0*(xf(i-1)+xf(i))
        end do


       !--------------------------
       ! Dev: Set up computational DG x-grid
       !--------------------------

       if(nmethod.eq.99) then
         DG_xec(1) = xEdge(1)+0.5d0*dxel
         do i = 2,nex
            DG_xec(i)=DG_xec(i-1)+dxel
         end do

        ! Transfer nodes to physical domain
         do i = 1,nex
            DG_x(1+(i-1)*(norder+1):i*(norder+1)) = DG_xec(i)+0.5d0*dxel*DG_nodes(0:norder)
         end do

        ! Fill in values in C-matrix
        ! C is used in transferring information between DG/PPM grids
	   IF(dodghybrid) THEN
		CALL Cmat_FILL(norder,DG_nodes,DG_wghts,dx(1),dxel,DG_C) ! Note that C assumes an evenly spaced FV sub-grid
		! Compute LU decomposition of C, stored in DG_LUC
		DG_LUC = DG_C
		DG_FOO = 0D0
		CALL DGESV(norder+1,1,DG_LUC,norder+1,IPIV,DG_FOO,norder+1,ierr)
		END IF
       end if

	   ! Set up computational DG y-grid
       if(nmethod2.eq.99) then
         DG_yec(1) = yEdge(1)+0.5d0*dyel
         do i = 2,ney
            DG_yec(i)=DG_yec(i-1)+dyel
         end do

        ! Transfer nodes to physical domain
         do i = 1,ney
            DG_y(1+(i-1)*(norder+1):i*(norder+1)) = DG_yec(i)+0.5d0*dyel*DG_nodes(0:norder)
         end do

		end if	

       ! set up y grid
!       do i = 1,ny
!          dy(i) = 1.d0 + ystr*COS(2.d0*pi*DBLE(i)/DBLE(ny))
!       end do
!       dy = dy/SUM(dy) ! normalize for unit length
       dy = (yEdge(2)-yEdge(1))/ny
       yf(0) = yEdge(1)
       do i = 1,ny
          yf(i) = yf(i-1) + dy(i)
          y(i) = 0.5d0*(yf(i-1)+yf(i))
       end do

       if((ntest.eq.7).or.(ntest.eq.8)) then
          dx = 2.d0*dx
          x = 2.d0*x
          xf = 2.d0*xf
          dy = 0.5d0*dy
          y = 0.5d0*y
          yf = 0.5d0*yf
       end if

       ! set up scalar and velocity fields
       xlambda = 0.d0
       xmonlimit = 0.d0
       ylambda = 0.d0
       ymonlimit = 0.d0
       q = 0.d0
       dqdt = 0.d0

	   ! Initialize q, velocity fields for PPM and output directory
       call init2d(ntest,nx,ny,q0,u,v,u2,v2,x,xf,y,yf,xEdge,yEdge,tfinal,cdf_out)

	   ! Initialize velocity fields at DG grid and element edges
       IF(nmethod .eq. 99 .or. nmethod2 .eq. 99) THEN 
			CALL DGinit2d(ntest,nex,ney,DG_x,DG_y,DG_xec,DG_yec,xf,yf,nx,ny,DGu0,DGuedge0,DGv0,DGvedge0,xEdge,yEdge)
       END IF

       q(1:nx,1:ny) = q0

       if (MAX(ABS(MAXVAL(dx)-MINVAL(dx)), &
            ABS(MAXVAL(dy)-MINVAL(dy))).gt.10.d0*EPSILON(dx)) then
          cdf_out = 'stretch_' // cdf_out
       end if

       cdf_out = outdir // cdf_out

       ! set up time step
       time = 0.d0
       
       if(transient3) then
          ncheck = 32
          nstep = 0
          tmp_umax = 0.d0
          tmp_vmax = 0.d0
          do i = 1,ncheck
             tmp_time = 2.d0*(dble(i)/dble(ncheck))*5.d0
             tmp_umax = MAX(tmp_umax,MAXVAL(ABS(u+u2*tfcn(tmp_time))))
             tmp_vmax = MAX(tmp_vmax,MAXVAL(ABS(v+v2*tfcn(tmp_time))))
          end do
       else
          tmp_umax = MAXVAL(ABS(u+u2))
          tmp_vmax = MAXVAL(ABS(v+v2))
       end if

       IF(dodghybrid) THEN
             dxm = dxel
             dym = dyel
         ELSE
             dxm = minval(dx)
             dym = minval(dy)
       ENDIF !dodghybrid

       if(p.eq.1) then
          if (noutput.eq.-1) then
			if(dodghybrid) then
			 nstep = ceiling(tfinal* &
                  MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl)
			else
             nstep = ceiling(tfinal* &
                  MAX(tmp_umax/MINVAL(dx),tmp_vmax/MINVAL(dy))/maxcfl)
			endif

             nout = nstep

          else
			if(dodghybrid) then
             nstep = noutput*ceiling(tfinal* &
                  MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl &
                  /DBLE(noutput))
			else
             nstep = noutput*ceiling(tfinal* &
                  MAX(tmp_umax/MINVAL(dx),tmp_vmax/MINVAL(dy))/maxcfl &
                  /DBLE(noutput))
			endif
             nout = noutput
          end if
       else
			if(dodghybrid) then
             nstep = nout*ceiling(tfinal* &
                  MAX(tmp_umax/dxel,tmp_vmax/dyel)/maxcfl &
                  /DBLE(noutput))
			else
          nstep = nout*ceiling(tfinal* &
               MAX(tmp_umax/MINVAL(dx),tmp_vmax/MINVAL(dy))/maxcfl &
               /DBLE(nout))
			endif
       end if

       if(dowenounsplit) nstep = nstep*2 ! USE HALF AS BIG CFL FOR UNSPLIT WENO
       dt = tfinal/dble(nstep)
       calculatedMu = MAX(tmp_umax/dxm,tmp_vmax/dym)*dt


       if (p.eq.1) call output2d(q(1:nx,1:ny),u,v,&
            xlambda,xmonlimit,ylambda,ymonlimit,nx,ny,x,xf,y,yf, &
            tfinal,-1,cdf_out,nout,calculatedMu,DG_nodes,DG_wghts,norder)
       
       call output2d(q(1:nx,1:ny),u,v,&
            xlambda,xmonlimit,ylambda,ymonlimit,nx,ny,x,xf,y,yf,time,0,cdf_out,p,calculatedMu,DG_nodes,DG_wghts,norder)

       !transform velocity into computational space, compute jacobian
       do j = 1,ny
          jcbn(:,j) = dx*dy(j)
          utilde(:,j) = u(:,j)*dy(j)                    ! Dev: Scales so that utilde = u*dy
          u2tilde(:,j) = u2(:,j)*dy(j)
       end do
       do i = 1,nx
          vtilde(i,:) = v(i,:)*dx(i)                    ! Dev: Scales so that vtilde = v*dx
          v2tilde(i,:) = v2(i,:)*dx(i)
       end do


       rho(1:nx,1:ny) = jcbn(1:nx,1:ny)                 ! Dev: Scales so that rho = rho*dx*dy ; 
       rhoq(1:nx,1:ny) = rho(1:nx,1:ny)*q(1:nx,1:ny)    ! (note that physical rho == 1, so rho = jac = dx*dy = rho*dx*dy)
       rhoprime(1:nx,1:ny) = rho(1:nx,1:ny)             ! rhov = rho*v*dx ; rhou = rho*u*dy

       tmp_qmin = MINVAL(q(1:nx,1:ny))
       tmp_qmax = MAXVAL(q(1:nx,1:ny))

       firststep = .true.
       do n = 1,nstep    
      

          if(mod(n,2).eq.1) then
             oddstep = .true.
          else
             oddstep = .false.
          end if
          rho(1:nx,1:ny) = jcbn(1:nx,1:ny)
          rhoq(1:nx,1:ny) = rho(1:nx,1:ny)*q(1:nx,1:ny) ! DEVIN: rhoq = rho*q
          rhoprime(1:nx,1:ny) = rho(1:nx,1:ny)          ! rhoprime = rho to begin

!!$          rho(1:nx,1:ny) = rhoprime(1:nx,1:ny)
!!$          if(dowenounsplit) then
!!$             call weno_rk3_2d(q,dqdt,utilde,vtilde,u2tilde,v2tilde, &
!!$                  nx,ny,npad,dt,jcbn)
!!$             xlambda = 0.d0
!!$             xmonlimit = 0
!!$             ylambda = 0.d0
!!$             ymonlimit = 0
!!$          else

             call skamstep_2d(q,dqdt,utilde,vtilde,u2tilde,v2tilde, &
                  rho,rhoq,rhoprime,nx,ny,npad,dt,jcbn,&
                  xlambda,xmonlimit,ylambda,ymonlimit,& 
                  DGu0,DGuedge0,DGv0,DGvedge0,nex,ney,dxel,dyel,norder,DG_wghts,DG_C,DG_LUC,IPIV,DG_L,DG_DL)

!!$          end if

          time = time + dt

          if ((mod(n,nstep/nout).eq.0).OR.(n.eq.nstep)) then
             if(transient) then
                uout = u*tfcn(time)
                vout = v*tfcn(time)
             elseif(transient2) then
                uout = u*tfcn(time)**2 + u2*(1.d0-tfcn(time)**2)
                vout = v*tfcn(time)**2 + v2*(1.d0-tfcn(time)**2)
             elseif(transient3) then
                tmp_time = 2.d0*time
                uout = u + u2*tfcn(tmp_time)
                vout = v + v2*tfcn(tmp_time)
             else
                uout = u
                vout = v
             end if
             call output2d(q(1:nx,1:ny),uout,vout,&
            xlambda,xmonlimit,ylambda,ymonlimit,nx,ny,x,xf,y,yf, &
                  time,1,cdf_out,p,calculatedMu,DG_nodes,DG_wghts,norder)
          end if

!          write(*,*) 'rho min, max = ', MINVAL(rho), MAXVAL(rho)
!          write(*,*) 'rhoprime min, max = ', MINVAL(rhoprime), MAXVAL(rhoprime)
!          write(*,*) 'rhoq min, max = ', MINVAL(rhoq), MAXVAL(rhoq)
!          write(*,*) 'q min, max = ', MINVAL(q(1:nx,1:ny)), MAXVAL(q(1:nx,1:ny))
!          write(*,*)

!!$          if(n.eq.3)     stop ' in back_traj'


          firststep = .false.

          ! KEEP TRACK OF MIN/MAX OVER ALL STEPS
          tmp_qmin = MIN(tmp_qmin,MINVAL(q(1:nx,1:ny)))
          tmp_qmax = MAX(tmp_qmax,MAXVAL(q(1:nx,1:ny)))

!		write(*,*) SUM(q(1:nx,1:ny)-q0)/DBLE(nx*ny)
!		q0 = q(1:nx,1:ny)
		

       end do
 
!       write(*,*) 'NOTE: SUM(ABS(q))=',SUM(ABS(q(1:nx,1:ny)))
!       write(*,*) ''
!       write(*,*) 'NOTE: q(1,:) =',q(1,:)
       e1(p) = SUM(ABS(q(1:nx,1:ny)-q0))/DBLE(nx)/DBLE(ny)
       e2(p) = SQRT(SUM((q(1:nx,1:ny)-q0)**2)/DBLE(nx)/DBLE(ny))
       ei(p) = MAXVAL(ABS(q(1:nx,1:ny)-q0))
       tf = etime(tend) - t0
       if (p.eq.1) then
          write(UNIT=6,FMT='(A119)') &
'  nx    ny       E1          E2         Einf   convergence rate  overshoot  undershoot   cons        cputime time step'
          cnvg1 = 0.d0
          cnvg2 = 0.d0
          cnvgi = 0.d0
       else
          cnvg1 = -log(e1(p)/e1(p-1))/log(dble(nscale))
          cnvg2 = -log(e2(p)/e2(p-1))/log(dble(nscale))
          cnvgi = -log(ei(p)/ei(p-1))/log(dble(nscale))
       end if
       write(*,990) nx, ny, e1(p), e2(p), ei(p), &
            cnvg1, cnvg2, cnvgi, &
            tmp_qmax-MAXVAL(q0), &
            MINVAL(q0)-tmp_qmin, &
            SUM(q(1:nx,1:ny)-q0)/DBLE(nx*ny), tf, nstep



       if (p.eq.nlev) &
            call output2d(q(1:nx,1:ny),u,v,&
            xlambda,xmonlimit,ylambda,ymonlimit,nx,ny,x,xf,y,yf,time,2,cdf_out,1,calculatedMu,DG_nodes,DG_wghts,norder)

       deallocate(q,q0,dqdt,u,v,utilde,vtilde,u2,v2,u2tilde,v2tilde, &
            uout,vout,jcbn,dx,x,xf,dy,y,yf,rho,rhoprime,rhoq,&
            xlambda,xmonlimit,ylambda,ymonlimit,STAT=ierr)
       deallocate(DGu0,DGuedge0,DGv0,DGvedge0,DG_x,DG_y, DG_C,DG_LUC, DG_FOO,DG_xec,DG_yec,STAT=ierr)

    end do
    deallocate(DG_nodes,DG_wghts, DG_L,DG_DL, STAT=ierr)
 end do

990    format(2i6,3e12.4,3f5.2,3e12.4,f8.2,i8)
991    format(2i6,3e12.4,3f5.2,2e12.4)
992    format(12f8.4)
  end subroutine test2dweno

  subroutine init2d(ntest,nx,ny,q,u,v,u2,v2,x,xf,y,yf,xEdge,yEdge,tfinal,cdf_out)
    implicit none
    integer, intent(in) :: ntest,nx,ny
    real(kind=8), dimension(1:nx,1:ny), intent(out) :: q
    real(kind=8), dimension(0:nx,1:ny), intent(out) :: u, u2
    real(kind=8), dimension(1:nx,0:ny), intent(out) :: v, v2
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), dimension(1:ny), intent(in) :: y
    real(kind=8), dimension(0:nx), intent(in) :: xf
    real(kind=8), dimension(0:ny), intent(in) :: yf
    REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: xEdge,yEdge
    real(kind=8), intent(out) :: tfinal
    character(len=40), intent(out) :: cdf_out 

    integer :: i, j, m, n
    real(kind=8) :: pi, tmpr, rmax, rmin, tmpx, tmpy, tmpth, xc,yc,xWidth,yWidth
    REAL(KIND=8),DIMENSION(1:2) :: domainCenter
    real(kind=8), dimension(1:nx,1:ny) :: r
    real(kind=8), dimension(0:nx,0:ny) :: psi, psi2

    pi = DACOS(-1D0)

    xWidth = xEdge(2)-xEdge(1)
    yWidth = yEdge(2)-yEdge(1)

    domainCenter(1) = xEdge(1)+xWidth/2D0
    domainCenter(2) = yEdge(1)+yWidth/2D0

    psi(:,:) = 0.d0
    psi2(:,:) = 0.d0

    select case(ntest)
    case (1:2)
       tfinal = 1.d0
       do j = 0,ny
          do i = 0,nx
             psi(i,j) = - xf(i) + yf(j) ! uniform flow, u=v=1
          end do
       end do
    case (3:4,10,101,201)
       tfinal = 1.d0
       do j = 0,ny
          do i = 0,nx
             ! solid body rotation
             tmpr = sqrt((xf(i)-domainCenter(1))**2 + (yf(j)-domainCenter(2))**2)
             psi(i,j) = pi*tmpr**2
          end do
       end do
    case(5:6,9,15,100)
       tfinal = 5.d0
       do j = 0,ny
          do i = 0,nx
             ! deformation flow from Leveque (1994) and Durran (1999)
             psi(i,j) = (1/pi)*sin(pi*xf(i))**2 * sin(pi*yf(j))**2
          end do
       end do
    case(7)
       tfinal = 5.d0
       do j = 0,ny
          do i = 0,nx
             ! first inversion test flow, modeled after the
             !   deformation flow from Leveque (1994) and Durran (1999)
             psi(i,j) = &
                  (0.05d0/pi)*sin(pi*(xf(i)-0.5d0))*sin(2.d0*pi*yf(j))
          end do
       end do
    case(8)
       tfinal = 5.d0
       do j = 0,ny
          do i = 0,nx
             ! second inversion test flow, modeled after the
             !   deformation flow from Leveque (1994) and Durran (1999)
             psi(i,j) = &
                  (0.125d0/pi)*sin(pi*(xf(i)-1.d0))*sin(4.d0*pi*yf(j)) &
                  + (0.03125d0/pi)*sin(2.d0*pi*xf(i))*sin(8.d0*pi*yf(j))
          end do
       end do
    case(11:12)
       tfinal = 5.d0
       do j = 0,ny
          do i = 0,nx
             ! second inversion test flow, modeled after the
             !   deformation flow from Leveque (1994) and Durran (1999)
             tmpr = sqrt((xf(i)-0.5d0)**2 + (yf(j)-0.5d0)**2)
             psi2(i,j) = (1.6d0*pi) &
                  *(log(1.d0 + 16.d0*tmpr**2)/96.d0 &
                  - log(1.d0 - 16.d0*tmpr**2 + 256.d0*tmpr**4)/192.d0 &
                  + sqrt(3.d0)*atan((32.d0*tmpr**2 - 1.d0)/sqrt(3.d0))/96.d0 )
             psi(i,j) = (1.6d0*pi)*tmpr**2/2.d0 - psi2(i,j)
          end do
       end do

!    case(14:17)
!
!       tfinal = 5.d0
!       do j = 0,ny
!          do i = 0,nx
!             tmpr = sqrt((xf(i)-0.5d0)**2 + (yf(j)-0.5d0)**2)
!             psi(i,j) = 0.4d0*pi*tmpr**2
!             psi2(i,j) = 0.8d0*pi &
!                  *(0.5d0*tmpr**2 &
!                    + (1.d0/96.d0)*log(1.d0 - 16.d0*tmpr**2 + 256.d0*tmpr**4) &
!                    - (1.d0/48.d0)*log(1.d0 + 16.d0*tmpr**2) &
!                    - sqrt(3.d0)/48.d0*atan((-1.d0 + 32.d0*tmpr**2)/sqrt(3.d0)))
!          end do
!       end do
    case(99)
       tfinal = 1.d0
       do j = 0,ny
          do i = 0,nx
             psi(i,j) =  - xf(i) + yf(j) !yf(j) ! uniform flow: u=1; v=1
          end do
       end do

    
    end select

    ! compute u velocity from streamfunction
    do j = 1,ny
       u(:,j) = (psi(:,j) - psi(:,j-1))/(yf(j) - yf(j-1))
       u2(:,j) = (psi2(:,j) - psi2(:,j-1))/(yf(j) - yf(j-1))
    end do

    ! compute v velocity from streamfunction
    do i = 1,nx
       v(i,:) = - (psi(i,:) - psi(i-1,:))/(xf(i) - xf(i-1))
       v2(i,:) = - (psi2(i,:) - psi2(i-1,:))/(xf(i) - xf(i-1))
    end do

!!!!! now, compute initial scalar field

    select case(ntest)
    case (1) ! sine wave advection by uniform velocity field
       cdf_out =  'weno2d_adv_sine.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_adv_sine.nc'
       do j = 1,ny
          q(:,j) = sin(2.d0*pi*x)*sin(2.d0*pi*y(j))
       end do
    case (2) ! sine^4 bump advection by uniform velocity field
       tfinal = 5.d0
       cdf_out =  'weno2d_adv_sine4.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_adv_sine4.nc'
       do j = 1,ny
          q(:,j) = (sin(pi*x)*sin(pi*y(j)))**4
       end do
    case (3) ! solid body rotation of a cone
       cdf_out =  'weno2d_rot_coneblock.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_rot_coneblock.nc'
       do j = 1,ny
          r(:,j) = sqrt((x-0.275d0)**2 + (y(j)-0.5d0)**2)/0.175d0
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = 1.d0-r
       end where
       do j = 1,ny
          r(:,j) = MIN((x-0.55d0)*(0.8d0-x),(y(j)-0.375d0)*(0.625d0-y(j)))
       end do
       where (r.gt.0.d0)
          q = 1.d0
       end where
    case (4) ! solid body rotation of a block
       cdf_out =  'weno2d_rot_block.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_rot_block.nc'
       do j = 1,ny
          r(:,j) = MIN((x-0.55d0)*(0.8d0-x),(y(j)-0.375d0)*(0.625d0-y(j)))
       end do
       q = 0.d0
       where (r.gt.0.d0)
          q = 1.d0
       end where
    case (5) ! deformation/return flow applied to cosine bell
       cdf_out =  'weno2d_def_cosinebell.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_def_cosinebell.nc'
       ! x = 0.25d0 , y = 0.25d0 originally

       xc = xEdge(1)+xWidth/4D0
       yc = yEdge(1)+yWidth/4D0

       do j = 1,ny
          r(:,j) = 4.d0*sqrt((x-xc)**2 + (y(j)-yc)**2)
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = (0.5d0*(1.d0 + cos(pi*r)))
       end where
    case (6) ! deformation/return flow applied to smoother cosine bell
       cdf_out =  'weno2d_def_smooth_cosinebell.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_def_smooth_cosinebell.nc'
        ! x = 0.4d0 , y = 0.4d0 originally

       xc = xEdge(1)+0.4D0*xWidth
       yc = yEdge(1)+0.4D0*yWidth
       do j = 1,ny
          r(:,j) = 3.d0*sqrt((x-xc)**2 + (y(j)-yc)**2)
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = (0.5d0*(1.d0 + cos(pi*r)))**3
       end where
    case (9) ! deformation/return flow applied to uniform scalar distribution
       cdf_out =  'weno2d_def_consistency_test.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_def_consistency_test.nc'
       q = 1.d0 
    case (10) ! rotating flow applied to cosine bell
       cdf_out =  'weno2d_rot_cosinebell.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_rot_cosinebell.nc'
       tfinal = 5.d0
       do j = 1,ny
          r(:,j) = sqrt((x-0.275d0)**2 + (y(j)-0.5d0)**2)/0.15d0
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = 0.5d0*(1.d0 + cos(pi*r))
       end where
    case (11) ! rotating flow applied to cosine bell
       cdf_out =  'weno2d_fdef_cosinebell.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_fdef_cosinebell.nc'
       tfinal = 5.d0
       do j = 1,ny
          r(:,j) = sqrt((x-0.3d0)**2 + (y(j)-0.5d0)**2)/0.2d0
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = 0.5d0*(1.d0 + cos(pi*r))
       end where
    case (12:14) ! rotating flow applied to smooth cosine bell
       cdf_out =  'weno2d_fdef_smooth_cosinebell.nc'
       if(dosemilagr) then
          cdf_out = 'sl_'// 'weno2d_fdef_smooth_cosbell.nc'
       end if

       if(ntest.eq.14) then
          cdf_out =  'weno2d_fdefvel_smooth_cosbell.nc'
          if(dosemilagr) cdf_out =  'sl_weno2d_fdefvel_smooth_bell.nc'
       end if
       tfinal = 5.d0
       do j = 1,ny
          r(:,j) = sqrt((x-0.3d0)**2 + (y(j)-0.5d0)**2)/0.2d0
       end do
       q = 0.d0 !1.d0
       where (r.lt.1.d0)
          q = q + (0.5d0*(1.d0 + cos(pi*r)))**2
       end where
    case (15) ! rotating flow applied to square wave
       cdf_out =  'weno2d_fdef_sqwave.nc'
       if(dosemilagr) then
          cdf_out = 'sl_'// 'weno2d_fdef_sqwave.nc'
       end if

       tfinal = 5.d0
       do j = 1,ny
          r(:,j) = max(abs(x-0.3d0),abs(y(j)-0.5d0))/0.15d0
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = 1.d0
       end where
    case (16) ! rotating flow applied to smooth cosine bell with background
       cdf_out =  'weno2d_fdefvel_smback_cosbell.nc'
       if(dosemilagr) then
          cdf_out = 'sl_'// 'weno2d_fdefvel_smback_bell.nc'
       end if
       tfinal = 5.d0
       do j = 1,ny
          r(:,j) = sqrt((x-0.3d0)**2 + (y(j)-0.5d0)**2)/0.2d0
       end do
       q = 1.d0
       where (r.lt.1.d0)
          q = q + (0.5d0*(1.d0 + cos(pi*r)))**2
       end where
    case (17) ! rotating flow applied to square wave
       cdf_out =  'weno2d_3xfdef_sqwave.nc'
       if(dosemilagr) then
          cdf_out = 'sl_'// 'weno2d_3xfdef_sqwave.nc'
       end if

       do j = 1,ny
          r(:,j) = max(abs(x-0.3d0),abs(y(j)-0.5d0))/0.15d0
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = 1.d0
       end where

    case(99) !thin cosinebell advection with horizontal wind
       cdf_out =  'weno2d_hadv_cosinebell.nc'
       if(dosemilagr) cdf_out =  'sl_'// 'weno2d_hadv_cosinebell.nc'
       do j = 1,ny
          r(:,j) = 4.d0*sqrt(10D0*(x-0.25d0)**2 + (y(j)-0.25d0)**2)
       end do
       q = 0.d0
       where (r.lt.1.d0)
          q = 0.5d0*(1.d0 + cos(pi*r))
       end where

	case(100) ! uniform field with deformation flow
	  cdf_out = 'weno2d_uniform.nc'
	  q = 1D0

	case(101) ! solid body rotation applied to cylinder
		cdf_out = 'weno2d_rot_cylinder.nc'
		q = 0d0
		do j=1,ny
			r(:,j) = sqrt((x-0.3d0)**2 + (y(j)-0.3d0)**2)
		enddo
		where(r .lt. 0.125d0)
			q = 1d0
		end where
      case(201) ! solid body rotation for cylinder (comparison to frank's code)
        cdf_out = 'weno2d_rot_cylinder_modified.nc'

        xc = xEdge(1)+xWidth/4D0
        yc = yEdge(1)+yWidth/2D0

        q = 0D0
        DO j=1,ny
            r(:,j) = sqrt((x-xc)**2 + (y(j)-yc)**2)
        ENDDO !j
        WHERE(r .lt. 0.25D0)
            q = 1D0
        END WHERE

    end select


    if(transient.or.transient2.or.transient3) tfinal = 5.d0
    if(ntest.eq.17) tfinal = 15.d0

  end subroutine init2d

  ! Initialize DG arrays DGu, DGv, and DGq
  SUBROUTINE DGinit2d(ntest,nex,ney,DG_x,DG_y,DG_xec,DG_yec,xf,yf,nx,ny,DGu,DGu_edge,DGv,DGv_edge,xEdge,yEdge)
	! Computes the initial condtions for the u and v velocity fields for use in modal DG methods
	! u and v must be evaluated at quadrature points, as well as at element interfaces, for each level to be swept through 
	! Outputs are DGu, DGu_edge, DGv, and DGv_edge.
	! DGu contains horizontal velocities evaluated at quadrature nodes in the x direction for each y level
	! DGu_edge contains horizontal velocity evalutated at right edge/interface of each element.
	! DGv and DGv_edge are symmetrically defined, but for y sweeps.
	! -------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: ntest,nex,ney,nx,ny
	REAL(KIND=8), DIMENSION(1:nx,1:ny),INTENT(OUT) :: DGu,DGv
	REAL(KIND=8), DIMENSION(1:nex,1:ny),INTENT(OUT) :: DGu_edge
	REAL(KIND=8), DIMENSION(1:nx,1:ney),INTENT(OUT) :: DGv_edge
	REAL(KIND=8), DIMENSION(0:nx), INTENT(IN) :: xf
	REAL(KIND=8), DIMENSION(0:ny), INTENT(IN) :: yf
	REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: DG_xec
	REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: DG_yec
	REAL(KIND=8), DIMENSION(1:nx), INTENT(IN) :: DG_x
	REAL(KIND=8), DIMENSION(1:ny), INTENT(IN) :: DG_y
    REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: xEdge,yEdge
	REAL(KIND=8), DIMENSION(1:nx,0:ny) :: psi1
	REAL(KIND=8), DIMENSION(0:nx,1:ny) :: psi2
	REAL(KIND=8), DIMENSION(1:nex,0:ny) :: psi1Edge
	REAL(KIND=8), DIMENSION(0:nx,1:ney) :: psi2Edge
	REAL(KIND=8) :: PI,dxe,dye,tmp,xc,yc,xWidth,yWidth
	INTEGER :: i,j,k,cur

    xWidth = xEdge(2)-xEdge(1)
    yWidth = yEdge(2)-yEdge(1)

	PI = DACOS(-1D0)
	dxe = DG_xec(2) - DG_xec(1)
	dye = DG_yec(2) - DG_yec(1)

    SELECT CASE(ntest)
		CASE(1:2) 
		! Uniform flow: u=v=1, no t dependence

		! Evaluate stream function for horizontal velocities
		DO j=0,ny
			DO i=1,nx
				psi1(i,j) = -DG_x(i) + yf(j)
			ENDDO
			DO i=1,nex
				psi1Edge(i,j) = -(DG_xec(i)+dxe/2D0)+yf(j)
			ENDDO
		ENDDO

		! Evaluate stream function for vertical velocities
		DO i=0,nx
			DO j=1,ny
				psi2(i,j) = -xf(i) + DG_y(j)
			ENDDO
			DO j=1,ney
				psi2Edge(i,j) = -xf(i) + (DG_yec(j)+dye/2D0)
			ENDDO
		ENDDO
		CASE(5:6,9,15,100)
		! LeVeque (1996) Deformation

		! Evaluate stream function for horizontal velocities
		DO j=0,ny
			DO i=1,nx
				psi1(i,j) = (1D0/PI)*(DSIN(PI*DG_x(i))**2)*(DSIN(PI*yf(j))**2)
			ENDDO
			DO i=1,nex
				psi1Edge(i,j) = (1D0/PI)*(DSIN(PI*(DG_xec(i)+dxe/2D0))**2)*(DSIN(PI*yf(j))**2)
			ENDDO
		ENDDO

		! Evaluate stream function for vertical velocities
		DO i=0,nx
			DO j=1,ny
				psi2(i,j) = (1D0/PI)*(DSIN(PI*xf(i))**2)*(DSIN(PI*DG_y(j))**2)
			ENDDO
			DO j=1,ney
				psi2Edge(i,j) = (1D0/PI)*(DSIN(PI*xf(i))**2)*(DSIN(PI*(DG_yec(j)+dye/2D0))**2)
			ENDDO
		ENDDO

		CASE(10,101,201)
		! Solid body rotation

        xc = xEdge(1)+xWidth/2D0
        yc = yEdge(1)+yWidth/2D0

		! Evaluate stream function for horizontal velocities
		DO j=0,ny
			DO i=1,nx
				tmp = sqrt((DG_x(i)-xc)**2 + (yf(j)-yc)**2)
				psi1(i,j) = PI*tmp**2
			ENDDO
			DO i=1,nex
				tmp = sqrt((DG_xec(i)+dxe/2d0-xc)**2 + (yf(j)-yc)**2)
				psi1Edge(i,j) = PI*tmp**2
			ENDDO
		ENDDO

		! Evaluate stream function for vertical velocities
		DO i=0,nx
			DO j=1,ny
				tmp = sqrt((xf(i)-xc)**2 + (DG_y(j)-xc)**2)
				psi2(i,j) = PI*tmp**2
			ENDDO
			DO j=1,ney
				tmp = sqrt((xf(i)-xc)**2 + (DG_yec(j)+dye/2d0-yc)**2)
				psi2Edge(i,j) = PI*tmp**2
			ENDDO
		ENDDO

!		CASE(14:17)
!	       do j = 0,ny
!     	      do i = 0,nx
!     	        tmpr = sqrt((x(i)-0.5d0)**2 + (y(j)-0.5d0)**2)
!     	        psi(i,j) = 0.4d0*pi*tmpr**2
!     	        psi2(i,j) = 0.8d0*pi &
!     	             *(0.5d0*tmpr**2 &
!     	               + (1.d0/96.d0)*log(1.d0 - 16.d0*tmpr**2 + 256.d0*tmpr**4) &
!     	               - (1.d0/48.d0)*log(1.d0 + 16.d0*tmpr**2) &
!     	               - sqrt(3.d0)/48.d0*atan((-1.d0 + 32.d0*tmpr**2)/sqrt(3.d0)))
!     	     end do
!     	  end do
		CASE(99)
		! Uniform flow: u=1,v=1, no t dependence

		! Evaluate stream function for horizontal velocities
		DO j=0,ny
			DO i=1,nx
				psi1(i,j) =  -DG_x(i) + yf(j)
			ENDDO
			DO i=1,nex
				psi1Edge(i,j) = -(DG_xec(i)+dxe/2D0)+yf(j)
			ENDDO
		ENDDO

		! Evaluate stream function for vertical velocities
		DO i=0,nx
			DO j=1,ny
				psi2(i,j) = -xf(i) + DG_y(j)
			ENDDO
			DO j=1,ney
				psi2Edge(i,j) = -xf(i) + (DG_yec(j)+dye/2D0)
			ENDDO
		ENDDO

	END SELECT

	! Compute u velocities from stream function
	DO j=1,ny
		DGu(:,j) = (psi1(:,j)-psi1(:,j-1))/(yf(j)-yf(j-1))
		DGu_edge(:,j) = (psi1Edge(:,j)-psi1Edge(:,j-1))/(yf(j)-yf(j-1))
	ENDDO

	! Compute v velocities from stream function
	DO i=1,nx
		DGv(i,:) = -1D0*(psi2(i,:)-psi2(i-1,:))/(xf(i)-xf(i-1))
		DGv_edge(i,:) = -1D0*(psi2Edge(i,:)-psi2Edge(i-1,:))/(xf(i)-xf(i-1))
	ENDDO
	
  END SUBROUTINE DGinit2d
  
  subroutine skamstep_2d(q,dqdt,u0,v0,u2,v2,rho_in,rhoq,rhoprime,nx,ny,npad,dt,jcbn,&
            xlambda,xmonlimit,ylambda,ymonlimit,&
			DGu0,DGuedge0,DGv0,DGvedge0,nex,ney,dxel,dyel,norder,DG_wghts,DG_C,DG_LUC,IPIV,DG_L,DG_DL)

    !  use forward-in-time Skamarock (2006) scheme

    implicit none

    integer, intent(in) :: nx, ny, npad
    real (kind=8), intent(in)  :: dt

    real (kind=8), dimension(1-npad:nx+npad,1-npad:ny+npad), intent(inout) :: q
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: dqdt ! extra tendency
    real (kind=8), dimension(0:nx,1:ny), intent(in) :: u0, u2
    real (kind=8), dimension(1:nx,0:ny), intent(in) :: v0, v2
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: jcbn
    real (kind=8), dimension(0:nx,1:ny), intent(out) :: xlambda,xmonlimit
    real (kind=8), dimension(1:nx,0:ny), intent(out) :: ylambda,ymonlimit

    ! local variables, two-dimensional arrays
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: rho_in
    real (kind=8), dimension(1:nx,1:ny), intent(inout) :: rhoq,rhoprime
    real (kind=8), dimension(0:nx,1:ny) :: u, uh
    real (kind=8), dimension(1:nx,0:ny) :: v, vh

!bloss    integer, parameter :: nmaxcfl = 20
    
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl, &
                             1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rho
    real (kind=8), dimension(0-npad-nmaxcfl:nx+npad+nmaxcfl, &
                             1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rhouh
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl, &
                             0-npad-nmaxcfl:ny+npad+nmaxcfl) :: rhovh

    ! local variables, one-dimensional arrays in x-direction
    real (kind=8), dimension(0:nx+1) :: rhoq1dx
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: q1dx
    real (kind=8), dimension(-2:nx+2) :: rhouh1d
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: rho1dx
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: rhoprime1dx
    

    ! local variables, one-dimensional arrays in y-direction
    real (kind=8), dimension(0:ny+1) :: rhoq1dy
    real (kind=8), dimension(1-npad-nmaxcfl:ny+npad+nmaxcfl) :: q1dy
    real (kind=8), dimension(-2:ny+2) :: rhovh1d
    real (kind=8), dimension(1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rho1dy
    real (kind=8), dimension(1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rhoprime1dy

    real(kind=8) :: scale
    integer, dimension(2) :: xbctype, ybctype
    real (kind=8), dimension(2) :: fxbc, fybc ! specified bndry flux 
    real (kind=8), dimension(0:nx  ,ny) :: xflx
    real (kind=8), dimension(nx  ,0:ny) :: yflx

    integer i,j, ii, npadq, npadu, npadrho, npadrp, npad2

    real(kind=8), external :: tfcn
    real(kind=8) :: t_temp, tmpflx(0:ny), tmplam(0:ny), tmpmon(0:ny)

	! DG parameters
	INTEGER, INTENT(IN) :: nex,ney,norder
    REAL(KIND=8), DIMENSION(1:nx,1:ny), INTENT(IN) :: DGu0, DGv0 ! u and v velocities at DG nodes
	REAL(KIND=8), DIMENSION(1:nex,1:ny), INTENT(IN) :: DGuedge0
	REAL(KIND=8), DIMENSION(1:nx,1:ney), INTENT(IN) :: DGvedge0
    REAL(KIND=8), DIMENSION(1:3,1:nx,1:ny) :: DGu,DGv ! u and v velocities at DG nodes at ssprk3 substep times
	REAL(KIND=8), DIMENSION(1:3,1:nex,1:ny) :: DGuedge
	REAL(KIND=8), DIMENSION(1:3,1:nx,1:ney) :: DGvedge

    REAL(KIND=8), DIMENSION(0:norder,0:norder), INTENT(IN):: DG_C,DG_LUC,DG_L,DG_DL
	INTEGER, DIMENSION(0:norder), INTENT(IN) :: IPIV
    REAL(KIND=8), DIMENSION(0:norder), INTENT(IN) :: DG_wghts
	REAL(KIND=8), INTENT(IN) :: dxel,dyel
	REAL(KIND=8), DIMENSION(1:nx) :: DG_rhoq1dx,DG_rhop1dx
    REAL(KIND=8), DIMENSION(1:3,1:nx) :: DGu1dx
	REAL(KIND=8), DIMENSION(1:3,1:nex) :: DGuedge1dx
	REAL(KIND=8), DIMENSION(1:ny) :: DG_rhoq1dy,DG_rhop1dy
    REAL(KIND=8), DIMENSION(1:3,1:ny) :: DGv1dy
	REAL(KIND=8), DIMENSION(1:3,1:ney) :: DGvedge1dy

	REAL(KIND=8) :: tstar
	LOGICAL :: dorhoupdate
	REAL(KIND=8), DIMENSION(1:nx) :: jcbn1dx
	REAL(KIND=8), DIMENSION(1:ny) :: jcbn1dy

    ! set up boundary types, fluxes, limiting flags
    if(nofluxew) then
       xbctype = 1 ! try open bc
    else
       xbctype = 2           ! periodic boundary conditions
    end if
    if(nofluxns) then
       ybctype = 1 ! try open bc
    else
       ybctype = 2           ! periodic boundary conditions
    end if

    fxbc = 0.d0          ! initialize in case fixed flux is used.
    fybc = 0.d0          ! initialize in case fixed flux is used.

!!$    if(nofluxew) xbctype = 3
!!$    if(nofluxns) xbctype = 3

    scale = 1.

    t_temp = time+0.5d0*dt
    uh = u0
    vh = v0
    IF(nmethod .eq. 99) THEN
		DGu(1,:,:) = DGu0(:,:)
		DGu(2,:,:) = DGu0(:,:)
		DGu(3,:,:) = DGu0(:,:)

		DGuedge(1,:,:) = DGuedge0(:,:)
		DGuedge(2,:,:) = DGuedge0(:,:)
		DGuedge(3,:,:) = DGuedge0(:,:)
    END IF
	IF(nmethod2 .eq. 99) THEN
		DGv(1,:,:) = DGv0(:,:)
		DGv(2,:,:) = DGv0(:,:)
		DGv(3,:,:) = DGv0(:,:)

		DGvedge(1,:,:) = DGvedge0(:,:)
		DGvedge(2,:,:) = DGvedge0(:,:)
		DGvedge(3,:,:) = DGvedge0(:,:)
	END IF

    if(transient) then
       uh = uh*tfcn(t_temp)
       vh = vh*tfcn(t_temp)
       IF(nmethod .eq. 99) THEN 
		tstar = time
		DGu(1,:,:) = DGu(1,:,:)*tfcn(tstar)
		DGuedge(1,:,:) = DGuedge(1,:,:)*tfcn(tstar)

		tstar = time+dt
		DGu(2,:,:) = DGu(2,:,:)*tfcn(tstar)
		DGuedge(2,:,:) = DGuedge(2,:,:)*tfcn(tstar)

		tstar = time+0.5D0*dt
		DGu(3,:,:) = DGu(3,:,:)*tfcn(tstar)
		DGuedge(3,:,:) = DGuedge(3,:,:)*tfcn(tstar)
       END IF
	   IF(nmethod2 .eq. 99) THEN
		tstar = time
		DGv(1,:,:) = DGv(1,:,:)*tfcn(tstar)
		DGvedge(1,:,:) = DGvedge(1,:,:)*tfcn(tstar)

		tstar = time+dt
		DGv(2,:,:) = DGv(2,:,:)*tfcn(tstar)
		DGvedge(2,:,:) = DGvedge(2,:,:)*tfcn(tstar)

		tstar = time+0.5D0*dt
		DGv(3,:,:) = DGv(3,:,:)*tfcn(tstar)
		DGvedge(3,:,:) = DGvedge(3,:,:)*tfcn(tstar)

	   END IF
    end if

    if(transient2) then
       uh = u0*tfcn(t_temp)**2 + u2*(1.d0-tfcn(t_temp)**2)
       vh = v0*tfcn(t_temp)**2 + v2*(1.d0-tfcn(t_temp)**2)
    end if

    if(transient3) then
       t_temp = 2.d0*t_temp
       uh = u0 + u2*tfcn(t_temp)
       vh = v0 + v2*tfcn(t_temp)
    end if
    ! set padding size for each variable
    npadrho = npad-1 + (nmaxcfl+1)/2
    npadu = 2 + (nmaxcfl+1)/2

    npadq = npad+nmaxcfl-1
    npadrp = nmaxcfl+2

    ! set up density array for back trajectory computation
    rho(:,:) = 0.d0 ! zero out
    rho(1:nx,1:ny) = rho_in(:,:)

    ! x-mass flux at time t + dt/2
    rhouh(:,:) = 0.d0 ! zero out
    rhouh(0:nx,1:ny) = uh(:,:)

    ! y-mass flux at time t + dt/2
    rhovh(:,:) = 0.d0 ! zero out
    rhovh(1:nx,0:ny) = vh(:,:)
    if(nofluxew) then ! no scalar flux at x boundaries
       do i = 1,npadrho
          rho(1-i,1:ny) = rho(1,1:ny)
          rho(nx+i,1:ny) = rho(nx,1:ny)
       end do
       do i = 1,npadu
          rhouh(-i,1:ny) = rhouh(0,1:ny)
          rhouh(nx+i,1:ny) = rhouh(nx,1:ny)
       end do
    else  ! periodic in x
       rho(1-npadrho:0,1:ny) = rho(nx+1-npadrho:nx,1:ny)
       rho(nx+1:nx+npadrho,1:ny) = rho(1:npadrho,1:ny)

       rhouh(-npadu:-1,1:ny) = rhouh(nx+1-npadu:nx,1:ny)
       rhouh(nx+1:nx+npadu,1:ny) = rhouh(1:npadu,1:ny)
    end if
    if(nofluxns) then ! no scalar flux at y boundaries
       do i = 1,npadrho
          rho(:,1-i) = rho(:,1)
          rho(:,ny+i) = rho(:,ny)
       end do
       do i = 1,npadu
          rhovh(:,-i) = rhovh(:,0)
          rhovh(:,ny+i) = rhovh(:,ny)
       end do
    else ! periodic in y
       rho(:,1-npadrho:0) = rho(:,ny+1-npadrho:ny)
       rho(:,ny+1:ny+npadrho) = rho(:,1:npadrho)
       rhovh(:,-npadu:-1) = rhovh(:,ny+1-npadu:ny)
       rhovh(:,ny+1:ny+npadu) = rhovh(:,1:npadu)
    end if
    q1dx = 0.d0
    rhoprime1dx = 0.d0
    rho1dx = 0.d0
    rhouh1d = 0.d0

    q1dy = 0.d0
    rhoprime1dy = 0.d0
    rho1dy = 0.d0
    rhovh1d = 0.d0
    if(oddstep) then

       ! perform sweeps in x-direction
       do j = 1,ny
          q1dx(1:nx)     = q(1:nx,j)
          rhoq1dx(1:nx)  = rhoq(1:nx,j)
          rhoprime1dx(1:nx) = rhoprime(1:nx,j)

          IF(nmethod .eq. 99) THEN
		   dorhoupdate = .TRUE.
		   jcbn1dx(1:nx) = jcbn(1:nx,j)

		   DGu1dx(1:3,1:nx) = DGu(1:3,1:nx,j)	
		   DGuedge1dx(1:3,1:nex) = DGuedge(1:3,1:nex,j)	  
 
		   ! Take values from rhoq and rhoprime arrays (cell averages)
		   DG_rhoq1dx(1:nx) = rhoq1dx(1:nx)
		   DG_rhop1dx(1:nx) = rhoprime1dx(1:nx)
		  END IF

          if(nofluxew) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dx(1-ii) = MAX(0.d0,q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2)))
                   q1dx(nx+ii) = MAX(0.d0,q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1)))
                end do
             else
                do ii = 1,npadq
                   q1dx(1-ii) = q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2))
                   q1dx(nx+ii) = q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dx(1-npadrp:0) = rhoprime1dx(1)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(nx)
             
             ! fill in rhoq values 
             rhoq1dx(0) = rhoprime1dx(0)*q1dx(0)
             rhoq1dx(nx+1) = rhoprime1dx(nx+1)*q1dx(nx+1)

          else ! periodic in x
             rhoq1dx(0) = rhoq1dx(nx)
             rhoq1dx(nx+1) = rhoq1dx(1)

             q1dx(1-npadq:0) = q1dx(nx+1-npadq:nx)
             q1dx(nx+1:nx+npadq) = q1dx(1:npadq)

             rhoprime1dx(1-npadrp:0) = rhoprime1dx(nx+1-npadrp:nx)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dx(1-npadrho:nx+npadrho) = rho(1-npadrho:nx+npadrho,j)
          rhouh1d(-2:nx+2) = rhouh(-2:nx+2,j)

          call ppmwrap(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx,xflx(0,j),dt, &
                       nx,npad,nmaxcfl,xbctype,fxbc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod,lambdamax,epslambda, &
                       xlambda(0,j),xmonlimit(0,j), &
                       DG_rhoq1dx,DG_rhop1dx,DGu1dx,DGuedge1dx,nex,dxel,norder,DG_wghts,DG_C,DG_LUC,IPIV,& 
					   DG_L,DG_DL,dorhoupdate,jcbn1dx)

          ! update solution
		  IF(nmethod .eq. 99) THEN
		   rhoq(1:nx,j) = DG_rhoq1dx(1:nx)
		   rhoprime(1:nx,j) = DG_rhop1dx(1:nx)
		   q(1:nx,j) = DG_rhoq1dx(1:nx)/DG_rhop1dx(1:nx)
		  ELSE
           q(1:nx,j) = rhoq1dx(1:nx)/rhoprime1dx(1:nx)
           rhoprime(1:nx,j) = rhoprime1dx(1:nx)
           rhoq(1:nx,j) = rhoq1dx(1:nx)
		  END IF

       end do

       ! perform sweeps in y-direction
       do i = 1,nx
          q1dy(1:ny)     = q(i,1:ny)
          rhoq1dy(1:ny)  = rhoq(i,1:ny)
          rhoprime1dy(1:ny) = rhoprime(i,1:ny)

          IF(nmethod2 .eq. 99) THEN
		   dorhoupdate = .FALSE.
		   jcbn1dy(1:ny) = jcbn(i,1:ny)

		   DGv1dy(1:3,1:ny) = DGv(1:3,i,1:ny)
		   DGvedge1dy(1:3,1:ney) = DGvedge(1:3,i,1:ney)

		   DG_rhoq1dy(1:ny) = rhoq1dy(1:ny)
	 	   DG_rhop1dy(1:ny) = rhoprime1dy(1:ny)
		  END IF

          if(nofluxns) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dy(1-ii) = MAX(0.d0,q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2)))
                   q1dy(ny+ii) = MAX(0.d0,q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1)))
                end do
             else
                do ii = 1,npadq
                   q1dy(1-ii) = q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2))
                   q1dy(ny+ii) = q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dy(1-npadrp:0) = rhoprime1dy(1)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(ny)
             
             ! fill in rhoq values 
             rhoq1dy(0) = rhoprime1dy(0)*q1dy(0)
             rhoq1dy(ny+1) = rhoprime1dy(ny+1)*q1dy(ny+1)

          else ! periodic in y
             rhoq1dy(0) = rhoq1dy(ny)
             rhoq1dy(ny+1) = rhoq1dy(1)

             q1dy(1-npadq:0) = q1dy(ny+1-npadq:ny)
             q1dy(ny+1:ny+npadq) = q1dy(1:npadq)

             rhoprime1dy(1-npadrp:0) = rhoprime1dy(ny+1-npadrp:ny)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dy(1-npadrho:ny+npadrho) = rho(i,1-npadrho:ny+npadrho)
          rhovh1d(-2:ny+2) = rhovh(i,-2:ny+2)

          call ppmwrap(rhoq1dy,q1dy,rhovh1d,rho1dy,rhoprime1dy,tmpflx,dt, &
                       ny,npad,nmaxcfl,ybctype,fybc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod2,lambdamax,epslambda,tmplam, tmpmon, &
                       DG_rhoq1dy,DG_rhop1dy,DGv1dy,DGvedge1dy,ney,dyel,norder,DG_wghts,DG_C,DG_LUC,IPIV, &
					   DG_L,DG_DL,dorhoupdate,jcbn1dy)

          yflx(i,:) = tmpflx(:)
          ylambda(i,:) = tmplam(:)
          ymonlimit(i,:) = tmpmon(:)

          ! update solution
		  IF(nmethod2 .eq. 99) THEN
		   rhoq(i,1:ny) = DG_rhoq1dy(1:ny)
		   rhoprime(i,1:ny) = DG_rhop1dy(1:ny)
		   q(i,1:ny) = DG_rhoq1dy(1:ny)/DG_rhop1dy(1:ny)
		  ELSE
           q(i,1:ny) = rhoq1dy(1:ny)/rhoprime1dy(1:ny)
           rhoprime(i,1:ny) = rhoprime1dy(1:ny)
           rhoq(i,1:ny) = rhoq1dy(1:ny)
		  END IF


       end do
    else

       ! perform sweeps in y-direction
       do i = 1,nx
          q1dy(1:ny)     = q(i,1:ny)
          rhoq1dy(1:ny)  = rhoq(i,1:ny)
          rhoprime1dy(1:ny) = rhoprime(i,1:ny)

          IF(nmethod2 .eq. 99) THEN
		   dorhoupdate = .TRUE.
		   jcbn1dy(1:ny) = jcbn(i,1:nx)

		   DGv1dy(1:3,1:ny) = DGv(1:3,i,1:ny)
		   DGvedge1dy(1:3,1:ney) = DGvedge(1:3,i,1:ney)

		   DG_rhoq1dy(1:ny) = rhoq1dy(1:ny)
	 	   DG_rhop1dy(1:ny) = rhoprime1dy(1:ny)
		  END IF

          if(nofluxns) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dy(1-ii) = MAX(0.d0,q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2)))
                   q1dy(ny+ii) = MAX(0.d0,q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1)))
                end do
             else
                do ii = 1,npadq
                   q1dy(1-ii) = q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2))
                   q1dy(ny+ii) = q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dy(1-npadrp:0) = rhoprime1dy(1)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(ny)
             
             ! fill in rhoq values 
             rhoq1dy(0) = rhoprime1dy(0)*q1dy(0)
             rhoq1dy(ny+1) = rhoprime1dy(ny+1)*q1dy(ny+1)

          else ! periodic in y
             rhoq1dy(0) = rhoq1dy(ny)
             rhoq1dy(ny+1) = rhoq1dy(1)

             q1dy(1-npadq:0) = q1dy(ny+1-npadq:ny)
             q1dy(ny+1:ny+npadq) = q1dy(1:npadq)

             rhoprime1dy(1-npadrp:0) = rhoprime1dy(ny+1-npadrp:ny)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dy(1-npadrho:ny+npadrho) = rho(i,1-npadrho:ny+npadrho)
          rhovh1d(-2:ny+2) = rhovh(i,-2:ny+2)
          call ppmwrap(rhoq1dy,q1dy,rhovh1d,rho1dy,rhoprime1dy,tmpflx,dt, &
                       ny,npad,nmaxcfl,ybctype,fybc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod2,lambdamax,epslambda,tmplam, tmpmon, &
                       DG_rhoq1dy,DG_rhop1dy,DGv1dy,DGvedge1dy,ney,dyel,norder,DG_wghts,DG_C,DG_LUC,IPIV,&
					   DG_L,DG_DL,dorhoupdate,jcbn1dy)


          yflx(i,:) = tmpflx(:)
          ylambda(i,:) = tmplam(:)
          ymonlimit(i,:) = tmpmon(:)

          ! update solution
		  IF(nmethod2 .eq. 99) THEN
		   rhoq(i,1:ny) = DG_rhoq1dy(1:ny)
		   rhoprime(i,1:ny) = DG_rhop1dy(1:ny)
		   q(i,1:ny) = DG_rhoq1dy(1:ny)/DG_rhop1dy(1:ny)
		  ELSE
           q(i,1:ny) = rhoq1dy(1:ny)/rhoprime1dy(1:ny)
           rhoprime(i,1:ny) = rhoprime1dy(1:ny)
           rhoq(i,1:ny) = rhoq1dy(1:ny)
		  END IF


       end do

       ! perform sweeps in x-direction
       do j = 1,ny
          q1dx(1:nx)     = q(1:nx,j)
          rhoq1dx(1:nx)  = rhoq(1:nx,j)
          rhoprime1dx(1:nx) = rhoprime(1:nx,j)

          IF(nmethod .eq. 99) THEN
		   dorhoupdate = .FALSE.
		   jcbn1dx(1:nx) = jcbn(1:nx,j)

		   DGu1dx(1:3,1:nx) = DGu(1:3,1:nx,j)	
		   DGuedge1dx(1:3,1:nex) = DGuedge(1:3,1:nex,j)	  
 
		   ! Take values from rhoq and rhoprime arrays (cell averages)
		   DG_rhoq1dx(1:nx) = rhoq1dx(1:nx)
		   DG_rhop1dx(1:nx) = rhoprime1dx(1:nx)
		  END IF

          if(nofluxew) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dx(1-ii) = MAX(0.d0,q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2)))
                   q1dx(nx+ii) = MAX(0.d0,q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1)))
                end do
             else
                do ii = 1,npadq
                   q1dx(1-ii) = q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2))
                   q1dx(nx+ii) = q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dx(1-npadrp:0) = rhoprime1dx(1)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(nx)
             
             ! fill in rhoq values 
             rhoq1dx(0) = rhoprime1dx(0)*q1dx(0)
             rhoq1dx(nx+1) = rhoprime1dx(nx+1)*q1dx(nx+1)

          else ! periodic in x
             rhoq1dx(0) = rhoq1dx(nx)
             rhoq1dx(nx+1) = rhoq1dx(1)

             q1dx(1-npadq:0) = q1dx(nx+1-npadq:nx)
             q1dx(nx+1:nx+npadq) = q1dx(1:npadq)

             rhoprime1dx(1-npadrp:0) = rhoprime1dx(nx+1-npadrp:nx)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dx(1-npadrho:nx+npadrho) = rho(1-npadrho:nx+npadrho,j)
          rhouh1d(-2:nx+2) = rhouh(-2:nx+2,j)

          call ppmwrap(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx,xflx(0,j),dt, &
                       nx,npad,nmaxcfl,xbctype,fxbc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod,lambdamax,epslambda, &
                       xlambda(0,j),xmonlimit(0,j), &
                       DG_rhoq1dx,DG_rhop1dx,DGu1dx,DGuedge1dx,nex,dxel,norder,DG_wghts,DG_C,DG_LUC,IPIV, &
					   DG_L,DG_DL,dorhoupdate,jcbn1dx)


          ! update solution
		  IF(nmethod .eq. 99) THEN
		   rhoq(1:nx,j) = DG_rhoq1dx(1:nx)
		   rhoprime(1:nx,j) = DG_rhop1dx(1:nx)
		   q(1:nx,j) = DG_rhoq1dx(1:nx)/DG_rhop1dx(1:nx)
		  ELSE
           q(1:nx,j) = rhoq1dx(1:nx)/rhoprime1dx(1:nx)
           rhoprime(1:nx,j) = rhoprime1dx(1:nx)
           rhoq(1:nx,j) = rhoq1dx(1:nx)
		  END IF

       end do
    end if

  end subroutine skamstep_2d

! QUICK JUMP

  subroutine output2d(q,u,v,xlam,xmnl,ylam,ymnl,&
       nx,ny,x,xf,y,yf,time_in,status,cdf_out,ilevel,&
       muIn,quadNode,quadWeight,norder)
    implicit none

    !     VMS include statement (on UW Vax machines)
    !     include 'local_root:[include]netcdf.inc'

    integer, intent(in) :: nx, ny, status, ilevel
    real (kind=8), intent(in) :: time_in
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: q
    real (kind=8), dimension(0:nx,1:ny), intent(in) :: u
    real (kind=8), dimension(1:nx,0:ny), intent(in) :: v
    real (kind=8), dimension(0:nx,1:ny), intent(in) :: xlam
    real (kind=8), dimension(1:nx,0:ny), intent(in) :: ylam
    real (kind=8), dimension(0:nx,1:ny), intent(in) :: xmnl
    real (kind=8), dimension(1:nx,0:ny), intent(in) :: ymnl
    real (kind=8), dimension(0:nx), intent(in) :: xf
    real (kind=8), dimension(0:ny), intent(in) :: yf
    real (kind=8), dimension(1:nx), intent(in) :: x
    real (kind=8), dimension(1:ny), intent(in) :: y
    INTEGER, INTENT(IN) :: norder
    REAL(KIND=8), INTENT(IN) :: muIn
    REAL(KIND=8), DIMENSION(0:norder) :: quadNode,quadWeight
    character(len=40), intent(in) :: cdf_out ! Name of the netCDF file

    integer msize,n2d,i,j,ierr,idq,idu,idv,idt,idxl,idxm,idyl,idym,idnode,idweight,idmu

    real (kind=4) :: time4
    real (kind=8) :: tfcn
    real (kind=4), dimension(nx)   :: temp
    real (kind=4), dimension(0:nx) :: temp1
    real (kind=4), dimension(:), allocatable :: temp2
    real (kind=4), dimension(1:ny) :: temp3
    real (kind=4), dimension(0:ny) :: temp4

    !  netCDF file declaration
    integer cdfid               ! ID for the netCDF file to be created
    character(len=8):: xname,xfname,nxname,nxfname,uname, &
         yname,yfname,nyname,nyfname,vname, &
         tname,ntname,qname,xmname,xlname,ymname,ylname,muname

    !  netCDF variables declaration

    integer :: &
    &    start(4) &             ! netCDF stuff;
    &   ,count(4) &             ! netCDF stuff;
    &   ,vdim(4)  &             ! netCDF stuff;
    &   ,ndim     &             ! netCDF stuff;
    &   ,idx,idy,idxf,idyf,idnx,idny,idnxf,idnyf,idnt,node_dimid
    data       &             
    &    start /1, 1, 1, 1/, count /1, 1, 1, 1/

    integer nxp1, nyp1
    save cdfid, idq, idu, idv, idt, idnt, idxl, idxm, idyl, idym

    !     UNIX include statement (on UW Unix machines)
    include 'netcdf.inc'

   !  Begin

    if (status.eq.-1) then
       !                  Create netCDF file
       ierr = NF_CREATE(cdf_out, NF_CLOBBER, cdfid)


       ! set up time dimension
       tname = 'time'
       ntname = 'nt'
       ierr = NF_REDEF(cdfid)
       ierr = NF_DEF_DIM(cdfid,TRIM(ntname),ilevel+1,idnt)
       ierr = NF90_DEF_DIM(cdfid, "nnodes", norder+1, node_dimid)

	   ierr = NF90_DEF_VAR(cdfid, "qweights",NF90_FLOAT, node_dimid, idweight)
	   ierr = NF90_DEF_VAR(cdfid, "qnodes",NF90_FLOAT, node_dimid, idnode)
       ierr = NF_DEF_VAR(cdfid,TRIM(tname), NF_FLOAT, 1, idnt, idt)
       ierr = NF_ENDDEF(cdfid)

       ALLOCATE(temp2(1:ilevel+1),STAT=ierr)
       temp2(1) = 0.d0
       do i = 1,ilevel
          temp2(i+1) = DBLE(i)*time_in/DBLE(ilevel)
       end do
       ierr = NF_PUT_VAR_REAL(cdfid, idt, temp2)
       ierr = NF90_PUT_VAR(cdfid,idweight,quadWeight)
       ierr = NF90_PUT_VAR(cdfid,idnode,quadNode)

       DEALLOCATE(temp2)
       
       return

    elseif (status.eq.0) then

       write(xname,'(a1,i1)') 'x', ilevel
       write(xfname,'(a2,i1)') 'xf', ilevel
       write(nxname,'(a2,i1)') 'nx', ilevel
       write(nxfname,'(a3,i1)') 'nxf', ilevel
       write(yname,'(a1,i1)') 'y', ilevel
       write(yfname,'(a2,i1)') 'yf', ilevel
       write(nyname,'(a2,i1)') 'ny', ilevel
       write(nyfname,'(a3,i1)') 'nyf', ilevel
       write(qname,'(a1,i1)') 'Q', ilevel
       write(uname,'(a1,i1)') 'U', ilevel
       write(vname,'(a1,i1)') 'V', ilevel
       write(xlname,'(a4,i1)') 'XLAM', ilevel
       write(xmname,'(a4,i1)') 'XMON', ilevel
       write(ylname,'(a4,i1)') 'YLAM', ilevel
       write(ymname,'(a4,i1)') 'YMON', ilevel
       WRITE(muname, '(a2,i1)') 'mu',ilevel

       nxp1=nx+1
       nyp1=ny+1

       ierr = NF_REDEF(cdfid)

       ierr = NF_DEF_DIM(cdfid,TRIM(nxname),nx,idnx)
       ierr = NF_DEF_DIM(cdfid,TRIM(nxfname),nxp1,idnxf)
       ierr = NF_DEF_DIM(cdfid,TRIM(nyname),ny,idny)
       ierr = NF_DEF_DIM(cdfid,TRIM(nyfname),nyp1,idnyf)

       ierr = NF_DEF_VAR(cdfid,TRIM(xname), NF_FLOAT, 1, idnx, idx)
       ierr = NF_DEF_VAR(cdfid,TRIM(xfname), NF_FLOAT, 1, idnxf, idxf)
       ierr = NF_DEF_VAR(cdfid,TRIM(yname), NF_FLOAT, 1, idny, idy)
       ierr = NF_DEF_VAR(cdfid,TRIM(yfname), NF_FLOAT, 1, idnyf, idyf)
       ierr = NF90_DEF_VAR(cdfid, TRIM(muname),NF90_FLOAT,idmu)

       ndim = 3
       vdim(1) = idnx
       vdim(2) = idny
       vdim(3) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(qname),NF_FLOAT,ndim,vdim,idq)

       ndim = 3
       vdim(1) = idnxf
       vdim(2) = idny
       vdim(3) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(uname),NF_FLOAT,ndim,vdim,idu)
       ierr = NF_DEF_VAR(cdfid,TRIM(xlname),NF_FLOAT,ndim,vdim,idxl)
       ierr = NF_DEF_VAR(cdfid,TRIM(xmname),NF_FLOAT,ndim,vdim,idxm)

       ndim = 3
       vdim(1) = idnx
       vdim(2) = idnyf
       vdim(3) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(vname),NF_FLOAT,ndim,vdim,idv)
       ierr = NF_DEF_VAR(cdfid,TRIM(ylname),NF_FLOAT,ndim,vdim,idyl)
       ierr = NF_DEF_VAR(cdfid,TRIM(ymname),NF_FLOAT,ndim,vdim,idym)

       ierr = NF_ENDDEF(cdfid)

       temp = x
       ierr = NF_PUT_VAR_REAL(cdfid, idx, temp)

       temp1 = xf
       ierr = NF_PUT_VAR_REAL(cdfid, idxf, temp1)

       temp3 = y
       ierr = NF_PUT_VAR_REAL(cdfid, idy, temp3)

       temp4 = yf
       ierr = NF_PUT_VAR_REAL(cdfid, idyf, temp4)

       ierr = NF90_PUT_VAR(cdfid,idmu,muIn)

       start(3) = 1

       !                 Open netCDF file:

!!$       idq=ncvid(cdfid,TRIM(qname), ierr)
!!$       idu=ncvid(cdfid,TRIM(uname), ierr)
!!$       idt=ncvid(cdfid,TRIM(tname) ,ierr)

    else if(status.eq.2) then
       !         Close netCDF
       call ncclos(cdfid, ierr)
       return
    end if

    !               Output fields
    !              concentration field
    count(1)=nx
    do j = 1,ny
       start(2) = j
       temp = q(:,j)
       call ncvpt(cdfid, idq , start, count, temp, ierr)
    end do

    !              u field
    count(1)=nx+1
    do j = 1,ny
       start(2) = j
       temp1 = u(:,j)
       call ncvpt(cdfid, idu , start, count, temp1, ierr)
    end do

    !              v field
    count(1)=nx
    do j = 0,ny
       start(2) = j+1
       temp = v(:,j)
       call ncvpt(cdfid, idv , start, count, temp, ierr)
    end do

    !              lambda field for x-advection
    count(1)=nx+1
    do j = 1,ny
       start(2) = j
       temp1 = xlam(:,j)
       call ncvpt(cdfid, idxl , start, count, temp1, ierr)
    end do

    !              lambda field for y-advection
    count(1)=nx
    do j = 0,ny
       start(2) = j+1
       temp = ylam(:,j)
       call ncvpt(cdfid, idyl , start, count, temp, ierr)
    end do

    !              monotonic limiting indicator  for x-advection
    count(1)=nx+1
    do j = 1,ny
       start(2) = j
       temp1 = xmnl(:,j)
       call ncvpt(cdfid, idxm , start, count, temp1, ierr)
    end do

    !              monotonic limiting indicator  for y-advection
    count(1)=nx
    do j = 0,ny
       start(2) = j+1
       temp = ymnl(:,j)
       call ncvpt(cdfid, idym , start, count, temp, ierr)
    end do

    start(3) = start(3) + 1

  end subroutine output2d


end program execute
