!---------------------------------------------------------------             
! solving the radial schroedinger equation by Numerov method
! Atomic (Ry) units
! by Y.D.Chen, 2015.11.03
!---------------------------------------------------------------
module materials
!---------------------------------------------------------------

    integer, parameter :: dp = selected_real_kind(18,250)
    real (dp) ,parameter :: a = 3.16_dp !in Angstrom
    real (dp) ,parameter :: Eg = 2.41_dp
    
    real (dp) ,parameter :: aB = 0.52917721092_dp !in Angstrom
    real (dp) ,parameter :: E1 = -13.60569253_dp !in eV
    
    real (dp), parameter :: PI = 3.141592653_dp
    real (dp), parameter :: ElectronCharge = 1.60217657000000e-19
    real (dp), parameter :: hbar = 1.05457172600000e-34
    
    !Coulumb potential
    real (dp), parameter :: x2D = 6.03_dp
    real (dp), parameter :: gamma = -0.1159315157_dp !gamma-ln2
    
    real (dp) ,parameter :: epsilon = 4.13_dp
    real (dp) ,parameter :: epsilon0 = 8.854187817620e-12
    
    real (dp), parameter :: Ea1s = 2.08695652173913
    
    real (dp), parameter :: me = 9.109383560000000e-31
    real (dp), parameter :: mu = 0.16 !in me
    real (dp), parameter :: mcm = 4 * mu
    real (dp), parameter :: hbar2mu = hbar**2/(mu*me)/ElectronCharge/(1e-20)
    real (dp), parameter :: hbar2M = hbar**2/(mcm*me)/ElectronCharge/(1e-20)
end module
!---------------------------------------------------------------
!---------------------------------------------------------------
program magneticexciton
!---------------------------------------------------------------
    use materials
    !USE mkl95_LAPACK    
    !USE mkl95_PRECISION, ONLY: WP=>SP
    USE lapack95
    USE f95_precision, ONLY: WP=>SP

    implicit none
  
    integer:: mesh
    integer n, m, m_abs, i, s, nindex
    !label n,m in 1 dimension
    integer, parameter :: Nnindex = 12
    integer, parameter :: Nnm = Nnindex ** 2

    ! grid initialization and rescale 
    real (dp) :: rmin, rmax, xmin, xmax, dx
    real (dp), allocatable :: r(:), sqr(:), r2(:), y(:)

    real (dp) :: e, e_const
    ! store e0, phi0
    real (dp), allocatable :: e0(:)
    real (dp), allocatable :: eP0(:)
    real (dp), allocatable :: phi0(:,:)
    real (dp) :: norm, perp
    
    ! Hn1m1n2m2 matrix for diagonalization
    integer :: s1, s2, n1, m1, n2, m2
    integer, allocatable :: nlabel(:), mlabel(:)
    real*8, allocatable :: Hnm(:,:)
    real*8, allocatable :: Htemp(:,:)!, HtempB(:,:) !store off-diagonal integrated elements
    real*8, allocatable :: Enm(:,:)
    real (dp) :: Hprime, Hc!, HprimeB, HcB
    integer :: info
    
    ! rearrange and output
    integer :: eigenindex, eigenindex_new, eigenindex_lw
    real*8 :: Eave, Eavelw, Hswap
    integer, parameter :: ntot = 3
    integer, parameter :: nout = ntot ** 2 !number of output energy levels
    
    ! skip point
    integer :: plotnull = 0
    
    ! magnetic field setup
    integer, parameter :: bN = 1 !bN is the total number of Bb
    integer :: b !b is the index of Bb
    real (dp) :: Bb, Bbmax, Bbmin, dBb
    
    ! initialize interval of P and etc.
    integer, parameter :: PN = 2
    integer p
    real (dp) :: Pp, Ppmin, Ppmax, dPp
    
    ! output filename variable
    character(len = 4) :: cTemp, cTemp1, cTempout

    Ppmin= 0.0D0
    Ppmax= 1.0D0/a*pi*0.25
    dPp= (Ppmax-Ppmin)/PN
    ! allocate logarithmic mesh
    !
    xmin = -8.0_dp
    rmax = 5*dble(Nnindex - 0.5)**2 * (Nnindex) + 50 * dble(Nnindex - 0.5)
    mesh = 2e3 * Nnindex
    
    allocate (e0(Nnm))
    allocate (eP0(Nnm))
    allocate (phi0(Nnm,0:mesh))
    allocate (nlabel(Nnm))
    allocate (mlabel(Nnm))
    allocate (r(0:mesh), sqr(0:mesh), r2(0:mesh), y(0:mesh))
    allocate (Hnm(Nnm,Nnm))
    allocate (Htemp(Nnm,Nnm))!, HtempB(Nnm,Nnm))
    allocate (Enm(PN+1,Nnm))
    !initialization

    r=0
    sqr=0
    r2=0
    
    dx = (log(rmax) - xmin)/dble(mesh)
    call do_mesh ( xmin, dx, mesh, r, sqr, r2 )   
    
    
    !scan magnetic field
 
        Bb = 100
	
	    !calculate H0
	    Hnm=0
        Htemp=0
        !HtempB=0
        Enm=0
        y=0
        e0=0
        phi0=0
        open(unit=15,file='error.dat')
        open(unit=14,file='convergence.dat')
        open(unit=16, file='Enm.dat')
        !$OMP PARALLEL DO
        s=1
        do nindex=1,Nnindex,1  
            do m_abs=0,nindex-1,1 ! read number of nodes (stop if nodes < 0)
                write(cTemp,'(i4)') int(m_abs) !label the file name to distinguish different B   !(f3.0)
                n = nindex-m_abs-1
                
                write(*,*)"n=",n,"m_abs=",m_abs, "mesh=",mesh
                write(cTemp1,'(i4)') int(n)
                !open(unit=13, file='Phinm'//trim(adjustl(cTemp1))//trim(adjustl(cTemp))//'.dat')
                
                call solve_sheq (n, m_abs, e, e_const, Bb, mesh, dx, r, sqr, r2, y, plotnull)
                write (16,'(i10,i10,e20.10)')n,m_abs,e*e_const

                if (plotnull == 1) then
                    write(*,*)"calculating error, n=",n," m_abs=",m_abs
                end if
                
                !write in +m
                m = m_abs
                e0(s)=e*e_const
	            do i=0,mesh,1
	                phi0(s,i)=y(i)
	            end do
	            nlabel(s)=n
	            mlabel(s)=m
	            s=s+1
                !write in -m
                if (m_abs /= 0) then
                    m = -m_abs
                    e0(s)=e*e_const
	                do i=0,mesh,1
	                    phi0(s,i)=y(i)
	                end do
	                nlabel(s)=n
	                mlabel(s)=m
                    s=s+1
                end if
            end do
        end do
        !$OMP END PARALLEL DO
    
close(12)    
!close(13)
close(14)
close(15)
close(16)
close(17)
close(20)
close(21)
deallocate ( r, sqr, r2, y )
!deallocate (Hnm, Enm, Htemp, HtempB, e0, eP0, phi0)
deallocate (Hnm, Enm, Htemp, e0, eP0, phi0)
!pause
end program magneticexciton
!
!--------------------------------------------------------------------
subroutine do_mesh ( xmin, dx, mesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  use materials
  implicit none
  integer, intent (in) :: mesh
  real (dp), intent (in) :: xmin, dx
  real (dp), intent (out) :: r(0:mesh), sqr(0:mesh), r2(0:mesh)
  !
  integer :: i
  real(dp) :: x
  !
  do i=0,mesh
     x = xmin + dx * i
     r(i)  = exp(x)
     sqr(i)= sqrt(r(i))
     r2(i) = r(i) * r(i)
  end do
  return
end subroutine do_mesh

!---------------------------------------------------------------------
subroutine solve_sheq ( n, m, e, e_const, Bb, mesh, dx, r, sqr, r2, y, plotnull )
!---------------------------------------------------------------------
  use materials
  
  implicit none
  real(dp), parameter :: eps0=1.0D-12
  real(dp) :: eps
  integer, intent(in) :: mesh, m, n
  !!!The intent(in) attribute of argument i means that i cannot be changed inside the function.
  real(dp), intent(in) :: Bb, dx, r(0:mesh), sqr(0:mesh), r2(0:mesh)
  real(dp), intent(out) :: e, e_const, y(0:mesh) !!!output

  integer :: maxiter, fix
  integer, intent(out) :: plotnull

  integer :: i, j, icl, nodes, ncross, kkk
  integer :: nindex
  
  real (dp) :: ddx12, x2l2, ycusp, dfcusp, fac, norm, de
  real (dp) :: f(0:mesh),vpot(0:mesh)
  real (dp) :: eup, elw, etempup, etemplw
  
  real(dp) :: m_2, vb
  
  nindex = n+abs(m)+1
  
  plotnull = 0
  
  m_2=dble(m)**2
  e_const = hbar**2/(2*mu*me)/((1/1e10)**2)/ElectronCharge

  vb = ElectronCharge / (8*mu*me) * Bb**2 * (1.0/1e10)**2
  vpot(:)= vb*r2(:)
  
  ddx12=dx*dx/12.0_dp
  x2l2 = 2*abs(m)+1.0_dp
  eps = eps0

  eup = vpot(mesh)/e_const
  elw = minval(m_2/r2(:)+vpot(:)/e_const)
  etempup = eup
  etemplw = elw

  e = 0.5_dp * (elw + eup)
  
  maxiter=200
  if ( m==0 ) then
      maxiter=300
  end if
  
  fix = 0
100 do kkk = 1,maxiter
   
     icl = -1
     do i=0,mesh,1
         f(i) = ddx12 * ( m_2 + r2(i) * (vpot(i) /e_const - e) )
     end do
     
     do i=1,mesh   
        if( f(i) == 0.0_dp ) f(i)=1.d-20
        if( f(i) /= sign( f(i), f(i-1)) ) icl=i
     end do  
     
     ! f function as required by numerov method
     f = 1.0_dp - f
     y = 0
     
     nodes = n
     y(0) = r(0)**abs(m)
     y(1) = r(1)**abs(m)
     !y(0) = Rnm(0)
     !y(1) = Rnm(1)
     
     ncross=0
     do i = 1, icl-1
        y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
        if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
     end do   
     
     if (icl < 0 .AND. kkk == maxiter .AND. fix == 0) then
         fix = 1
         eup = etempup
         elw = etemplw
         e = (elw + eup)/2
         goto 100
     end if
     
     if (kkk < maxiter) then
         if ( icl < 0 .or. icl >= mesh-2 ) then
            if (fix == 0) then
                eup = e
                e = 0.5_dp * (eup+elw)
                cycle
            end if
            if (fix == 1) then
                elw = e
                e = 0.5_dp * (eup+elw)
                cycle
            end if      
         end if
     end if
     
     if (icl < 0) then
         plotnull = 1
         write(*,*)'null calculation at (nm)=','(',n,',',m,')'
     end if
     if (kkk == maxiter)exit
      
     fac = y(icl)
     
     if ( ncross /= nodes ) then
        if ( ncross > nodes ) then
           eup = e
        else
           elw = e
        end if
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     
     y(mesh) = dx
     y(mesh-1) = (12.0_dp-10.0_dp*f(mesh))*y(mesh)/f(mesh-1)
     
     ! inward integration 
     do i = mesh-1, icl+1, -1
        y(i-1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)       
        if (y(i-1) > 1.0d10) then
           do j=mesh,i-1,-1
              y(j) = y(j)/y(i-1)
           end do
        end if
     end do
     
     fac = fac/y(icl)
     y (icl:) = y(icl:)*fac
     
     !do i=0,mesh
     !    if (abs(y(i)) < 1.0D-200) then
     !       y(i)=0
     !    end if
     !end do
     norm = 0.d0
     do i=0,mesh-1
         norm = norm + y(i)*y(i) * r2(i) * dx
     end do
     norm = sqrt(norm)
     y = y / norm

     i = icl
     ycusp = (y(i-1)*f(i-1)+f(i+1)*y(i+1)+10.0_dp*f(i)*y(i)) / 12.0_dp
     dfcusp = f(i)*(y(i)/ycusp - 1.0_dp)

     de = dfcusp/ddx12 * ycusp*ycusp * dx !multiplied by 1/constant of E in f(x) (regardless of 1-f and ddx12)
     if (de > 0.0_dp) elw = e
     if (de < 0.0_dp) eup = e

     e = max( min (e+de,eup),elw)

     if (abs(de) < eps) exit    
  end do
  ! ---- convergence has been achieved -----
     write (14,'("n=",i2,"m=",i2, "convergence achieved at iter #",i4," de = ", e11.4)') n,m,kkk,de
  !-------------------------------------------
     
     if (nindex <= 3) then  
        do i=0,mesh
            write(20,'(3e20.10)')r(i),y(i)
        end do
        write(20,*)
     end if
  return
end subroutine solve_sheq