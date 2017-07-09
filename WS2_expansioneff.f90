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
    integer, parameter :: Nnindex = 20
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
    integer, parameter :: bN = 2 !bN is the total number of Bb
    integer :: b !b is the index of Bb
    real (dp) :: Bb, Bbmax, Bbmin, dBb
    
    ! initialize interval of P and etc.
    integer, parameter :: PN = 5e3
    integer p
    real (dp) :: Pp, Ppmin, Ppmax, dPp
    
    real (dp) :: aB0_ad, a0, l_ad
    real (dp) :: a0min, a0max
    integer :: meshtimes
    
    real (dp) :: r0
    real (dp) :: r0lw, r0up, dr0, er0lw, er0up
    real (dp), parameter :: dr0eps = 1e-10
  
    ! output filename variable
    character(len = 4) :: cTemp, cTemp1, cTempout

    Ppmin= 0.0D0
    Ppmax= 1.0D0/a*pi*0.25
    dPp= (Ppmax-Ppmin)/PN
    ! allocate logarithmic mesh
    !
    xmin = -8.0_dp
    mesh = 2e3 * Nnindex
    
    !initialize interval of Bb and etc.
    Bbmin= 0.0_dp
    Bbmax= 10.0_dp
    dBb= (Bbmax-Bbmin)/bN
    
    !revise mesh
    call constants ( Bbmin, aB0_ad, a0max, l_ad )
    call constants ( Bbmax, aB0_ad, a0min, l_ad )
    a0 = a0min
    meshtimes = int(a0max / a0min)
    if (a0max < a0min) write(*,*)"error"
    !initialize mesh to guarantee the interval of wavefunction is sufficient
    rmax = (5*dble(Nnindex - 0.5)**2 * (Nnindex) + 50 * dble(Nnindex - 0.5))*meshtimes
    mesh = mesh * meshtimes
    
    write(*,*)"a0=",a0
    
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
    do b=0,bN,1
 
        Bb = Bbmin + b*dBb
        write(cTempout,'(i4)') int(b)
        write(*,*)"magnetic field=",Bb,"T"
	
	    r0lw = 0
        r0up = 100
        
        !calculate r0
        n=0
        m=0
        r0=x2D*2*PI
        dr0 = 5
        call solve_sheq (n, m, e, e_const, Bb, r0, mesh, a0, dx, r, sqr, r2, y, plotnull)
        if ((e*e_const+Eg) < Ea1s) then
            r0lw = r0
            r0up = r0 + dr0
        else
            r0up = r0
            r0lw = r0 - dr0
        end if
        do
            call solve_sheq (n, m, er0lw, e_const, Bb, r0lw, mesh, a0, dx, r, sqr, r2, y, plotnull)
            call solve_sheq (n, m, er0up, e_const, Bb, r0up, mesh, a0, dx, r, sqr, r2, y, plotnull)
            if (r0lw < dr0) then
                dr0 = r0lw /10
            end if
            if ((er0up*e_const+Eg) < Ea1s) then
                r0lw = r0lw + dr0
                r0up = r0up + dr0
            else if ((er0lw*e_const+Eg) > Ea1s) then
                r0lw = r0lw - dr0
                r0up = r0up - dr0
            else if ((er0up*e_const+Eg) >= Ea1s .and. (er0lw*e_const+Eg) <= Ea1s) then
                dr0 = dr0/2
                r0up = r0lw + dr0
            else
                write(*,*)"error finding correct r0"
                stop
            end if
            r0 = 0.5 * (r0lw + r0up)
            !write(*,*)"r0=",r0
            if (dr0 < dr0eps) then
                write(*,*)"r0=",r0
                exit
            end if
        end do
        
	    !calculate H0
	    Hnm=0
        Htemp=0
        !HtempB=0
        Enm=0
        y=0
        e0=0
        phi0=0
        open(unit=15,file='error'//trim(adjustl(cTempout))//'.dat')
        open(unit=14,file='convergence'//trim(adjustl(cTempout))//'.dat')
        open(unit=12, file='Orthonormal'//trim(adjustl(cTempout))//'.dat')
        open(unit=13, file='Perpendicular'//trim(adjustl(cTempout))//'.dat')
        open(unit=16, file='Enm'//trim(adjustl(cTempout))//'.dat')
        open(unit=20, file='phicompare'//trim(adjustl(cTempout))//'.dat')
        open(unit=21, file='enmcompare'//trim(adjustl(cTempout))//'.dat')
        !$OMP PARALLEL DO
        s=1
        do nindex=1,Nnindex,1  
            do m_abs=0,nindex-1,1 ! read number of nodes (stop if nodes < 0)
                write(cTemp,'(i4)') int(m_abs) !label the file name to distinguish different B   !(f3.0)
                n = nindex-m_abs-1
                
                write(*,*)"n=",n,"m_abs=",m_abs, "mesh=",mesh
                write(cTemp1,'(i4)') int(n)
                !open(unit=13, file='Phinm'//trim(adjustl(cTemp1))//trim(adjustl(cTemp))//'.dat')
                
                call solve_sheq (n, m_abs, e, e_const, Bb, r0, mesh, a0, dx, r, sqr, r2, y, plotnull)
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

        
        !check if orthonormal        
        do s1 = 1, Nnm
            do s2 = 1, Nnm
                n1 = nlabel(s1)
                m1 = mlabel(s1)
                n2 = nlabel(s2)
                m2 = mlabel(s2)
                perp = 0
                norm = 0

                if (n1==n2 .and. m1==m2) then
                    do i=0,mesh,1
                        norm = norm + phi0(s1,i)*phi0(s2,i)*r2(i)*dx
                    end do
                    write (12,'(2I5,2e20.10)')n1,m1,norm
                else if (m1==m2) then
                    do i=0,mesh,1
                        perp = perp + phi0(s1,i)*phi0(s2,i)*r2(i)*dx
                    end do
                    write(13,'(4I5,2e20.10)')n1,m1,n2,m2,perp
                end if
            
            end do
        end do
           
        Htemp=0
        !HtempB=0
        !$OMP PARALLEL DO
        do s1 = 1, Nnm
            do s2 = 1, Nnm
                n1 = nlabel(s1)
                m1 = mlabel(s1)
                n2 = nlabel(s2)
                m2 = mlabel(s2)
                if (abs(m1-m2)==1) then                  
                    Hprime=0
                    do i=0,mesh,1
                        Hprime = Hprime + phi0(s1,i)*phi0(s2,i) * (r(i))**3.0 * dx
                    end do
                    Htemp(s1,s2) = Hprime             
                !else if (abs(m1-m2)==0) then
                !    HprimeB=0
                !    do i=0,mesh,1
                !        HprimeB = HprimeB + phi0(s1,i)*phi0(s2,i) * (r(i))**4.0 * dx
                !    end do
                !    HtempB(s1,s2) = HprimeB
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        Enm=0
        !$OMP PARALLEL DO
        !calculate matrix elements of H for each Pp
        do p=0,PN,1
            Pp = Ppmin + p*dPp
            Hnm=0
            Hc = hbar*Pp*Bb*a0 / (2*mcm*me) !assume M=4mu
            !HcB = ElectronCharge / (8*mu*me) * Bb**2 * (a0/1e10)**2
            !write(*,*)"Hc=",Hc,"Hcb=",HcB,'a0=',a0
            do s1 = 1, Nnm
                do s2 = 1, Nnm
                    n1 = nlabel(s1)
                    m1 = mlabel(s1)
                    n2 = nlabel(s2)
                    m2 = mlabel(s2)
                    !diagonal block matrix consideration
                    if (abs(m1-m2)==0) then
                        !Hnm(s1,s2) = HtempB(s1,s2) * HcB
                        if (abs(n1-n2)==0) then
                            Hnm(s1,s1)= Hnm(s1,s1) + e0(s1) + hbar2M * Pp**2 / 2.0_dp
                        end if
                    !off-diagonal block matrix consideration
                    else if (abs(m1-m2)==1) then
                        Hnm(s1,s2) = Htemp(s1,s2) * Hc
                    end if
                end do
            end do
            
            call syev(Hnm,Enm(p+1,:))
            !Enm(p+1,:) = Enm(p+1,:) + hbar2M * Pp**2 / 2.0_dp
            if (p==0) eP0(:)=Enm(p+1,:)
            !output different energy level with the same Pp
        end do !momentum p loop
        !$OMP END PARALLEL DO
        
        !do p=1, PN-1, 1
        !    do eigenindex=1,Nnm
        !        Eavelw = abs((Enm(p,eigenindex)+Enm(p+2,eigenindex))/2 - Enm(p+1,eigenindex))
        !        eigenindex_lw = eigenindex
        !        do eigenindex_new=1,Nnm
        !            Eave = abs((Enm(p,eigenindex)+Enm(p+2,eigenindex_new))/2 - Enm(p+1,eigenindex))
        !            if (Eave < Eavelw) then
        !                Eavelw = Eave
        !                eigenindex_lw = eigenindex_new
        !            end if
        !        end do
        !        Hswap=Enm(p+2,eigenindex)
        !        Enm(p+2,eigenindex)=Enm(p+2,eigenindex_lw)
        !        Enm(p+2,eigenindex_lw)=Hswap
        !    end do
        !end do
                    
        !output
        open(unit=17, file='Eout'//trim(adjustl(cTempout))//'.dat')
        do s=1,nout,1
           do p=0,PN,1
              Pp = Ppmin + p*dPp
              write(17,'(3e20.10)')Pp,Enm(p+1,s)+Eg,Enm(p+1,s)-eP0(s)
            end do
            write(17,*)
        end do
                        
    end do !magnetic field loop
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
    
!--------------------------------------------------------------------     
subroutine constants ( Bb, aB0_ad, a0, l_ad)
!--------------------------------------------------------------------
  !
  !setting characteristic length of the system
  !
  use materials  
  implicit none
  real(dp), intent(in) :: Bb
  real(dp), intent(out) :: aB0_ad, a0, l_ad
  integer ,parameter :: NB = 1
  
  l_ad = sqrt(ElectronCharge*Bb/hbar)/1.0e10_dp

  aB0_ad = -E1*aB/(epsilon*hbar2mu)
  a0 = 1.0_dp/(sqrt(aB0_ad**2+l_ad**2)+aB0_ad)

end subroutine constants

!---------------------------------------------------------------------
subroutine solve_sheq (n, m, e, e_const, Bb, r0, mesh, a0, dx, r, sqr, r2, y, plotnull)
!---------------------------------------------------------------------
  use materials
  
  implicit none
  real(dp), parameter :: eps0=1.0D-12
  real(dp) :: eps
  integer, intent(in) :: mesh, m, n
  !!!The intent(in) attribute of argument i means that i cannot be changed inside the function.
  real(dp), intent(in) :: Bb, a0, r0, dx, r(0:mesh), sqr(0:mesh), r2(0:mesh)
  real(dp), intent(out) :: e, e_const, y(0:mesh) !!!output

  integer :: maxiter, fix
  integer, intent(out) :: plotnull

  integer :: i, j, icl, nodes, ncross, kkk
  integer :: nindex
  
  real (dp) :: ddx12, x2l2, ycusp, dfcusp, fac, norm, de
  real (dp) :: f(0:mesh),vpot(0:mesh)
  real (dp) :: eup, elw, etempup, etemplw
  

  
  real(dp) :: m_2, va, vb
  
  nindex = n+abs(m)+1
  
  plotnull = 0
  

  m_2=dble(m)**2
  e_const = hbar**2/(2*mu*me)/((a0/1e10)**2)/ElectronCharge
  
  va = -2*aB*E1/r0
  vb = ElectronCharge / (8*mu*me) * Bb**2 * (a0/1e10)**2
  do i =0,mesh
     vpot(i)=va*(log(a0*r(i)/(a0*r(i)+r0))+gamma*exp(-a0*r(i)/r0)) + vb*r2(i)
  end do
  
  ddx12=dx*dx/12.0_dp
  x2l2 = 2*abs(m)+1.0_dp
  eps = eps0

  eup = vpot(mesh)/e_const
  elw = minval(m_2/r2(:)+vpot(:)/e_const)
  etempup=eup
  etemplw=elw

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
  return
end subroutine solve_sheq
    
