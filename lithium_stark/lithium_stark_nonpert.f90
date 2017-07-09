!---------------------------------------------------------------             
! solving hydrogenic Stark effect by Numerov method
! Atomic (Ry) units
! by Y.D.Chen, 2015.11.03
!---------------------------------------------------------------
module materials
!---------------------------------------------------------------

    integer, parameter :: dp = selected_real_kind(18,250)
    !integer, parameter :: dp = real *8
    
    real (dp), parameter :: aB = 0.52917721092_dp !in Angstrom
    real (dp), parameter :: E1 = -13.60569253_dp !in eV
    
    real (dp), parameter :: PI = 3.141592653_dp
    real (dp), parameter :: ElectronCharge = 1.60217657000000e-19
    real (dp), parameter :: hbar = 1.05457172600000e-34
    
    !Coulumb potential
    real (dp), parameter :: epsilon0 = 8.854187817620e-12
    real (dp), parameter :: Z = 3.0_dp
    real (dp), parameter :: au2ev = 27.2107_dp
    real (dp), parameter :: Eau3s = -0.07418_dp
    real (dp), parameter :: Ea3s = Eau3s * au2ev !3s energy in eV
    
    integer :: tau=1 !+1 corresponds to K valley while -1 to -K valley
    integer :: ms=1 !ms is the spin projection, tau is used to index the valleys
    
    real (dp), parameter :: me = 9.109383560000000e-31
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
    integer n, l, i, s, nl, m, nindex
    
    !label n,l in 1 dimension
    integer :: nm
    integer, parameter :: Nnindex = 17
    integer, parameter :: Nnm = Nnindex * (Nnindex + 1) * (2*Nnindex + 1) / 6 - 1

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
    integer :: s1, s2, lm1, lm2, nindex1, nindex2, nl1, nl2, n1, l1, n2, l2, m1, m2
    integer, allocatable :: nlabel(:), llabel(:), mlabel(:)
    
    !double precision, allocatable :: Hnm(:,:)
    !double precision, allocatable :: Htemp(:,:), Htempv(:,:) !store off-diagonal integrated elements
    !double precision, allocatable :: Enm(:,:)
    real*8, allocatable :: Hnm(:,:)
    real*8, allocatable :: Htemp(:,:)!, Htempv(:,:) !store off-diagonal integrated elements
    real*8, allocatable :: Enm(:,:)
    real (dp) :: Hprime, Hc!, Hprimev, Hcv
    real (dp) :: p_int(Nnindex**2, Nnindex**2), pmn_int
    
    real (dp) :: phi0_int((1+Nnindex)*Nnindex/2-1,(1+Nnindex)*Nnindex/2-1)
    
    !integer,parameter :: LDA = Nnm
    !integer,parameter :: LWMAX = 200*Nnm
    !integer :: info, lwork
    !double precision :: work(LWMAX) 
    
    ! rearrange and output
    integer :: eigenindex, eigenindex_new, eigenindex_lw
    real*8 :: Eave, Eavelw, Hswap
    integer, parameter :: ntot = Nnindex
    !integer, parameter :: nout = ntot * (1+ntot) / 2 !number of output energy levels
    integer, parameter :: nout = ntot * (ntot + 1) * (2*ntot + 1) / 6 - 1
    real (dp) :: ediff
    
    ! skip point
    integer :: plotnull = 0
    
    ! initialize interval of Electric field and etc.
    integer, parameter :: EfN = 1e3
    integer Efindex
    real (dp) :: Ef, Efmin, Efmax, dEf
    
    real (dp) :: a0, a
    real (dp) :: alw, aup, da, ealw, eaup
    real (dp), parameter :: daeps = 1e-6
    integer :: cycletimes
    real (dp) :: EffectiveCharge
    
    ! output filename variable
    character(len = 4) :: cTemp, cTemp1, cTempout

    ! Electric field setup
    Efmin= 0.0D0
    Efmax= 5e5
    dEf= (Efmax-Efmin)/EfN
    ! allocate logarithmic mesh
    !
    xmin = -8.0_dp
    mesh = 3e3 * Nnindex / 2
    !revise mesh
    a0 = aB
    
    !initialize mesh to guarantee the interval of wavefunction is sufficient
    rmax = 5*dble(Nnindex - 0.5)**2 * (Nnindex) + 50 * dble(Nnindex - 0.5)
    
    write(*,*)"a0=",a0
    
    allocate (e0(Nnm))
    allocate (eP0(Nnm))
    allocate (phi0((1+Nnindex)*Nnindex/2-1,0:mesh))
    allocate (nlabel(Nnm), llabel(Nnm), mlabel(Nnm))
    allocate (r(0:mesh), sqr(0:mesh), r2(0:mesh), y(0:mesh))
    allocate (Hnm(Nnm,Nnm))
    allocate (Htemp(Nnm,Nnm))!,Htempv(Nnm,Nnm)
    allocate (Enm(EfN+1,Nnm))
    !initialization
    e0=0
    phi0=0
    Hnm=0
    Htemp=0
    !Htempv=0
    Enm=0
    r=0
    sqr=0
    r2=0
    y=0
    
    dx = (log(rmax) - xmin)/dble(mesh)
    call do_mesh ( xmin, dx, mesh, r, sqr, r2 )
    
    !determine a
    !evaluate by 3s energy
    n=2
    l=0
    a = ElectronCharge**2 * me / (hbar**2) * Z
    da = 1.0_dp
    call solve_sheq ( n, l, e, e_const, mesh, a0, a, dx, r, sqr, r2, y, plotnull )
    if ((e*e_const) < Ea3s) then
        alw = a
        aup = a + da
    else
        aup = a
        alw = a - da
    end if
    cycletimes = 1
    do
        call solve_sheq ( n, l, ealw, e_const, mesh, a0, alw, dx, r, sqr, r2, y, plotnull )  
        call solve_sheq ( n, l, eaup, e_const, mesh, a0, aup, dx, r, sqr, r2, y, plotnull )
        
        if (alw < da) then
            da = alw /10
        end if
        if ((eaup*e_const) < Ea3s) then
            alw = alw + da
            aup = aup + da
        else if ((ealw*e_const) > Ea3s) then
            alw = alw - da
            aup = aup - da
        else if ((eaup*e_const) >= Ea3s .and. (ealw*e_const) <= Ea3s) then
            da = da/2
            aup = alw + da
        else
            write(*,*)"error finding correct a"
            stop
        end if
        a = 0.5 * (alw + aup)
        !write(*,*)"r0=",r0
        if (da < daeps) then
            write(*,*)"a=",a
            exit
        end if
        if (cycletimes > 5e3) then
            write(*,*)"error finding correct a"
            stop
        end if
        cycletimes = cycletimes + 1
    end do
    EffectiveCharge = a / (ElectronCharge**2 * me / (hbar**2))
    write(*,*)"Screened core charge=",EffectiveCharge
    
    !calculate H0
    e0=0
    phi0=0
    open(unit=15,file='error.dat')
    open(unit=14,file='convergence.dat')
    open(unit=12, file='Orthonormal.dat')
    open(unit=13, file='Perpendicular.dat')
    open(unit=16, file='Enm.dat')
    open(unit=20, file='phicompare.dat')
    open(unit=21, file='enmcompare.dat')
    !$OMP PARALLEL DO
    s=1
    do nindex=2,Nnindex,1  
        do l=0,nindex-1,1 ! read number of nodes (stop if nodes < 0)
            write(cTemp,'(i4)') int(l) !label the file name to distinguish different B   !(f3.0)
            n = nindex-l-1
            
            nl = nindex * (nindex - 1) / 2 + l
            
            write(*,*)"n=",n,"l=",l, "mesh=",mesh
            write(cTemp1,'(i4)') int(n)
            !open(unit=13, file='Phinm'//trim(adjustl(cTemp1))//trim(adjustl(cTemp))//'.dat')
                
            call solve_sheq ( n, l, e, e_const, mesh, a0, a, dx, r, sqr, r2, y, plotnull )
            
            !do i=0,mesh,1
            !    write(13,'(2e20.10)')r(i),y(i)
            !end do
	    
	        do i=0,mesh,1
            	phi0(nl,i)=y(i)
            end do
            write (16,'(i10,i10,e20.10)')n,l,e*e_const

            if (plotnull == 1) then
                write(*,*)"calculating error, n=",n," l=",l
            end if
            
            do m = -l,l,1
                e0(s)=e*e_const
                nlabel(s)=n
                llabel(s)=l
                mlabel(s)=m
                s=s+1
            end do
        end do
    end do
    !$OMP END PARALLEL DO
        
    !check if orthonormal        
    do nindex1=2,Nnindex,1  
        do l1=0,nindex1-1,1
            do nindex2=2,Nnindex,1  
                do l2=0,nindex2-1,1
                    n1 = nindex1-l1-1
                    n2 = nindex2-l2-1
                    nl1 = nindex1 * (nindex1 - 1) / 2 + l1
                    nl2 = nindex2 * (nindex2 - 1) / 2 + l2
                    perp = 0
                    norm = 0
                    if (nindex1==nindex2 .and. l1==l2) then
                        do i=0,mesh,1
                            norm = norm + phi0(nl1,i)*phi0(nl2,i)*r2(i)*dx
                        end do
                        write (12,'(2I5,2e20.10)')n1,l1,norm
                    else
                        do i=0,mesh,1
                            perp = perp + phi0(nl1,i)*phi0(nl2,i)*r2(i)*dx
                        end do
                        write (13,'(4I5,2e20.10)')n1,l1,n2,l2,perp  
                    end if  
                end do
            end do
        end do
    end do
    
    write(*,*)"=====p_int calculation begin====="
    !$OMP PARALLEL DO
    p_int = 0
    do l1 = 0, Nnindex - 1
        do m1 = -l1, 0, 1
            do l2 = 0, Nnindex - 1
                do m2 = -l2, 0, 1
                    if (m1 == m2 .and. abs(l1-l2)==1) then
                        call legendre_int (l1, l2, m1, m2, pmn_int)
                        lm1 = l1**2 + m1 + l1 + 1
                        lm2 = l2**2 + m2 + l2 + 1
                        p_int(lm1,lm2) = pmn_int
                        write(*,*)"(l1,m1,l2,m2)=","(",l1,",",m1,",",l2,",",m2,")"!,"p_int=",p_int(lm1,lm2)
                    end if
                end do
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    write(*,*)"=====p_int calculation end====="
    
    do l1 = 0, Nnindex - 1
        do m1 = -l1, l1, 1
            do l2 = 0, Nnindex - 1
                do m2 = -l2, l2, 1
                    lm1 = l1**2 + m1 + l1 + 1
                    lm2 = l2**2 + m2 + l2 + 1
                    p_int(lm1,lm2) = p_int(l1**2 - abs(m1) + l1 + 1, l2**2 - abs(m2) + l2 + 1)
                end do
            end do
        end do
    end do
    
    write(*,*)"=====phi0*r3_int calculation begin====="
    !$OMP PARALLEL DO
    phi0_int = 0
    do nl1 = 1,(1+Nnindex)*Nnindex/2-1
        do nl2 = 1, (1+Nnindex)*Nnindex/2-1
            Hprime=0
            do i=0,mesh,1
                Hprime = Hprime + phi0(nl1,i)*phi0(nl2,i) * (r(i))**3.0 * dx
                !write(*,*) Hprime
            end do
            phi0_int(nl1, nl2) =  Hprime
        end do
    end do
    !$OMP END PARALLEL DO
    write(*,*)"=====phi0*r3_int calculation end====="
                    
    Htemp=0
    !$OMP PARALLEL DO
    do s1 = 1, Nnm
        do s2 = 1, Nnm
            n1 = nlabel(s1)
            l1 = llabel(s1)
            m1 = mlabel(s1)
            n2 = nlabel(s2)
            l2 = llabel(s2)
            m2 = mlabel(s2)
            nindex1 = n1 + l1 + 1
            nindex2 = n2 + l2 + 1
            if (abs(l1-l2)==1) then
                if (m1 == m2) then
                    nl1 = nindex1 * (nindex1 - 1) / 2 + l1
                    nl2 = nindex2 * (nindex2 - 1) / 2 + l2
                    lm1 = l1**2 + m1 + l1 + 1
                    lm2 = l2**2 + m2 + l2 + 1
                    Htemp(s1,s2) = phi0_int(nl1, nl2) * p_int(lm1,lm2)
                end if
    !	    else if (l1==l2) then
    !	    	if (m1 == m2) then
    !		    Hprimev = 0
    !		    do i=0,mesh,1
    !                	Hprimev = Hprimev + phi0(s1,i)*phi0(s2,i) * r(i) * exp(-alpha*r(i)) * dx
    !		    end do
    !		    Htempv(s1,s2) = Hprimev
            end if
        end do
    end do
    !$OMP END PARALLEL DO
    
    !open(unit=13, file='Hnm.dat')
    Enm=0
    !$OMP PARALLEL DO
    !calculate matrix elements of H for each Ef
    do Efindex=0,EfN,1
        Ef = Efmin + Efindex*dEf
        Hnm=0
        Hc = Ef*(1e-10)*a0
	    !Hcv = 2*aB*E1/a0 * zeta
        do s1 = 1, Nnm
            do s2 = 1, Nnm
                n1 = nlabel(s1)
                l1 = llabel(s1)
                m1 = mlabel(s1)
                n2 = nlabel(s2)
                l2 = llabel(s2)
                m2 = mlabel(s2)
                !diagonal block matrix consideration
                if (m1 == m2) then
                    if (abs(l1-l2)==0) then
		                    !Hnm(s1,s2) = Hnm(s1,s2) + Htempv(s1,s2) * Hcv
                        if (abs(n1-n2)==0) then
                            Hnm(s1,s1)= Hnm(s1,s1) + e0(s1)!+(hbar*Ef*1e10)**2/(2*me*mcm)/ElectronCharge
                        end if
                    !off-diagonal block matrix consideration
                    else if (abs(l1-l2)==1) then
                        Hnm(s1,s2) = Htemp(s1,s2) * Hc
                    end if
                end if
                !if (Efindex == 1) write(13,'(e20.10"    "\)')Hnm(s1,s2)
            end do
            !if (Efindex == 1) write(13,*)
        end do
            
        !lwork=-1
        !call dsyev('V','U',Nnm,Hnm,LDA,Enm(Efindex+1,:),work,lwork,info)
        !lwork = min(LWMAX, int(work(1)))
        !call dsyev('V','U',Nnm,Hnm,LDA,Enm(Efindex+1,:),work,lwork,info)
        call syev(Hnm,Enm(Efindex+1,:))
        
        if (Efindex==0) eP0(:)=Enm(Efindex+1,:)
        !output different energy level with the same Ef
    end do !momentum Efindex loop
    !$OMP END PARALLEL DO
                   
    !output
    open(unit=17, file='Eout.dat')
    do s=2,nout,1
        do Efindex=0,EfN,1
            Ef = Efmin + Efindex*dEf
            write(17,'(3e20.10)')Ef,Enm(Efindex+1,s),Enm(Efindex+1,s)-eP0(s)
        end do
        write(17,*)
    end do
    
close(12)    
!close(13)
close(14)
close(15)
close(16)
close(17)
close(20)
close(21)
deallocate ( r, sqr, r2, y )
deallocate (Hnm, Enm, Htemp, e0, eP0, phi0)
deallocate (nlabel, llabel, mlabel)
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
subroutine solve_sheq ( n, l, e, e_const, mesh, a0, a, dx, r, sqr, r2, y, plotnull )
!---------------------------------------------------------------------
  use materials
  
  implicit none
  real(dp), parameter :: eps0=1.0D-12
  real(dp) :: eps
  integer, intent(in) :: mesh, l, n
  !!!The intent(in) attribute of argument i means that i cannot be changed inside the function.
  real(dp), intent(in) :: a0, a, dx, r(0:mesh), sqr(0:mesh), r2(0:mesh)
  real(dp), intent(out) :: e, e_const, y(0:mesh) !!!output

  integer :: maxiter, fix
  integer, intent(out) :: plotnull

  integer :: i, j, icl, nodes, ncross, kkk
  
  real (dp) :: ddx12, x2l2, ycusp, dfcusp, fac, norm, de
  real (dp) :: f(0:mesh),vpot(0:mesh)
  real (dp) :: eup, elw, etempup, etemplw
  
  !real(dp), external :: factorial
  !real*8 :: confluentF(0:mesh)
  !real(dp) :: const
  !real(dp) :: Rnm(0:mesh)
  !real(dp) :: e_compare
  
  real(dp) :: m_2, va
  
    real(dp), parameter :: zeta = Z - 1.0_dp
    !real(dp), parameter :: a = 2.13_dp
    
    real(dp) :: alpha
    
    alpha = a0 * a
  
  !nindex = n+l+1
  !const = 2 * mu / epsilon / aB /(nindex - 0.50_dp) * a0
  !e_compare = E1 * mu / (epsilon*(nindex-0.50_dp))**2
  !call hypergeoF(const,r,l,nindex,mesh,confluentF)
  !do i=0,mesh
  !    Rnm(i) = const/factorial(2*l)*sqrt(factorial(nindex+l-1)/factorial(nindex-l-1)/(2*nindex-1))*(const*r(i))**l*exp(-const*r(i)/2)*confluentF(i)
  !end do
  
  plotnull = 0
  
  va = 2*aB*E1/a0
  
  m_2= (l+0.5_dp)**2
  e_const = hbar**2/(2*me)/((a0/1e10)**2)/ElectronCharge
  
  !vpot(:)=va/r(:) + va * zeta / r(:) * exp(-alpha*r(:))
  vpot(:)=va/r(:) + va*zeta/r(:) * (1 + alpha*a0*r(:)) * exp(-2*alpha*a0*r(:))
  !vpot(:)= va * zeta / r(:) * (1 + alpha * r(:)) * exp(-2*alpha*r(:))
  
  ddx12=dx*dx/12.0_dp
  x2l2 = 2*l+2.0_dp
  eps = eps0
  
  !e_landau = beta*(2*dble(n)+abs(dble(l))+1)

  eup = vpot(mesh)/e_const
  elw = minval(m_2/r2(:)+vpot(:)/e_const)
  etempup = eup
  etemplw = elw

  e = 0.5_dp * (elw + eup)
  
  maxiter=200
  if ( l==0 ) then
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
     y(0) = r(0)**(l+0.5)*(1.0_dp + Z*va/e_const*r(0)/x2l2)
     y(1) = r(1)**(l+0.5)*(1.0_dp + Z*va/e_const*r(1)/x2l2)     
     !y(0) = r(0)**(l+0.5)*(1.0_dp + zeta*va/e_const*r(0)/x2l2)
     !y(1) = r(1)**(l+0.5)*(1.0_dp + zeta*va/e_const*r(1)/x2l2)
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
         !write(*,*)n,l,e
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
         write(*,*)'null calculation at (nl)=','(',n,',',l,')'
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
     write (14,'("n=",i2,"l=",i2, "convergence achieved at iter #",i4," de = ", e11.4)') n,l,kkk,de
  !-------------------------------------------
     
     !if (n<3 .and. l<2) then  
     !   norm = 0.d0
     !   do i=0,mesh
     !       write(20,'(3e20.10)')r(i),y(i),Rnm(i)
     !       norm = norm + Rnm(i)*Rnm(i) * r2(i) * dx
     !   end do
     !   write(*,*)norm
     !   write(20,*)
     !   write(21,'(2e20.10)')e_compare,e_const*e
     !end if
     !use analytical results for expansion
     !y(:)=Rnm(:)
     !e=e_compare / e_const
  return
end subroutine solve_sheq
    
function factorial(n) !()!

    use big_integer_module

    integer :: n
    real(dp) :: factorial

    integer :: i
    type(big_integer) :: fact
    real*8 :: time_begin,time_end

    !call cpu_time(time_begin)
    fact=1
    do i=1,n
      fact=fact*i
    end do

    factorial=big_integer_to_real(fact)
    if (n==0 .or. n==1) factorial=1

    !call cpu_time(time_end)

    !write(*,'(I6"!="\)') n
    !call print_big(fact)
    !
    !write(*,*)
    !write(*,*) 'computatonal time=',time_end-time_begin,'seconds'
    !write(*,*) 'number of digits=',len_trim(char(fact))
end function
    
!laguerre initialization
subroutine hypergeoF(const,r,l,nindex,mesh,confluentF)
    !laguerre
    use materials
    integer, intent(in) :: l, nindex, mesh
    real(dp), intent(in) :: const, r(0:mesh)
    real*8, intent(out) :: confluentF(0:mesh)
    real(dp), external :: factorial
    
    integer :: i
    integer :: n1,l1
    real*8, allocatable :: v(:,:)
    real*8, allocatable :: x(:)
    
    n1=nindex-l-1
    l1=2*l
    allocate ( v(mesh+1,0:n1) )
    allocate ( x(mesh+1) )
    
    do i=1,mesh+1
        x(i)=r(i-1)*const
    end do
    
    call lm_polynomial ( mesh+1, n1, l1, x, v )
    do i=1,mesh+1
        confluentF(i-1)=v(i,n1)*factorial(n1)*factorial(l1)/factorial(n1+l1)
    end do
    deallocate ( v )
    return
end subroutine
    
subroutine legendre_int (l1, l2, m1, m2, pmn_int)

use materials
implicit none
integer, intent(in) :: l1,l2,m1,m2
real(dp), intent(out) :: pmn_int
real(dp) :: pval_int, eps0, h
real(dp), external :: factorial

integer :: xdivide = 1e5+1
real(dp) :: xlw = -1.0
real(dp) :: xup = 1.0
real(dp) :: x, dx, Pmn1(0:l1), Pmn2(0:l2)
integer :: i

dx = (xup - xlw) / (xdivide-1)

pmn_int = 0
do i = 1, xdivide-1
    x = xlw + (i-1) * dx
    call pmn_polynomial_value( l1, abs(m1), x, Pmn1 )
    call pmn_polynomial_value( l2, abs(m2), x, Pmn2 )
    eps0 = 1.0_dp
    h = dx
    call intout ( x, x+dx, Pmn1(l1)*Pmn2(l2)*x*dx, pval_int, eps0, h, l1, l2, m1, m2 )
    pmn_int = pmn_int + pval_int
    !pmn_int = pmn_int + Pmn1(l1)*Pmn2(l2)*x*dx
end do 

return
end subroutine legendre_int

subroutine pmn_polynomial_value ( n, m, x, cx )

  use materials
  implicit none

  integer ( kind = 4 ) n

  real (dp) cx(0:n)
  real (dp) factor
  integer j
  integer m
  real (dp), external :: factorial
  real (dp) x

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop 1
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop 1
  end if

  cx(0:n) = 0.0D+00

  if ( m <= n ) then
    cx(m) = 1.0D+00
    factor = 1.0D+00
    do j = 1, m
      cx(m) = - cx(m) * factor * sqrt ( 1.0D+00 - x**2 )
      factor = factor + 2.0D+00
    end do
  end if

  if ( m + 1 <= n ) then
    cx(m+1) = x * real ( 2 * m + 1, kind = 8 ) * cx(m)
  end if

  do j = m + 2, n
    cx(j) = ( real ( 2 * j     - 1, kind = 8 ) * x * cx(j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(j-2) ) &
                 / real (     j - m,     kind = 8 )
  end do
!
!  Normalization.
!
  do j = m, n
    factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * factorial ( j - m ) ) &
      / ( 2.0D0 * factorial ( j + m ) ) )
    cx(j) = cx(j) * factor
  end do

  return
end
    
subroutine intout ( xval_lw, xval_up, pval_ini_int, pval_int, eps0, h, l1, l2, m1, m2 )
    use materials
    implicit none
    integer, intent(in) :: l1,l2,m1,m2
    real(dp), intent(in) :: xval_lw, xval_up, pval_ini_int, eps0, h
    real(dp), intent(out) :: pval_int
    real(dp), parameter :: eps = 3e-8
    real(dp) :: x, dx, Pmn1(0:l1), Pmn2(0:l2), inttemp, epstemp, pval_intold
    real*8 :: time_begin,time_end
    
    epstemp = eps0
    dx = h
    pval_intold = pval_ini_int
    
    do
        x = xval_lw
        inttemp = 0
        do
            call pmn_polynomial_value( l1, abs(m1), x, Pmn1 )
            call pmn_polynomial_value( l2, abs(m2), x, Pmn2 )
            inttemp = inttemp + Pmn1(l1)*Pmn2(l2)*x*dx
            x = x + dx
            if (x>=xval_up) exit
        end do
        pval_int = inttemp
        epstemp = abs(pval_int - pval_intold)
        pval_intold = pval_int
        if (epstemp < eps) exit
        dx = dx / 2.0_dp
    end do
    return
    
end subroutine intout   