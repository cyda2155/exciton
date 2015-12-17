!---------------------------------------------------------------             
! solving the radial schroedinger equation by Numerov method
! Atomic (Ry) units
! by Y.D.Chen, 2015.11.03
!---------------------------------------------------------------
module materials
!---------------------------------------------------------------
    integer, parameter :: dp = selected_real_kind(18,250)
    real (dp) ,parameter :: Z=23.474178403755868_dp,B0=12.356318010184221_dp
    real (dp) ,parameter :: En0=0.001395077869591_dp
    real (dp) ,parameter :: delta=1.66_dp,lambda=0.075_dp
    
    integer,parameter :: tau = 1 !+1 corresponds to K valley while -1 to -K valley
    integer,parameter :: ms = 1 !ms is the spin projection, tau is used to index the valleys
    integer,parameter :: m = 2
    integer,parameter :: mr = -1
    real (dp),parameter :: PI=3.41592653_dp
    real (dp) :: eg
    contains
        subroutine eg_initial
        eg=delta-dble(ms*tau)*lambda
        end subroutine
end module
!
!---------------------------------------------------------------
program magneticexciton
!---------------------------------------------------------------
use materials

implicit none
!!!SELECTED_REAL_KIND(P,R) returns the kind value of a real data type with decimal precision of at least P digits, exponent range of at least R, and with a radix of RADIX.
  
integer mesh
integer meshold
integer n, i
integer s
integer,parameter :: nn=3, nr=2
integer b
integer,parameter :: bN = 1000
integer :: plotnull
integer :: fix=1
  
real (dp) :: zeta, zmesh, rmax, xmin, dx, e
real (dp) :: er
real (dp) :: rmaxold, dxold
real (dp), allocatable :: r(:), sqr(:), r2(:), y(:), vpot(:)
real (dp) :: Bb, Bbmax, Bbmin, dBb, beta
!real (dp) :: e0, egood, ebad
  
!character(len = 2) :: cTemp

call eg_initial
    ! initialize atomic charge (Z)
    !
    zeta=Z/2.0_dp
    ! initialize logarithmic mesh
    !
    zmesh=  zeta 
    rmax = 1000.0_dp
    xmin = -8.0_dp
    dx  =  0.01_dp
    rmaxold = rmax
    dxold = dx
    !
    mesh = (log(zmesh*rmax) - xmin)/dx
    meshold = mesh
    allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh) )
    !
    call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
    !
    ! initialize the potential
    !
    call init_pot ( zeta, r, mesh, vpot)
    !
    !initialize interval of Bb and etc.
    Bbmin= 0.0_dp
    Bbmax= 300.0_dp
    dBb= (Bbmax-Bbmin)/bN
    open(unit=15,file='error.dat')
    open(unit=12, file='En.dat') !output a list of En(m)
    do n=1,nn,1 ! read number of nodes (stop if nodes < 0)
        do s=0,nr-1
            write(*,*)m,n,s
            ! scan magnetic field
            Bb=Bbmin
            do b=0,bN,1
                beta=Bb/B0
                ! solve the schroedinger equation in radial coordinates by Numerov method
                call solve_sheq (n, e, beta, mesh, dx, r, sqr, r2, vpot, zeta, y, plotnull, fix )
                !if (b == 0) then
                !    e0 = e
                !end if
                !if ( plotnull == 0 ) then
                !    egood = e
                !end if
                if ( plotnull == 1 ) then
                    !ebad = e
                    !if (abs(ebad-egood)>abs(egood-e0+0.01*e0)) then
                        !rescale
    200                 deallocate ( r, sqr, r2, vpot, y )
                        rmax = rmax+500.0_dp
                        mesh = (log(zmesh*rmax) - xmin)/dx
                        write(15,*) 'mesh rescaled, mesh=',mesh
                        allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh) )
                        call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
                        call init_pot ( zeta, r, mesh, vpot)
                        call solve_sheq (n, e, beta, mesh, dx, r, sqr, r2, vpot, zeta, y, plotnull, fix )
                        fix=fix+1
                        if ( plotnull == 1 .AND. fix<=5 ) then
                            !dx=dx-dx/2.0_dp
                            write(15,*) 'fix times=',fix
                            goto 200
                        end if
                        !recover
                        deallocate ( r, sqr, r2, vpot, y )
                        rmax = rmaxold
                        mesh = meshold
                        dx = dxold
                        allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh) )
                        call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
                        call init_pot ( zeta, r, mesh, vpot)
                        !!if problem not solved, give up the point
                        if ( fix >2 ) then
                            Bb=Bb+dBb
                            fix=1
                            cycle
                        end if
                        fix=1
                    !end if
                end if
                ! write out the eigenvalue energy
                er=beta*(2*dble(s)+abs(dble(mr))+dble(tau)*dble(mr)+1.0_dp)
                write (12,'(e20.10,e20.10)') Bb, En0*(e+er)+eg
                !automatic pulsing or minusing 
                Bb=Bb+dBb
            end do
            write (12,*) 
            !
        end do
    end do
close(12)
close(14)
close(15)
close(16)
deallocate ( r, sqr, r2, vpot, y )
end program magneticexciton
!
!--------------------------------------------------------------------
subroutine do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  use materials
  implicit none
  integer, intent (in) :: mesh
  real (dp), intent (in) :: zmesh, xmin, dx
  real (dp), intent (out) :: r(0:mesh), sqr(0:mesh), r2(0:mesh)
  !
  integer :: i
  real(dp) :: x
  !
  do i=0,mesh
     x = xmin + dx * i
     r(i)  = exp(x)/zmesh
     sqr(i)= sqrt(r(i))
     r2(i) = r(i) * r(i)
  end do
  !write(*,'(/" radial grid information:")')
  !write(*,'(" dx =",f10.6,", xmin =",f10.6,", zmesh =",f10.6)') &
  !     dx,xmin,zmesh
  !write(*,'(" mesh =",i6,", r(0) =",f10.6,", r(mesh) =",f10.6)') &
  !     mesh, r(0), r(mesh)
  !write(*,*) 
  !
  return
end subroutine do_mesh
!---------------------------------------------------------------------
subroutine solve_sheq ( n, e, beta, mesh, dx, r, sqr, r2, vpot, zeta, y, plotnull, fix )
  !---------------------------------------------------------------------
  !
  ! solve the schroedinger equation in radial coordinates on a 
  ! logarithmic grid by Numerov method - atomic (Ry) units
  !
  use materials
  
  implicit none
  real(dp), parameter :: eps0=1.0D-12
  real(dp) :: eps
  integer, intent(in) :: mesh, n, fix
  !!!The intent(in) attribute of argument i means that i cannot be changed inside the function.
  real(dp), intent(in) :: beta, dx, r(0:mesh), sqr(0:mesh), r2(0:mesh), & 
       vpot(0:mesh), zeta
  real(dp), intent(out) :: e, y(0:mesh) !!!output
  integer, intent(out) :: plotnull
  
  integer :: maxiter, epsfix
  integer, parameter :: epsfixtor = 4
  integer :: i, j, icl, nodes, ncross, kkk
  
  real (dp) :: ddx12, sqlhf, x2l2, ycusp, dfcusp, fac, norm, eup, elw, de
  real (dp) :: f(0:mesh)
  real (dp) :: Amesh,temp,temp1,temp2
  real (dp) :: ytemp_1
  real (dp) :: etempup,etemplw
  
  integer*8,external :: factorial !()!
  integer :: a1,u
  
  open(14,file='convergence.dat')
  open(15,file='error.dat')
  plotnull=0
  epsfix=0
  
  ddx12=dx*dx/12.0_dp
  sqlhf = dble(m)**2.0_dp
  x2l2 = 2*abs(m)+1.0_dp
  eps = eps0
  
  ! set (very rough) initial lower and upper bounds to the eigenvalue
  !eup = vpot(mesh)
  !elw = minval (sqlhf/r2(:)+vpot(:))
  !eup = vpot(mesh)+beta**2.0_dp/4.0_dp*r2(mesh)
  !elw = minval (sqlhf/r2(:)+vpot(:))
  eup = maxval(vpot(:)+beta**2.0_dp/4.0_dp*r2(:)+dble(tau)*beta*dble(m))
  elw = minval(vpot(:)+sqlhf/r2(:)+dble(tau)*beta*dble(m))
  etempup=eup
  etemplw=elw
  
  if (eup-elw < eps) then
     write (*,*) elw, eup
     write (*,*) 'solve_sheq: lower and upper bounds are equal'
     stop
  end if
  e = 0.5_dp * (elw + eup)
  
  maxiter=200
  if ( m==0 ) then
      maxiter=500
  end if
  
  
100  do kkk = 1,maxiter
     ! set up the f-function and determine the position of its last
     ! change of sign
     ! f < 0 (approximately) means classically allowed region
     ! f > 0 forbidden
     icl = -1
     f(0) = ddx12 *( sqlhf + beta**2.0D0 /4.0D0 * r2(0)**2 + r2(0) * (vpot(0)+dble(tau)*beta*dble(m)-e) )
     
     do i=1,mesh
        f(i) = ddx12 *( sqlhf + beta**2.0D0 /4.0D0 * r2(i)**2 + r2(i) * (vpot(i)+dble(tau)*beta*dble(m)-e) )
        !
        ! beware: if f(i) is exactly zero the change of sign is not observed
        ! the following line is a trick to prevent missing a change of sign 
        ! in this unlikely but not impossible case:
        !
        if( f(i) == 0.0_dp ) f(i)=1.d-20
        if( f(i) /= sign(f(i),f(i-1)) ) icl=i
     end do

     if ( icl < 0 .or. icl >= mesh-2 ) then
        !
        ! classical turning point not found or too far away
        ! no panic: it may follow from a bad choice of eup in
        ! the first iterations. Update e and eup and re-try
        eup = e
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     !
     ! f function as required by numerov method
     f = 1.0_dp - f
     y = 0
     !
     ! determination of the wave-function in the first two points 
     ! (asymptotic behaviour - second term depends upon the potential)
     !
     nodes = n-1
     y(0) = r(0)**(abs(m)) *(1.0_dp + vpot(0)*r2(0)/x2l2)
     y(1) = r(1)**(abs(m)) *(1.0_dp + vpot(1)*r2(1)/x2l2)
     !if ( fix > 5 ) then
     !    y(0) = r(0)**(x2l2/2) / sqr(0)
     !    y(1) = r(1)**(x2l2/2) / sqr(1)
     !    !eup = etempup-(1.0e-5)*etempup
     !    !elw = etemplw+(1.1e-5)*etemplw
     !    !e = 0.5_dp * (eup+elw)
     !end if
     !
     ! outward integration, count number of crossings
     !
     ncross=0
     do i = 1, icl-1
        y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
        if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
     end do
     fac = y(icl) 
     !
     ! check number of crossings
     !
     if ( ncross /= nodes ) then
        if ( ncross > nodes ) then
           eup = e
        else
           elw = e
        end if
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     !
     ! determination of the wave-function in the last two points 
     ! assuming y(mesh+1) = 0
     !
     ! calculate Amesh to obtain y(mesh)
     Amesh = 0
     temp1 = ( dble(factorial(n-1)) * dble(factorial(n-1+abs(m))) )**0.5_dp
     ul: do u = 0, n-1
        temp2 = dble(factorial(n-1-u)) * dble(factorial(abs(m)+u)) * dble(factorial(u))
        temp = (-1)**u * temp1 / temp2 * (beta/2.0_dp)**(dble(abs(m)+1)/2.0_dp+dble(u)) * r(mesh)**(dble(abs(m))-0.5_dp+dble(2*u))
        Amesh = Amesh + temp
     end do ul  
     y(mesh) =  PI**(-0.5_dp) * exp( -beta /4.0_dp * r(mesh)**2.0_dp) * Amesh
     
     Amesh = 0
     temp1 = ( dble(factorial(n-1)) * dble(factorial(n-1+abs(m))) )**0.5_dp
     do u = 0, n-1
        temp2 = dble(factorial(n-1-u)) * dble(factorial(abs(m)+u)) * dble(factorial(u))
        temp = (-1)**u * temp1 / temp2 * (beta/2.0_dp)**(dble(abs(m)+1)/2.0_dp+dble(u)) * r(mesh-1)**(dble(abs(m))-0.5_dp+dble(2*u))
        Amesh = Amesh + temp
     end do
     y(mesh-1) = dx + PI**(-0.5_dp) * exp( -beta /4.0_dp * r(mesh-1)**2.0_dp) * Amesh
     !y(mesh-1) = (12.0_dp-10.0_dp*f(mesh))*y(mesh)/f(mesh-1)
     !
     ! inward integration 
     !
     do i = mesh-1, icl+1, -1
        y(i-1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)
         
        if (y(i-1) > 1.0d10) then
           do j=mesh,i-1,-1
              y(j) = y(j)/y(i-1)
           end do
        end if
     end do
     !
     ! rescale function to match at the classical turning point (icl)
     !
     fac = fac/y(icl)
     y (icl:) = y(icl:)*fac
     !
     ! normalization - note the change of variable:
     !  \int f(r)dr => \sum_i f_i r_i Delta x
     !
     norm = 0.d0
     do i=1,mesh
         norm = norm + y(i)*y(i) * r2(i) * dx
     end do
     norm = sqrt(norm)
     y = y / norm
     !
     ! find the value of the cusp at the matching point (icl)
     !
     i = icl
     ycusp = (y(i-1)*f(i-1)+f(i+1)*y(i+1)+10.0_dp*f(i)*y(i)) / 12.0_dp
     dfcusp = f(i)*(y(i)/ycusp - 1.0_dp)
     !
     ! eigenvalue update using perturbation theory
     !
     de = dfcusp/ddx12 * ycusp*ycusp * dx 
     if (de > 0.0_dp) elw = e
     if (de < 0.0_dp) eup = e
     !
     ! prevent e to go out of bounds, i.e. e > eup or e < elw 
     ! (might happen far from convergence)
     !
     e = max( min (e+de,eup),elw)
     !
     ! convergence check
     !
     if (abs(de) < eps) exit
     if (kkk >= maxiter/2 .AND. epsfix <= epsfixtor ) then
         eps = eps*10
         eup = etempup
         elw = etemplw
         e = 0.5_dp * (eup+elw)
         epsfix = epsfix+1
         goto 100
     end if
     if ( epsfix > epsfixtor ) then
         plotnull=1
         write(15,*)'wrong occur,''elw=',elw,'eup=',eup,'e=',e
         exit
     else
         plotnull=0
     end if
     !
  end do
  ! ---- convergence has been achieved -----
     write (14,'("n=",i2,"m=",i2,"beta=",f6.3," convergence achieved at iter #",i4," de = ", e11.4)') n,m,beta,kkk,de
  !--------------------------------------------
  return
end subroutine solve_sheq
!
!--------------------------------------------------------------------
subroutine init_pot ( zeta, r, mesh, vpot )
!--------------------------------------------------------------------
  !
  ! initialize potential
  !
  use materials
  
  implicit none
  integer, intent (in) :: mesh
  real (dp), intent(in) :: zeta, r(0:mesh)
  real (dp), intent(out):: vpot(0:mesh)
  integer :: i

  open (16,file='pot.dat',status='unknown',form='formatted')
  do i =0,mesh
     vpot (i) = - 2.0_dp*zeta/r(i)
     write (16,*) r(i),vpot(i)
  end do
  write (16,*)
  return
end subroutine init_pot
    
function factorial(n) !()!
    integer :: n,i
    integer*8 :: factorial
    factorial=1
    do i=1,n
        factorial=factorial*i
    end do
    return
end function
