!---------------------------------------------------------------
! Find a solution with given n, m for an hydrogenic atom                
! solving the radial schroedinger equation by Numerov method
! Atomic (Ry) units
! by Y.D.Chen, 2015.11.03
!---------------------------------------------------------------
module materials
!---------------------------------------------------------------
    !material determined parameters
    integer, parameter :: dp = selected_real_kind(14,200)
    real*8,parameter :: En0=13.950778695908895_dp
    real*8,parameter :: Z=0.234741784037559_dp
    
    real (dp) ,parameter :: ratio=3.197132292564991_dp
    real (dp) ,parameter :: V0=0.312780300748119_dp
    integer :: tau=1 !+1 corresponds to K valley while -1 to -K valley
    integer :: ms !ms is the spin projection, tau is used to index the valleys
    real (dp),parameter :: PI=3.41592653_dp
end module
!
!---------------------------------------------------------------
program exciton
!---------------------------------------------------------------
use materials
implicit none
integer mesh, mout
integer n, m, i
integer,parameter :: mm = 5
  
real (dp) ::  zeta, zmesh, rmax, xmin, dx, e
real (dp), allocatable :: r(:), sqr(:), r2(:), y(:), vpot(:)
real (dp) :: En(mm)
  
character(len = 2) :: cTemp

    ! initialize atomic charge (Z)
    !
    zeta=Z/2.0_dp
    ! initialize logarithmic mesh
    !
    zmesh=  zeta 
    rmax = 100.0_dp
    xmin = -8.0_dp
    dx   =  0.001_dp
    !
    mesh = (log(zmesh*rmax) - xmin)/dx
    allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh) )
    !
    call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
    !
    ! initialize the potential
    !
    call init_pot ( zeta, r, mesh, vpot)
    !
    !
    ! read number of nodes (stop if nodes < 0)
    !
    do n=1,5
    write(cTemp,'(i2)') n !label the file name to distinguish different m   
    open(unit=12, file='En'//trim(adjustl(cTemp))//'.dat') !output a list of En(m)
    open(unit=11, file='wfc'//trim(adjustl(cTemp))//'.dat',status='unknown',form='formatted')  !output wavefunction
    
        do m=0,mm,1
        !
        ! solve the schroedinger equation in radial coordinates by Numerov method
        !
        call solve_sheq (n, m, mout, e, mesh, dx, r, sqr, r2, vpot, zeta,y )
        if (mout==1) exit
        !
        ! write out the eigenvalue energy to be compared with the external potential
        !
        write (12,'(" n=", I5, ", m=", I5, ", eigenvalue =",f17.8)') &
            n,m,En0*e
        write(11,"('#       r             R(r)          e')")
            do i=0,mesh
             write (11,*) r(i),y(i)*sqr(i), e
            end do
        write (11,'(/)')
  
        end do
    end do
close(12)
close(11)
pause
end program exciton
!
!--------------------------------------------------------------------
subroutine do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
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
  write(*,'(/" radial grid information:")')
  write(*,'(" dx =",f10.6,", xmin =",f10.6,", zmesh =",f10.6)') &
       dx,xmin,zmesh
  write(*,'(" mesh =",i6,", r(0) =",f10.6,", r(mesh) =",f10.6)') &
       mesh, r(0), r(mesh)
  write(*,*) 
  !
  return
end subroutine do_mesh
!---------------------------------------------------------------------
subroutine solve_sheq ( n, m, mout, e, mesh, dx, r, sqr, r2, vpot, zeta, y )
  !---------------------------------------------------------------------
  !
  ! solve the schroedinger equation in radial coordinates on a 
  ! logarithmic grid by Numerov method - atomic (Ry) units
  !
  use materials
  implicit none
  integer, parameter :: maxiter=1000
  real(dp), parameter :: eps=1.0D-10
  integer :: n
  integer, intent(in) :: mesh, m
  !!!The intent(in) attribute of argument i means that i cannot be changed inside the function.
  real(dp), intent(in) :: dx, r(0:mesh), sqr(0:mesh), r2(0:mesh), & 
       vpot(0:mesh), zeta
  real(dp), intent(out) :: e, y(0:mesh) !!!output
  integer, intent(out) :: mout
  
  integer :: i, j, icl, nodes, ncross, kkk
  real (dp) :: ddx12, sqlhf, x2l2, ycusp, dfcusp, fac, norm, eup, elw, de
  real (dp), allocatable :: f(:)
  real (dp) :: Amesh,temp,temp1,temp2
  integer :: a1,u
  
  integer, external :: factorial !()!
  
  allocate ( f(0:mesh) )
  
  ddx12=dx*dx/12.0_dp
  sqlhf = dble(m)**2.0_dp
  x2l2 = 2*abs(m)+1.0_dp
  
  ! set (very rough) initial lower and upper bounds to the eigenvalue
  !eup = vpot(mesh)
  !elw = minval (sqlhf/r2(:)+vpot(:))
  eup = vpot(mesh)
  elw = minval (sqlhf/r2(:)+vpot(:))
  mout = 0
  
  if (eup-elw < eps) then
     write (*,*) 'solve_sheq: lower and upper bounds are equal'
     mout = 1
     goto 100
  end if
  e = 0.5_dp * (elw + eup)

  do kkk = 1, maxiter
     ! set up the f-function and determine the position of its last
     ! change of sign
     ! f < 0 (approximately) means classically allowed region
     ! f > 0 forbidden
     icl = -1
     f(0) = ddx12 *( sqlhf + r2(0) * (vpot(0)-e) )
     
     do i=1,mesh
        f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
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
     !y(0) = r(0)**(x2l2/2) *(1.0_dp - 2.0_dp*zeta*r(0)/x2l2) / sqr(0)
     !y(1) = r(1)**(x2l2/2) *(1.0_dp - 2.0_dp*zeta*r(1)/x2l2) / sqr(1)
     y(0) = r(0)**(x2l2/2) *(1.0_dp + vpot(0)*r2(0)/x2l2) / sqr(0)
     y(1) = r(1)**(x2l2/2) *(1.0_dp + vpot(1)*r2(1)/x2l2) / sqr(1)
     !y(0) = r(0)**(x2l2/2) / sqr(0)
     !y(1) = r(1)**(x2l2/2) / sqr(1)
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
     !Amesh = 0
     !temp1 = ( dble(factorial(n)) * dble(factorial(n+m)) )**0.5
     !ul: do u = 0, n
     !   temp2 = dble(factorial(n-u)) * dble(factorial(m+u)) * dble(factorial(u))
     !   temp = (-1)**u * temp1 / temp2 * (beta/2.0D0)**(dble(m+1)/2.0D0+dble(u)) * r(mesh)**(dble(m)-0.5D0+dble(2*u))
     !   Amesh = Amesh + temp
     !end do ul
     
     y(mesh) = dx !+ PI**(-0.5) * exp( -beta /4.0D0 * r(mesh)**2) * Amesh
     y(mesh-1) = (12.0_dp-10.0_dp*f(mesh))*y(mesh)/f(mesh-1)
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
     !
  end do
  ! was convergence achieved ?
  !
  !if ( abs(de) > 1.d-10 ) then
  !   if ( ncross /= nodes ) then
  !      write (*,*) e, elw, eup, ncross, nodes, icl
  !   else
  !      write (*,*) e, de
  !   end if
  !   write (*,*) ' error in solve_sheq: too many iterations'
  !   stop
  !else 
  ! ---- convergence has been achieved -----
     write (*,'(" convergence achieved at iter #",i3," de = ", e11.4)') kkk,de
  !end if
100  return
end subroutine solve_sheq
!
!--------------------------------------------------------------------
subroutine init_pot ( zeta, r, mesh, vpot)
!--------------------------------------------------------------------
  !
  ! initialize potential
  !
  use materials
  implicit none
  integer, intent (in) :: mesh
  real (dp), intent(in) :: zeta, r(0:mesh)
  real (dp), intent(out):: vpot(0:mesh)
  real (dp), parameter :: gamma=-0.115931515659945_dp
  integer :: i

  open (7,file='pot.dat',status='unknown',form='formatted')
  do i =0,mesh
     vpot (i) = V0 * ( log( r(i)/(r(i)+ratio) ) + gamma * exp( -r(i)/ratio ) ) 
     write (7,*) r(i),vpot(i)
  end do
  close(7)
  return
end subroutine init_pot
    
function factorial(n) !()!
    integer :: n,i
    integer :: factorial
    factorial=1
    do i=1,n
        factorial=factorial*i
    end do
    return
end function
