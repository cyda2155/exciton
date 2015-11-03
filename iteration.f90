program main

    USE mkl95_LAPACK    
    USE mkl95_PRECISION, ONLY: WP=>SP
    
    implicit none
    
    !equation defined variables
    real*8,parameter :: Z=2.347417840375587e+02,B0=0.123563180101842D0
    real*8,parameter :: E0=1.395077869590889e-05
    integer*8 :: tau=1 !+1 corresponds to K valley while -1 to -K valley
    integer*8 :: m !ANGULAR MOMENTUM
    
    real*8,parameter :: PI=3.141593D0
    real*8 :: A00
    real*8,allocatable :: A(:,:)
    
    !user defined variables

    integer*8 :: N=16 !DIMENSION OF THE MATRIX
    integer*8 :: k,s
    integer*8 :: N1 !ACCURACY OF MAGNETIC FIELD
    integer*8 :: Nout=3 !number of output eigen states
    real*8 :: beta,Bb,dBb,Bb0,Bb1
    
    !program defined variables
    real*8 :: Atemp
    real*8 :: Htemp
    integer*8 :: i
    integer*8,external :: accumulation !()£¡
    real*8,external :: accumulation1 !()£¡!
    real*8,allocatable :: H(:,:) !H relates to m1 whilst H1 relates to -m1 
    real*8,allocatable :: En(:,:)
    integer*8 :: c,b
    character(len = 2) :: cTemp
    
    allocate(A((N+1),(N+1)))
    
    !magnetic field
    Bb0=1.0D0
    Bb1=100.0D0
    dBb=0.1D0 !precision
    N1=(Bb1-Bb0)/dBb+1
    allocate(H(N+1,N+1))
    allocate(En(N1,N+1))
    
    !LOOP OF m
    do m=0,2
        !clear En
        En=0
        En1=0 
        
        A=0
        !!Aks
        !!INITIALIZATION
        A00=sqrt(PI)
    
        Atemp=1
        if (m/=0) then
            do i=1,m
                Atemp=Atemp*accumulation1(dble(i),dble(m))/dble(accumulation(i,m))
            end do
        end if
        A(1,1)=A00*Atemp
    
        do k=1,N
            s=k
            A(k+1,1)=1.0D0/sqrt(dble((k+m)*k))*(dble(k)-0.5D0)*A(k,1)
            A(1,s+1)=A(k+1,1)
        end do
    
        do k=1,N
            s=k
            A(k+1,2)=1.0D0/sqrt(dble(1+m))*((0.5D0-dble(k))*A(k+1,1)+sqrt(dble(k*(k+m)))*A(k,1))
            A(2,s+1)=A(k+1,2)
        end do
        !!ITERATION
        do s=2,N
            do k=s,N
                A(k+1,s+1)=(dble(m+k+s-(k-s)**2)-5.0D0/4.0D0)/sqrt(dble(k*(k+m)*s*(s+m)))*A(k,s) &
                +sqrt(dble((k-1+m)*(s-1+m)*(k-1)*(s-1))/dble(k*(k+m)*s*(s+m)))*A(k-1,s-1) &
                +sqrt(dble((k-1+m)*(k-1))/dble(k*(k+m)))*((dble(k-s)-0.5D0)/sqrt(dble(s*(s+m))))*A(k-1,s) &
                +sqrt(dble((s-1+m)*(s-1))/dble(s*(s+m)))*((dble(s-k)-0.5D0)/sqrt(dble(k*(k+m))))*A(k,s-1)
                A(s+1,k+1)=A(k+1,s+1)
            end do
        end do
        
        !!LOOP of b
        Bb=Bb0
        do b=1,N1
            !clear H
            H=0
            H1=0
            beta=Bb/B0
            !open(unit=9, file='Entemp.dat') !output a list of En(m)
            !!HAMILTONIAN Hks
            do s=0,N
                do k=0,N
                    if (k==s) then
                        H(k+1,s+1)=beta*(dble(2*s+m+1))-Z*(beta/2)**0.5*A(k+1,s+1)
                    else
                        H(k+1,s+1)=-Z*(beta/2)**0.5*A(k+1,s+1)
                    end if
                    !if (tau<0) then
                    !    Htemp=H(k+1,s+1)
                    !    H(k+1,s+1)=H1(k+1,s+1)
                    !    H1(k+1,s+1)=Htemp
                    !end if
                end do
            end do
            !!EIGEN ENERGY

            !diagonalization
            call syev(H,En(b,:),'V')
            
            !write(9,*)En(b,:)
        !scan B
        Bb=Bb+dBb
        end do !b loop is closed, now turn to another m1
        !!!OUTPUT
        !open(unit=10, file='Aks.txt')
        !do k=1,N+1
        !    write(10,*)A(k,:)
        !end do
        
        !output eigen energy 
        write(cTemp,'(i2)') m !label the file name to distinguish different m   
        open(unit=10, file='En'//trim(adjustl(cTemp))//'.dat') !output a list of En(m)
        do c=1,Nout
            Bb=Bb0
            do b=1,N1
                 beta=Bb/B0
                 write(10,*)Bb,E0*(En(b,c)-beta*(2*(c-1)+m+1))
                 Bb=Bb+dBb
            end do
            write(10,*)
        end do
        !output another data file
        close(10)
    end do !m loop is closed
    
    deallocate(H)
    deallocate(En)
    !check status
    write(*,*)"succeed"
end program
    
function accumulation(j,n) !i+ until reach the upper limit
    integer*8 :: j,n
    integer*8 :: accumulation
    if (j<=n) then
        accumulation=j
    else
        accumulation=1
    end if
    return
end function

function accumulation1(j,n) !i+ until reach the upper limit
    real*8 :: j,n
    real*8 :: accumulation1
	if (j<=n) then
        accumulation1=j-0.5D0
    else
        accumulation1=1
    end if
    return
end function
