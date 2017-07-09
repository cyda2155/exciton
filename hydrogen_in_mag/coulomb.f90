!!exciton energy for MoS2
    
module materials
    !material determined parameters
    real*8,parameter :: Z=2.347417840375587e+02,B0=0.123563180101842
    real*8,parameter :: En0=1.395077869590890e-05
    integer*8 :: tau=1 !+1 corresponds to K valley while -1 to -K valley
end module   
    
    
program main
    
    use materials
    
    implicit none
    
    !equation defined variables
    integer*8 :: m !ANGULAR MOMENTUM
    
    real*8,parameter :: PI=3.141593D0
    
    !user defined variables
    integer*8 :: s
    integer*8 :: N=5 !
    real*8,allocatable :: E(:,:)
    
    integer*8 :: N1 !ACCURACY OF MAGNETIC FIELD
    real*8 :: beta,Bb,dBb,Bb0,Bb1
    
    !program defined variables
    real*8,allocatable :: A(:),E0(:)
    real*8 :: temp
    integer*8 :: i
    integer*8 :: b
    character(len = 2) :: cTemp
    
    
    !magnetic field
    Bb0=1.0D0
    Bb1=100.0D0
    dBb=0.1D0 !precision
    N1=(Bb1-Bb0)/dBb+1
    beta=0
    
    allocate(E(N1,N+1))
    allocate(A(N+1),E0(N+1))
    E=0
    
    !LOOP OF m
    do m=0,2
        A=0
        E0=0
        !!Ann
        do s=0,N
            A(s+1)=dble((1+2*m+2*s)**2*(3+5*s*(1+s)+5*m+10*s*m+2*m**2))/(8.0D0*Z**2)
            !write(*,*)A(s)
            E0(s+1)=-Z**2/dble((1+2*m+2*s)**2)
        end do
        
        !!LOOP of b
        Bb=Bb0
        do b=1,N1
            beta=Bb/B0
            do s=0,N
                temp=dble(m)*beta+beta**2.0*A(s+1)
                !!consistency of the 1st order perturbation
                if (temp>0.1D0*abs(E0(s+1)) .AND. b>0.1*N1) then
                E(b,s+1)=0
                else
                E(b,s+1)=E0(s+1)+temp
                end if
                !E(b,s+1)=E0(s+1)+temp
            end do
        !scan B
        Bb=Bb+dBb
        end do !b loop is closed, now turn to another m1

        !output eigen energy 
        write(cTemp,'(i2)') m !label the file name to distinguish different m   
        open(unit=10, file='En'//trim(adjustl(cTemp))//'.dat') !output a list of En(m)
        do s=0,N
            Bb=Bb0
            do b=1,N1
                 beta=Bb/B0
                 if (E(b,s+1)==0) then
                 cycle
                 else
                 write(10,*)Bb,En0*(E(b,s+1)-beta*(2*s+2*m+1))
                 end if
                 !write(*,*)beta,E(b,s),beta*(2*s+2*m+1)
                 Bb=Bb+dBb
            end do
            write(10,*)
        end do
        !output another data file
        close(10)
        
    end do !m loop is closed
    
    !check status
    write(*,*)"succeed"
    PAUSE
end program