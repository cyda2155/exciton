program main
!this part is used to calculate the Hamitonian matrix H spanned by phi_s
    USE mkl95_LAPACK    
    USE mkl95_PRECISION, ONLY: WP=>SP
    
    implicit none
    
    real*8,parameter :: Z=0.234741784037559,B0=1.235631801018422e+05
    real*8,parameter :: E0=13.950778695908895
    
    real*8,parameter :: PI=3.141593
    
    real*8 :: beta,Bb,dBb,Bb0,Bb1
    integer*8 :: m !m is the projection of angular momentum, we choose 0,1,2,... and the results of -m can be derived from m 
    integer*8 :: tau,ms !ms is the spin projection, tau is used to index the valleys
    
    integer*8 :: N,N1,Nt,b,c,k,s,j,u,t !N is the dimension of H, or the number of basis, N1 is the total number of Bb, b is the index of Bb, t is the index of all the eigen energy
    integer*8 :: Nout=3 !number of output eigen states
    real*8 :: Gks,temp,temp1,temp2,Htemp
    integer*8 :: max,a1 !a1 used to index temp
    
    integer*8,external :: accumulation !()£¡
    real*8,external :: accumulation1 !()£¡!
    
    real*8,allocatable :: H(:,:) !H relates to m whilst H1 relates to -m 
    real*8,allocatable :: En(:,:)
    real*8,allocatable :: Aks(:,:)
   
    character(len = 2) :: cTemp
 

    N=16
    tau=1 !+1 corresponds to K valley while -1 to -K valley
    ms=1 !+1 relates to spin up electrons while -1 to spin down(Sz)
    
    
    !magnetic field
    Bb0=1D0
    Bb1=100D0
    dBb=0.1D0 !precision
    N1=(Bb1-Bb0)/dBb+1
    
    allocate(H(N+1,N+1))
    allocate(En(N1,N+1))
    allocate(Aks(N+1,N+1))
    
    !scan magnetic field and projection of angular momentum
    do m=0,2
        
        En=0
        !firstly figure out Aks, deal with denominator and numerator simultanuously
        Aks=0   
        sl: do s=0,N
            kl: do k=0,N
                max=m+s+k
                jl: do j=0,k
                    ul: do u=j,s+j
                            temp=(PI)**0.5
                            al: do a1=1,max
                                temp1=(dble(accumulation(a1,k))*dble(accumulation(a1,s))*dble(accumulation(a1,k+m))*dble(accumulation(a1,s+m)))**0.5
                                temp2=dble(accumulation(a1,k-j))*dble(accumulation(a1,s-u+j))*dble(accumulation(a1,m+j))*dble(accumulation(a1,m+u-j))*dble(accumulation(a1,j))*dble(accumulation(a1,u-j))
                                temp=temp*temp1/temp2*accumulation1(dble(a1),dble(m+u))
                            end do al
                        temp=temp*(-1)**u
                        Aks(k+1,s+1)=Aks(k+1,s+1)+temp
                    end do ul
                end do jl
            end do kl
        end do sl
    
        Bb=Bb0
        do b=1,N1
            beta=Bb/B0
            !write(*,*)b !show calculating point
    
            H=0

             !calculate each element of H
            do s=0,N
                do k=0,N
                    if (k==s) then
                        H(k+1,s+1)=(beta**2+1.0D0)/2.0D0*(dble(2*s+m+1))-Z*(1.0D0/2.0D0)**0.5*Aks(k+1,s+1)
                    else
                        if (s==k+1 .OR. s==k-1) then
                            Gks=1.0D0
                            do a1=1,s+k+m
                                Gks=Gks*(dble(accumulation(a1,s))/dble(accumulation(a1,k))*dble(accumulation(a1,s+m))/dble(accumulation(a1,k+m)))**0.5D0
                            end do
                            H(k+1,s+1)=-(beta**2-1.0D0)/2.0D0*Gks-Z*(1.0D0/2.0D0)**0.5*Aks(k+1,s+1)
                        else
                            H=-Z*(1.0D0/2.0D0)**0.5*Aks(k+1,s+1)
                        end if
                    end if
                   !write(*,*)H(k+1,s+1) !check matrix elements
                end do
            end do
            
            !diagonalization
            call syev(H,En(b,:),'V')
    
            !scan B
        Bb=Bb+dBb
        end do !b loop is closed, now turn to another m
    
        !output eigen energy 
        write(cTemp,'(i2)') m !label the file name to distinguish different m   
        open(unit=10, file='En'//trim(adjustl(cTemp))//'.dat') !output a list of En(m)
            do c=1,Nout
                Bb=Bb0
                do b=1,N1
                     beta=Bb/B0
                     write(10,*)Bb,E0*(En(b,c)-beta*(2*(c-1)+m+1))!,E0*En(b,c)
                     !write(10,*)Bb,E0*En(b,c)
                     Bb=Bb+dBb
                end do
                write(10,*)
            end do
        !output another data file
        close(10)
    end do
    
    deallocate(H)
    deallocate(En)
    deallocate(Aks)
    !check status
    !write(*,*)"N=",N
end program main

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
    
