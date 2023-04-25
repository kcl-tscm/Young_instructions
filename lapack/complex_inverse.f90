program complex_matrix_inverse
  implicit none
  complex*16,dimension(:,:),allocatable ::  A,A1
  complex*16,dimension(:),allocatable ::  WORK
  integer,dimension(:),allocatable :: IPIV
  integer :: INFO,LWORK
  integer :: i,j,k,err,N
  complex*16 :: c,error
  
!.......... read the matrix and dimension
  
      open(1,file='test.dat',form='formatted',status='old')
      read(1,*) N 
      
      allocate(A(N,N),A1(N,N),WORK(N),IPIV(N),stat=err)
      if (err.ne.0)then
         print *, "error: not enough memory"
         stop
      end if

      do i=1,N
         read(1,*) (A(i,j),j=1,N)
         write(*,'(i2,x,10(a,f10.5,x,f10.5,a,x))') i, ('(',A(i,j),')',j=1,N)
      end do
      close (1)
      A1=A
      
!......... run einverse matrix routine

!_______ determine first the dimension of WORK

      call ZGETRF(N,N,A,N,IPIV,info)
      
      if(info .eq. 0) then
         write(*,*)"succeded"
      else
         write(*,*)"failed"
      end if
      
 
      call ZGETRI(N,A,N,IPIV,WORK,N,INFO)
      if (err.ne.0)then
         print *,"error:fail to release"
         stop
      end if
      
!........ check eigenvectors & eigenvalues

      write(*,'(/a,i3)') '....... CHECKING MATRIX via A1*A-E=0 .......'

      error=(0.0d0,0.0d0)
      do i=1,N
         do j=1,N
            c=(0.0d0,0.0d0)
            do k=1,N
               c=c+A1(i,k)*A(k,j)
            end do
            if(i.eq.j) c=c-1.0d0
            err=err+c
            write(*,*) i,j,c
         end do
         write(*,'(a,2(e12.6,x))')' error= ',error
      end do

      deallocate(A,A1,IPIV,WORK,stat=error)
      
end program complex_matrix_inverse
