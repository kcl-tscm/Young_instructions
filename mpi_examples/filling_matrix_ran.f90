program matrix_fill
  use mpi
  implicit none

  integer :: my_rank, comm_size, ierr, n, m, i, j, rows_per_proc
  real, allocatable :: local_matrix(:,:), matrix(:,:)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr)

  if (my_rank == 0) then
    print *, "Enter the dimensions of the matrix (n x m): "
    read *, n, m
    allocate(matrix(n, m)) ! Allocate matrix only on root process
  end if

  call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! Divide rows evenly among processes
  rows_per_proc = n / comm_size
  allocate(local_matrix(rows_per_proc, m))

  ! Scatter matrix blocks to processes, ensuring all rows are scattered
  call MPI_Scatter(matrix, rows_per_proc * m, MPI_REAL, local_matrix, &
                   rows_per_proc * m, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

  ! Fill local matrix blocks
  call fill_local_array(local_matrix, rows_per_proc, m, my_rank)

  ! Gather filled blocks back to root process
  call MPI_Gather(local_matrix, rows_per_proc * m, MPI_REAL, matrix, &
                   rows_per_proc * m, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

  ! If root process, print the filled matrix
  if (my_rank == 0) then
    print *, "Filled matrix:"
    call print_matrix(matrix, n, m)
  end if

  deallocate(local_matrix)
  if (my_rank == 0) deallocate(matrix)

  call MPI_Finalize(ierr)

contains

  subroutine fill_local_array(local_matrix, rows, cols, my_rank)
    real, intent(inout) :: local_matrix(rows, cols)
    integer, intent(in) :: rows, cols, my_rank
    integer :: i, j
    integer :: seed_size

    seed_size = cols
    call random_seed(size = seed_size)
    do i = 1, rows
      call random_number(local_matrix(i,:))
      local_matrix(i,:) = local_matrix(i,:) + my_rank
    end do
  end subroutine fill_local_array

  subroutine print_matrix(matrix, n, m)
    integer :: i, j, n, m
    real, intent(in) :: matrix(n,m)

    do i = 1, n
      write(*,*) (matrix(i, j), j = 1, m)
      write(*,*)
    end do
 end subroutine print_matrix

end program matrix_fill
