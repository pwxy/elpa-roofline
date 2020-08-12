! Brandon Cook's ELPA driver that was published in article Cook, Kurth, Deslippe, Carrier, Hill WichMann "Eigensolver performance comparison on Cray XC systems," Concurrency Computat Pract Exper. 2019, 31, https://doi.org/10.1002/cpe.4997
! Updated by Paul Lin July 2020 for ELPA 2020.05.001

module helper
  
   implicit none
  
  contains

    subroutine transpose_matrix(uplo, ln, a_ref, a_sym, desca)
      
      
      character, intent(in)  :: uplo
      integer, intent(in) :: ln
      complex(kind=8), intent(in) :: a_ref(:,:)
      complex(kind=8), intent(out) :: a_sym(:,:)
      integer, intent(in) :: desca(:)
      
      complex(kind=8) :: alpha, beta
      
      logical :: lsame
      
      alpha = (1.0, 0.0)
      beta = (0.0, 0.0)
      
      call pztranc(ln, ln, alpha, a_ref, 1, 1, desca, beta, a_sym, 1, 1, desca)
      
      if (lsame(uplo, "L")) then
         call pztrmr2d("L", "N", ln, ln, a_ref, 1, 1, desca, a_sym, 1, 1, desca, desca(2))
      else 
         call pztrmr2d("U", "N", ln, ln, a_ref, 1, 1, desca, a_sym, 1, 1, desca, desca(2))
      end if
      

    end subroutine transpose_matrix

    elemental subroutine convert_to_int(str, i, stat)
      
      character(len=*), intent(in) :: str
      integer, intent(out) :: i, stat

      read(str, *, iostat=stat) i
      
    end subroutine convert_to_int

end module helper

program benchmark_elpa

  use elpa
!  use elpa_driver   ! not in ELPA 2020.05.001
  use helper

  implicit none

  include 'mpif.h'
  
  integer :: ctxt_sys, my_blacs_ctxt
  integer :: rank, size, i, j
  integer :: ln, lnprow, lnpcol, lnbrow, lnbcol
  integer :: llda, ml, nl, lsize, numEle
  integer :: myrow, mycol, info, err
  real(kind=8) :: t1, t2

  class(elpa_t), pointer :: e

!  integer, allocatable :: desca(:)
  integer :: desca(9)
  complex(kind=8), allocatable :: a_ref(:,:), a_sym(:,:), z(:,:)
  real(kind=8), allocatable :: w(:)
  integer, dimension(4) :: iseed
  character(len=10), dimension(4) :: argc

  integer numroc, mpi_comm_rows, mpi_comm_cols

  integer, parameter :: ELS_TO_PRINT  = 6             ! arbitrary els to print out

  integer(kind=c_int)          :: error_elpa


  ! Initialize MPI and BLACS
  call mpi_init(err)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, err)
  call mpi_comm_size(MPI_COMM_WORLD, size, err)

  do i = 1, 4
     call get_command_argument(i, argc(i))
  end do

  call convert_to_int(argc(1), ln, err) 
  call convert_to_int(argc(2), lnprow, err)
  call convert_to_int(argc(3), lnpcol, err)
  call convert_to_int(argc(4), lnbrow, err)
  lnbcol = lnbrow

  if (rank == 0) then
     print *, "ln =", ln 
     print *, "lnprow =", lnprow
     print *, "lnpcol =", lnpcol
     print *, "lnbropw =", lnbrow
  end if

!  call blacs_get(0, 0, ctxt_sys)
!  call blacs_gridinit(ctxt_sys, "C", lnprow, lnpcol)
!  call blacs_gridinit(ctxt_sys, "R", lnprow, lnpcol)
!  call blacs_gridinfo(ctxt_sys, lnprow, lnpcol, myrow, mycol)

  my_blacs_ctxt = MPI_COMM_WORLD
  call blacs_gridinit(my_blacs_ctxt, 'R', lnprow, lnpcol)
  call blacs_gridinfo(my_blacs_ctxt, lnprow, lnpcol, myrow, mycol)


  if (myrow .eq. -1) then
     print *, "Failed to properly initialize MPI and/or BLACS!"
     call MPI_FINALIZE(err)
     stop
  end if

  ! Explicitly get and set the row and column communicators, as the API seems to be
  ! failing to initialize them as they should
  ! NOTE: elpa_get_communicators not in ELPA 2020.05.001
!  err = elpa_get_communicators(MPI_COMM_WORLD, myrow, mycol, mpi_comm_rows, mpi_comm_cols)


  ! Allocate my matrices now
  ml = numroc(ln, lnbrow, myrow, 0, lnprow)
  nl = numroc(ln, lnbcol, mycol, 0, lnpcol)
  llda = ml

  allocate(a_ref(ml,nl))
  allocate(a_sym(ml,nl))
  allocate(z(ml,nl))
  allocate(w(ln))
!  allocate(desca(9))

  ! Create my blacs descriptor for transposing the matrix
!    call descinit(desca, ln, ln, lnbrow, lnbcol, 0, 0, ctxt_sys, llda, info)  
  call descinit(desca, ln, ln, lnbrow, lnbcol, 0, 0, my_blacs_ctxt, llda, info)

  iseed(1) = myrow
  iseed(2) = mycol
  iseed(3) = mycol + myrow*lnpcol
  iseed(4) = 1
  if (iand(iseed(4), 2) == 0) then
     iseed(4) = iseed(4) + 1
  end if

  ! Try initializing and allocating elpa
  if (elpa_init(20170403) /= elpa_ok) then
     print *, "ELPA API not supported"
     stop
  end if

  e => elpa_allocate(err)
  if (err /= ELPA_OK) then
    print *,"crap!  elpa_allocate failed"
  endif

  numEle = ml*nl
!  write(6,*) iseed, numEle

  call zlarnv(1, iseed, numEle, a_ref)
  call transpose_matrix("L", ln, a_ref, a_sym, desca) 

! end do
! do i=1, ml
!    do j=1, ml
!       write(6,*) "rank=", rank, "a_ref(",i,j,")=", a_ref(i,j)
!    end do
! end do
  deallocate(a_ref) ! We don't need it anymore, used only for constructing a_sym

t1 = MPI_WTIME()

  call e%set("na", ln, err)
  call e%set("nev", ln, err)
  call e%set("nblk", lnbrow, err)
  call e%set("local_nrows", ml, err)
  call e%set("local_ncols", nl, err)
  call e%set("process_row", myrow, err)
  call e%set("process_col", mycol, err)
  call e%set("mpi_comm_parent", MPI_COMM_WORLD, err)

   call e%set("solver", ELPA_SOLVER_2STAGE, err)
!   call e%set("gpu", 0, err)   ! no gpu offload
   call e%set("gpu", 1, err)
!  call e%set("complex_kernel", ELPA_2STAGE_COMPLEX_AVX2_BLOCK1, err)
!  call e%set("complex_kernel", ELPA_2STAGE_COMPLEX_AVX512_BLOCK1, err)
   call e%set("complex_kernel", ELPA_2STAGE_COMPLEX_GPU, err)

   err=e%setup()

  call e%eigenvectors(a_sym, w, z, err)
  call elpa_deallocate(e)
  call elpa_uninit()
  
  t2 = MPI_WTIME()
  
  if (rank == 0) then
     print *, "ELPA time: ", t2 - t1
  end if
  
  ! deallocate stuff

   if (rank == 0 ) then
      do j=1, ELS_TO_PRINT
         write(6,*) "w(", j, ")=", w(j), "z(", j, ")=", z(j,1)
      end do
   end if

  deallocate(a_sym)
  deallocate(z)
  deallocate(w)
!  deallocate(desca)

  ! Finish MPI
!  call MPI_FINALIZE(err)

  call blacs_gridexit(my_blacs_ctxt)
  call blacs_exit(0)

  ! Finish MPI
!  call MPI_FINALIZE(err)

end program benchmark_elpa
