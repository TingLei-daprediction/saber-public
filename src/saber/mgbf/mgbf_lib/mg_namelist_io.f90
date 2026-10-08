module mg_namelist_io
!$$$  module documentation block
!                .      .    .                                       .
! module:   mg_namelist_io
!   prgmmr: lei              org:                     date: 2026-10-07
!
! abstract:  Read a namelist file on one rank and broadcast its text, so
!            that thousands of ranks do not open the same file at once.
!            Each rank then parses the namelist from the returned lines
!            with an internal-file read (Fortran 2003):
!
!              call read_namelist_lines(filename, comm, lines)
!              read(lines, nml=my_group)
!
! module history log:
!   2026-10-07  lei     - initial version
!
! Subroutines Included:
!   read_namelist_lines -
!
! attributes:
!   language: f2003
!
!$$$ end documentation block

  use mpi
  implicit none

  private

  public :: read_namelist_lines

  integer, parameter, public :: nml_line_len = 512

contains

subroutine read_namelist_lines(filename, comm, lines)
!***********************************************************************
!                                                                      !
! Rank 0 of comm reads filename line by line; all ranks of comm        !
! receive the same lines. Aborts if the file cannot be opened or a     !
! line is longer than nml_line_len, which would otherwise be           !
! truncated silently.                                                  !
!                                                                      !
!***********************************************************************
character(len=*), intent(in) :: filename
integer, intent(in) :: comm
character(len=nml_line_len), allocatable, intent(out) :: lines(:)

! One extra character so that an over-long line can be detected
character(len=nml_line_len+1) :: buffer
integer :: mype, nlines, iline, myunit, ios, ierr
!-----------------------------------------------------------------------

call MPI_Comm_rank(comm, mype, ierr)

nlines = 0
if (mype == 0) then
  open(newunit=myunit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(6,*) 'read_namelist_lines: cannot open ', trim(filename)
    call flush(6)
    call MPI_Abort(comm, 1, ierr)
  end if

  ! First pass: count lines and check their length
  do
    read(myunit, '(A)', iostat=ios) buffer
    if (ios /= 0) exit
    if (len_trim(buffer) > nml_line_len) then
      write(6,*) 'read_namelist_lines: line ', nlines+1, ' of ', trim(filename), &
                 ' is longer than ', nml_line_len, ' characters'
      call flush(6)
      call MPI_Abort(comm, 1, ierr)
    end if
    nlines = nlines + 1
  end do

  ! Second pass: store the lines
  allocate(lines(nlines))
  rewind(myunit)
  do iline = 1, nlines
    read(myunit, '(A)') lines(iline)
  end do
  close(myunit)
end if

call MPI_Bcast(nlines, 1, MPI_INTEGER, 0, comm, ierr)
if (mype /= 0) allocate(lines(nlines))
if (nlines > 0) then
  call MPI_Bcast(lines, nml_line_len*nlines, MPI_CHARACTER, 0, comm, ierr)
end if

end subroutine read_namelist_lines

end module mg_namelist_io
