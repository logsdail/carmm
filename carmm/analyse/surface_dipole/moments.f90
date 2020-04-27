! AJL 02/2013
!
! I'm putting all the functions in here to make my life "easier"
!
  subroutine print_second_moment(nlength, nwidth, xyzq, centre_of_mass, the_second_moment)
!
  integer :: i, j
  integer :: nlength, nwidth
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), the_second_moment(3,3)
!
  real    :: second_moment
!
  the_second_moment = 0.0d0
!
  write(*,*)
  write(*,*)'SECOND MOMENT'
  write(*,*)'==================='
!
! Loop over all directions
!
  do i=1,3
    do j=1,3
!
! Calculate Second Moment
!
      the_second_moment(i,j) = second_moment(nlength, nwidth, i, j, xyzq, centre_of_mass)
!
! Print Second Moment
!
      write(*,*) i, ',', j, ' = ', the_second_moment(i,j), ' e.Ang^2'
    enddo
  enddo

  return
  end

  subroutine print_third_moment(nlength, nwidth, xyzq, centre_of_mass, the_third_moment)
!
  integer :: i, j
  integer :: nlength, nwidth
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), the_third_moment(3,3,3)
!
  real    :: third_moment
!
  the_third_moment = 0.0d0
!
  write(*,*)
  write(*,*)'THIRD MOMENT'
  write(*,*)'==================='
!
! Loop over all directions
!
  do i=1,3
    do j=1,3
      do k=1,3
!
! Calculate Third Moment
!
        the_third_moment(i,j,k) = third_moment(nlength, nwidth, i, j, k, xyzq, centre_of_mass)
!
! Print Third Moment
!
        write(*,*) i, ',', j, ',', k, ' = ', the_third_moment(i,j,k), ' e.Ang^3'
      enddo
    enddo
  enddo

  return
  end

