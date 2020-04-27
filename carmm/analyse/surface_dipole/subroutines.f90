! AJL 02/2013
!
! I'm putting all the functions in here to make my life "easier"
!
  subroutine calculate_centre_of_mass(nmax, nlength, nwidth, xyzq, minimums, maximums, centre_of_mass) 
  integer :: i, j
  integer :: nmax, nlength, nwidth
  real    :: xyzq(nmax,nwidth), minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth)
!
! Clean arrays
!
  minimums = 10000.0d0
  maximums = -10000.0d0
  centre_of_mass = 0.0d0
!
! Loop over all values
!
  do i=1,nlength
    do j=1,nwidth
!
! Sum totals
!
      centre_of_mass(j) = centre_of_mass(j) + xyzq(i,j)
!
! Check minimums
!
      if (xyzq(i,j) < minimums(j)) minimums(j) = xyzq(i,j)
!
! Check maximums
!
      if (xyzq(i,j) > maximums(j)) maximums(j) = xyzq(i,j)
    enddo
  enddo

!
! Averages for centre_of_mass
!
  do j=1,nwidth
    centre_of_mass(j) = centre_of_mass(j) / nlength
  enddo
!
  return
  end

  subroutine print_coords(title, nmax, nlength, nwidth, xyzq)
  character(*) :: title
  integer        :: i, nmax, nlength, nwidth
  real           :: xyzq(nmax,nwidth)
!
  write(*,*)
  write(*,*)
  write(*,*) title
  write(*,*) '==================='
  do i=1,nlength
    write(*,*) xyzq(i,1), xyzq(i,2), xyzq(i,3), xyzq(i,4)
  enddo

  return
  end

  subroutine print_min_max_com(nwidth, minimums, maximums, centre_of_mass)
  integer :: nwidth
  real    :: minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth)
!
  write(*,*)'MINIMA'
  write(*,*)'==================='
  write(*,*) minimums(1), minimums(2), minimums(3), minimums(4)
  write(*,*)'MAXIMA'
  write(*,*)'==================='
  write(*,*) maximums(1), maximums(2), maximums(3), maximums(4)
  write(*,*)'CENTRE OF MASS'
  write(*,*)'==================='
  write(*,*) centre_of_mass(1), centre_of_mass(2), centre_of_mass(3), centre_of_mass(4)
!
  return
  end

  subroutine print_dipoles(nmax, nlength, nwidth, xyzq, centre_of_mass, the_dipole)
!
  integer :: i
  integer :: nmax, nlength, nwidth
  real    :: xyzq(nmax,nwidth), centre_of_mass(nwidth), the_dipole(3)
!
  real    :: dipole
!
  the_dipole = 0.0d0
!
  write(*,*)'DIPOLE'
  write(*,*)'==================='
!
! Loop over all directions
!
  do i=1,3
!
! Calcule dipole
!
    the_dipole(i) = dipole(nmax, nlength, nwidth, i, xyzq, centre_of_mass)
!
! Print Dipole
!
    write(*,*) i, ' = ', the_dipole(i)
  enddo

  return
  end

  subroutine print_second_moment(nmax, nlength, nwidth, xyzq, centre_of_mass, the_second_moment)
!
  integer :: i, j
  integer :: nmax, nlength, nwidth
  real    :: xyzq(nmax,nwidth), centre_of_mass(nwidth), the_second_moment(3,3)
!
  real    :: second_moment
!
  the_second_moment = 0.0d0
!
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
      the_second_moment(i,j) = second_moment(nmax, nlength, nwidth, i, j, xyzq, centre_of_mass)
!
! Print Second Moment
!
      write(*,*) i, ',', j, ' = ', the_second_moment(i,j)
    enddo
  enddo

  return
  end

  subroutine print_third_moment(nmax, nlength, nwidth, xyzq, centre_of_mass, the_third_moment)
!
  integer :: i, j
  integer :: nmax, nlength, nwidth
  real    :: xyzq(nmax,nwidth), centre_of_mass(nwidth), the_third_moment(3,3,3)
!
  real    :: third_moment
!
  the_third_moment = 0.0d0
!
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
        the_third_moment(i,j,k) = third_moment(nmax, nlength, nwidth, i, j, k, xyzq, centre_of_mass)
!
! Print Third Moment
!
        write(*,*) i, ',', j, ',', k, ' = ', the_third_moment(i,j,k)
      enddo
    enddo
  enddo

  return
  end

  subroutine print_quadrupole(nmax, nlength, nwidth, xyzq, centre_of_mass, the_quadrupole)
!
  integer :: i, j
  integer :: nmax, nlength, nwidth
  real    :: xyzq(nmax,nwidth), centre_of_mass(nwidth), the_quadrupole(3,3), trace
!
  real    :: quadrupole
!
  the_quadrupole = 0.0d0
!
  write(*,*)'QUADRUPOLE'
  write(*,*)'==================='
!
! Loop over all directions
!
  do i=1,3
    do j=1,3
!
! Calculate Quadrupole
!
      the_quadrupole(i,j) = quadrupole(nmax, nlength, nwidth, i, j, xyzq, centre_of_mass)
!
! Print Quadrupole
!
      write(*,*) i, ',', j, ' = ', the_quadrupole(i,j)
    enddo
  enddo

  trace = the_quadrupole(1,1) + the_quadrupole(2,2) + the_quadrupole(3,3)

  write(*,*)'TRACE = ', trace

  return
  end

  subroutine print_quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole, the_quadrupole_decomposed)
!
  integer :: nwidth
  real    :: the_quadrupole(3,3), minimums(nwidth), maximums(nwidth), the_quadrupole_decomposed
!
  the_quadrupole_decomposed = 0.0d0
!
  write(*,*)'QUADRUPOLE DECOMPOSED'
  write(*,*)'==================='
!
  the_quadrupole_decomposed = quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole) 

  write(*,*)'DIPOLAR CHARGE = ', the_quadrupole_decomposed

  return
  end

  subroutine print_voltage_offset(nwidth, minimums, maximums, centre_of_mass, unit_cell, the_quadrupole_decomposed, &
                                  the_voltage_offset)
!
  integer :: nwidth, i, j
  real    :: minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth), unit_cell(3,3)
  real    :: the_quadrupole_decomposed, the_voltage_offset
  real    :: unit_cell_sum
!
  the_voltage_offset = 0.0d0
!
! Check unit cell is defined
!
  unit_cell_sum = 0.0d0
  do i=1,3
    do j=1,3
      if (unit_cell(i,j) > 0.0d0) then
        unit_cell_sum = unit_cell_sum + unit_cell(i,j)
      endif
    enddo
  enddo
!
! If so then calculate voltage offset
!
  if (unit_cell_sum > 0.0d0) then
!
    write(*,*)'VOLTAGE OFFSET'
    write(*,*)'==================='
!
    the_voltage_offset = voltage_offset(nwidth, minimums, maximums, centre_of_mass, unit_cell, the_quadrupole_decomposed)

    write(*,*)'VOLTAGE OFFSET = ', the_voltage_offset
!
  endif

  return
  end

  subroutine print_diagonalised_quadrupole(A)
!
  real    :: A(3,3)
  real    :: W(3)
  integer :: N, LDA, LWMAX, INFO, LWORK
  integer :: i
  parameter  ( N = 3 )
  parameter  ( LDA = N )
  parameter  ( LMAX = 1000 )
  real    :: WORK(LMAX)
!
  write(*,*)'EIGENVALUES OF QUADRUPOLE'
  write(*,*)'==================='
!
  W = 0.0d0
  INFO = 0.0d0
  WORK = 0.0d0
  LWORK = 3 * N
!
  call dsyev('V', 'Upper', N, A, LDA, W, WORK, LWORK, INFO)
!
  write(*,*) ' VALUES = ', W(1), ',', W(2), ',', W(3)

  do i=1,3
    write(*,*) 'VECTORS = ', A(i,1), ' , ',  A(i,2), ' , ', A(i,3)
  enddo
!
  return
  end

  subroutine print_octupole(nmax, nlength, nwidth, xyzq, centre_of_mass, the_octupole)
!
  integer :: i, j
  integer :: nmax, nlength, nwidth
  real    :: xyzq(nmax,nwidth), centre_of_mass(nwidth), the_octupole(3,3,3), trace
!
  real    :: octupole
!
  the_octupole = 0.0d0
!
  write(*,*)'OCTUPOLE'
  write(*,*)'==================='
!
! Loop over all directions
!
  do i=1,3
    do j=1,3
      do k=1,3
!
! Calculate Octupole
!
        the_octupole(i,j,k) = octupole(nmax, nlength, nwidth, i, j, k, xyzq, centre_of_mass)
!
! Print Octupole
!
        write(*,*) i, ',', j, ',', k, ' = ', the_octupole(i,j,k)
      enddo
    enddo
  enddo

! As taken from "The Theory of Intermolecular Forces", Anthony Stone (Page 17)

  trace = the_octupole(1,1,1) + the_octupole(1,2,2) + the_octupole(1,3,3)

  write(*,*)'TRACE = ', trace

  return
  end

  subroutine add_new_q_values_to_coordinates(nmax, nlength, nwidth, xyzq, minimums, maximums, centre_of_mass, &
  the_dipole_neutralising_layer, the_dipole_neutralising_q, the_dipole_neutralising_density, unit_cell)
!
  integer :: i, j, q_squared
  integer :: nmax, nlength, nwidth, the_dipole_neutralising_density
  real    :: xyzq(nmax,nwidth), minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth)
  real    :: the_dipole_neutralising_layer, the_dipole_neutralising_q
  real    :: grid_steps(2)
  integer :: temp(2)
  real    :: unit_cell(3,3), unit_cell_sum
!  logical :: fractional
!
! Check unit cell is defined
!
  unit_cell_sum = 0.0d0
  do i=1,3
    do j=1,3
      if (unit_cell(i,j) > 0.0d0) then
        unit_cell_sum = unit_cell_sum + unit_cell(i,j)
      endif
    enddo
  enddo
!
! If so, redefine minimums and maximums
!
  if (unit_cell_sum > 0.0d0) then
    do i=1,3
      minimums(i) = 0.0d0
      maximums(i) = 1.0d0 
    enddo
  endif
!
! Calculate distance for each step
!
  do j=1,2
    grid_steps(j) = (maximums(j)-minimums(j))/the_dipole_neutralising_density
  enddo
!
! Add points above surface
!
  q_squared = the_dipole_neutralising_density**2
  temp = 0
  do i=1,q_squared
    if (the_dipole_neutralising_density > 1) then
!      if (unit_cell_sum > 0.0d0 .and. .not. fractional) then
      if (unit_cell_sum > 0.0d0) then
        xyzq(nlength+i,1) = minimums(1)+((temp(1)*grid_steps(1))*unit_cell(1,1))+((temp(2)*grid_steps(2))*unit_cell(2,1))
        xyzq(nlength+i,2) = minimums(2)+((temp(1)*grid_steps(1))*unit_cell(1,2))+((temp(2)*grid_steps(2))*unit_cell(2,2))
!      else
!        xyzq(nlength+i,1) = minimums(1)+temp(1)*grid_steps(1)
!        xyzq(nlength+i,2) = minimums(2)+temp(2)*grid_steps(2)
      endif
    else
      xyzq(nlength+i,1) = centre_of_mass(1)
      xyzq(nlength+i,2) = centre_of_mass(2)
    endif
    xyzq(nlength+i,3) = centre_of_mass(3)+the_dipole_neutralising_layer
    xyzq(nlength+i,4) = the_dipole_neutralising_q/q_squared

    temp(1) = temp(1) + 1
    if (temp(1) >= the_dipole_neutralising_density) then
      temp(1) = 0
      temp(2) = temp(2) + 1
    endif
  enddo

  nlength = nlength + q_squared
!
  temp = 0
  do i=1,q_squared
    if (the_dipole_neutralising_density > 1) then
!      if (unit_cell_sum > 0.0d0 .and. .not. fractional) then
      if (unit_cell_sum > 0.0d0) then
        xyzq(nlength+i,1) = minimums(1)+((temp(1)*grid_steps(1))*unit_cell(1,1))+((temp(2)*grid_steps(2))*unit_cell(2,1))
        xyzq(nlength+i,2) = minimums(2)+((temp(1)*grid_steps(1))*unit_cell(1,2))+((temp(2)*grid_steps(2))*unit_cell(2,2))
!      else
!        xyzq(nlength+i,1) = minimums(1)+temp(1)*grid_steps(1)
!        xyzq(nlength+i,2) = minimums(2)+temp(2)*grid_steps(2)
      endif
    else
      xyzq(nlength+i,1) = centre_of_mass(1)
      xyzq(nlength+i,2) = centre_of_mass(2)
    endif
    xyzq(nlength+i,3) = centre_of_mass(3)-the_dipole_neutralising_layer
    xyzq(nlength+i,4) = -the_dipole_neutralising_q/q_squared

    temp(1) = temp(1) + 1
    if (temp(1) >= the_dipole_neutralising_density) then
      temp(1) = 0
      temp(2) = temp(2) + 1
    endif
  enddo

  nlength = nlength + q_squared

  return
  end

