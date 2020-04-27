! AJL 02/2013
!
! I'm putting all the functions in here to make my life "easier"
!
  subroutine print_dipoles(nlength, nwidth, xyzq, centre_of_mass, the_dipole)
!
  integer :: i
  integer :: nlength, nwidth
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), the_dipole(3)
!
  real    :: dipole
!
  the_dipole = 0.0d0
!
  write(*,*)
  write(*,*)'DIPOLE'
  write(*,*)'==================='
!
! Loop over all directions
!
  do i=1,3
!
! Calcule dipole
!
    the_dipole(i) = dipole(nlength, nwidth, i, xyzq, centre_of_mass)
!
! Print Dipole
!
    write(*,*) i, ' = ', the_dipole(i), ' e.Ang'
  enddo

  return
  end

  subroutine print_dipole_voltage_offset(unit_cell, the_dipole, the_voltage_offset)
!
  integer :: i, j
  real    :: unit_cell(3,3)
  real    :: the_dipole, the_voltage_offset
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
  write(*,*)
    write(*,*)'DIPOLE VOLTAGE OFFSET'
    write(*,*)'==================='
!
    the_voltage_offset = voltage_offset(unit_cell, the_dipole)

    write(*,*)'VOLTAGE OFFSET = ', the_voltage_offset, ' V'
!
  endif

  return
  end

  subroutine print_quadrupole(nlength, nwidth, xyzq, centre_of_mass, the_quadrupole)
!
  integer :: i, j
  integer :: nlength, nwidth
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), the_quadrupole(3,3), trace
!
  real    :: quadrupole
!
  the_quadrupole = 0.0d0
!
  write(*,*)
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
      the_quadrupole(i,j) = quadrupole(nlength, nwidth, i, j, xyzq, centre_of_mass)
!
! Print Quadrupole
!
      write(*,*) i, ',', j, ' = ', the_quadrupole(i,j), ' e.Ang^2'
    enddo
  enddo

  trace = the_quadrupole(1,1) + the_quadrupole(2,2) + the_quadrupole(3,3)

  write(*,*)'TRACE = ', trace, ' e.Ang^2'

  return
  end

  subroutine print_quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole, the_quadrupole_decomposed)
!
  integer :: nwidth
  real    :: the_quadrupole(3,3), minimums(nwidth), maximums(nwidth), the_quadrupole_decomposed
!
  the_quadrupole_decomposed = 0.0d0
!
  write(*,*)
  write(*,*)'QUADRUPOLE DECOMPOSED'
  write(*,*)'==================='
!
  the_quadrupole_decomposed = quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole) 

  write(*,*)'DIPOLAR CHARGE = ', the_quadrupole_decomposed, ' e.Ang'

  return
  end

  subroutine print_quadrupole_voltage_offset(nwidth, minimums, maximums, centre_of_mass, unit_cell, the_quadrupole_decomposed, &
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
  write(*,*)
    write(*,*)'QUADRUPOLE VOLTAGE OFFSET'
    write(*,*)'==================='
!
    the_voltage_offset = quadrupole_offset(nwidth, minimums, maximums, centre_of_mass, unit_cell, the_quadrupole_decomposed)

    write(*,*)'VOLTAGE OFFSET = ', the_voltage_offset, ' V'
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
  write(*,*)
  write(*,*)'EIGENVALUES OF QUADRUPOLE'
  write(*,*)'==================='
!
  W = 0.0d0
  INFO = 0.0d0
  WORK = 0.0d0
  LWORK = 3 * N
!
  write(*,*) 'DIAGONALISATION DISABLED AS LAPACK IS NOT PRESENT!'
!  call dsyev('V', 'Upper', N, A, LDA, W, WORK, LWORK, INFO)
!
  write(*,*) ' VALUES = ', W(1), ',', W(2), ',', W(3)

  do i=1,3
    write(*,*) 'VECTORS = ', A(i,1), ' , ',  A(i,2), ' , ', A(i,3)
  enddo
!
  return
  end

  subroutine print_octupole(nlength, nwidth, xyzq, centre_of_mass, the_octupole)
!
  integer :: i, j
  integer :: nlength, nwidth
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), the_octupole(3,3,3), trace
!
  real    :: octupole
!
  the_octupole = 0.0d0
!
  write(*,*)
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
        the_octupole(i,j,k) = octupole(nlength, nwidth, i, j, k, xyzq, centre_of_mass)
!
! Print Octupole
!
        write(*,*) i, ',', j, ',', k, ' = ', the_octupole(i,j,k), ' e.Ang^3'
      enddo
    enddo
  enddo

! As taken from "The Theory of Intermolecular Forces", Anthony Stone (Page 17)

  trace = the_octupole(1,1,1) + the_octupole(1,2,2) + the_octupole(1,3,3)

  write(*,*)'TRACE = ', trace, ' e.Ang^3'

  return
  end

  subroutine get_new_q_values_to_coordinates(nlength, nwidth, temp_xyzq, minimums, maximums, centre_of_mass, &
  the_dipole_neutralising_layer, the_dipole_neutralising_q, the_dipole_neutralising_density, unit_cell)
!
  integer :: i, j, q_squared
  integer :: nlength, nwidth, the_dipole_neutralising_density
  real    :: temp_xyzq(nlength,nwidth), minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth)
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
  q_squared = the_dipole_neutralising_density**2
!
! Add points above surface
!
  temp = 0
  do i=1,q_squared
    if (the_dipole_neutralising_density > 1) then
!      if (unit_cell_sum > 0.0d0 .and. .not. fractional) then
      if (unit_cell_sum > 0.0d0) then
        temp_xyzq(i,1) = minimums(1)+((temp(1)*grid_steps(1))*unit_cell(1,1))+((temp(2)*grid_steps(2))*unit_cell(2,1))
        temp_xyzq(i,2) = minimums(2)+((temp(1)*grid_steps(1))*unit_cell(1,2))+((temp(2)*grid_steps(2))*unit_cell(2,2))
!      else
!        temp_xyzq(i,1) = minimums(1)+temp(1)*grid_steps(1)
!        temp_xyzq(i,2) = minimums(2)+temp(2)*grid_steps(2)
      endif
    else
      temp_xyzq(i,1) = centre_of_mass(1)
      temp_xyzq(i,2) = centre_of_mass(2)
    endif
    temp_xyzq(i,3) = centre_of_mass(3)+the_dipole_neutralising_layer
    temp_xyzq(i,4) = the_dipole_neutralising_q/q_squared

    temp(1) = temp(1) + 1
    if (temp(1) >= the_dipole_neutralising_density) then
      temp(1) = 0
      temp(2) = temp(2) + 1
    endif
  enddo
!
  temp = 0
  do i=1,q_squared
    if (the_dipole_neutralising_density > 1) then
!      if (unit_cell_sum > 0.0d0 .and. .not. fractional) then
      if (unit_cell_sum > 0.0d0) then
        temp_xyzq(q_squared+i,1) = minimums(1)+((temp(1)*grid_steps(1))*unit_cell(1,1))+((temp(2)*grid_steps(2))*unit_cell(2,1))
        temp_xyzq(q_squared+i,2) = minimums(2)+((temp(1)*grid_steps(1))*unit_cell(1,2))+((temp(2)*grid_steps(2))*unit_cell(2,2))
!      else
!        temp_xyzq(q_squared+i,1) = minimums(1)+temp(1)*grid_steps(1)
!        temp_xyzq(q_squared+i,2) = minimums(2)+temp(2)*grid_steps(2)
      endif
    else
      temp_xyzq(q_squared+i,1) = centre_of_mass(1)
      temp_xyzq(q_squared+i,2) = centre_of_mass(2)
    endif
    temp_xyzq(q_squared+i,3) = centre_of_mass(3)-the_dipole_neutralising_layer
    temp_xyzq(q_squared+i,4) = -the_dipole_neutralising_q/q_squared

    temp(1) = temp(1) + 1
    if (temp(1) >= the_dipole_neutralising_density) then
      temp(1) = 0
      temp(2) = temp(2) + 1
    endif
  enddo

!
! Debug: Print To Screen
!
!  call print_coords('TO BE ADDED COORDINATES', nlength, nwidth, temp_xyzq)

  return
  end

