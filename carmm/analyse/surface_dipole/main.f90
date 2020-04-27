program surface_dipole_sling
!
! Test sling for dipole correction
!
  implicit none
!
! Set constants
!
  integer, parameter :: nwidth = 4 
!
! Allocate variables
!
  integer :: i, j
  integer :: nlength, new_nlength, chg_length
  integer :: stat
  character(32) :: arg
!
! Object carriers
!
  real, allocatable    :: xyzq(:,:)
  real, allocatable    :: fractional_charge(:)
  real, allocatable    :: charge(:,:)
  real, allocatable    :: temp_xyzq(:,:)
  real    :: minimums(nwidth)
  real    :: maximums(nwidth)
  real    :: geometric_centre(nwidth)
  real    :: centre_of_mass(nwidth) 
  real    :: the_dipole(3) ! 1 = x-, 2 = y-, 3 = z-
  real    :: the_dipole_voltage_offset
  real    :: the_second_moment(3,3) ! 1 = x-, 2 = y-, 3 = z-
  real    :: the_quadrupole(3,3) ! 1 = x-, 2 = y-, 3 = z-
  real    :: the_quadrupole_decomposed ! Just Z- axes
  real    :: the_quadrupole_voltage_offset ! Just Z- axes
  real    :: the_third_moment(3,3,3) ! 1 = x-, 2 = y-, 3 = z-
  real    :: the_octupole(3,3,3) ! 1 = x-, 2 = y-, 3 = z-
  real    :: unit_cell(3,3)
  logical :: file_exists
  logical :: fractional(3)
!
  real    :: the_dipole_neutralising_layer
  real    :: the_dipole_neutralising_q
  integer, parameter :: the_dipole_neutralising_density = 4
  integer :: repeat
!
  logical :: verbose 
  logical :: calculate
  logical :: compensate
  logical :: splice
  logical :: splice_twice
  logical :: splice_top
!
! Functions
!
  real    :: remove_dipole
  integer :: atomslength
  integer :: first_line_length
  integer :: chglength
!
! Zero Arrays
!
  unit_cell = 0.0d0
!
  verbose = .false.
  calculate = .false.
  compensate = .false.
  fractional = .true.
  splice = .false.
  splice_twice = .false.
  splice_top = .false.
!
! Set test variables 
!
  i = 0
  do
    call get_command_argument(i, arg)
    if (len_trim(arg) == 0) exit
    if (trim(arg) == "--calculate") calculate = .true.
    if (trim(arg) == "--compensate") compensate = .true.
    if (trim(arg) == "--verbose") verbose = .true.
    if (trim(arg) == "--splice") splice = .true.
    if (trim(arg) == "--splice_twice") splice_twice = .true.      
    if (trim(arg) == "--splice_top") then
      splice = .true. 
      splice_top = .true.
    endif
    i = i+1
  enddo

  inquire(file="CHGCAR", exist=file_exists)
  if (file_exists) then
    nlength = atomslength("CHGCAR")
    allocate(xyzq(nlength,nwidth))
    chg_length = chglength("CHGCAR",nlength)
    allocate(charge(chg_length,nwidth))
    xyzq = 0.0d0
    charge = 0.0d0   
    call read_chgcar(nlength, nwidth, xyzq, chg_length, charge, unit_cell, fractional)
  else 
    inquire(file="xyzq.txt", exist=file_exists)
    if (file_exists) then
      nlength = atomslength("xyzq.txt")
      allocate(xyzq(nlength,nwidth));
      allocate(fractional_charge(nlength))
      xyzq = 0.0d0
      fractional_charge = 0.0d0
      call read_xyzq(nlength, nwidth, xyzq, fractional, fractional_charge)
!
! If no fractional charges are present then we'll deallocate fractional_charge straight away
! In fact, for now we'll just deallocate it anyway. We'll just scale the charges correctly
!
      if (sum(fractional_charge) .ne. 0.0d0) then
        write(*,*) 'FRACTIONAL CHARGES FOUND: PERFORMING SCALING OPERATION'
        do i=1,nlength
          xyzq(i,4) = xyzq(i,4) * fractional_charge(i)
        enddo
      endif
      deallocate(fractional_charge)
!
!
    else
      write(*,*) 'ERROR: XYZQ.TXT AND CHGCAR DO NOT EXIST'
      stop
    endif

    inquire(file="uc.txt", exist=file_exists)
    if (file_exists) then
      call read_uc(unit_cell)
    else
!
! Cannot have fractionals with no unit cell
!
      fractional = .false.
    endif
  endif
!
! Are the x- and y- coordinates fractional?
!
! fractional = .true.
! fractional = .false.
  if (fractional(1) .and. fractional(2)) then
!
! DEBUG: PRINT COORDS
!
    if (verbose) then 
      call print_coords('FRACTIONAL COORDINATES', nlength, nwidth, xyzq)
      write(*,*)
      write(*,*) "Fractional in X- and Y-: TRUE"
      write(*,*) "Converting Coordinates to Cartesian"
    endif
!
    call convert_to_cartesian(nlength, nwidth, xyzq, fractional, unit_cell)
    if (allocated(charge)) call convert_to_cartesian(chg_length, nwidth, charge, fractional, unit_cell)
  else
    if (verbose) then
      write(*,*) "Fractional in X- and Y-: FALSE"
    endif
  endif
!
! Debug: Print To Screen 
!
  if (verbose) then 
    call print_coords('INITIAL UNIT CELL', 3, 3, unit_cell)
    call print_coords('INITIAL COORDINATES', nlength, nwidth, xyzq)
  endif
!
! Calculate COM
!
  call calculate_centre_of_mass(nlength, nwidth, xyzq, minimums, maximums, centre_of_mass, geometric_centre)
!
! Debug: Print To Screen
!
  call print_min_max_com(nwidth, minimums, maximums, centre_of_mass, geometric_centre)
  if (allocated(charge)) then
    call restructure_charge(chg_length, nwidth, charge, unit_cell, 3, geometric_centre(3))
!
! Debug process to check the charge is correctly altered
!
!    call calculate_centre_of_mass(chg_length, nwidth, charge, minimums, maximums, centre_of_mass, geometric_centre)
!    call print_min_max_com(nwidth, minimums, maximums, centre_of_mass, geometric_centre)
!
! Reset the arrays I've just ruined!
!
!    call calculate_centre_of_mass(nlength, nwidth, xyzq, minimums, maximums, centre_of_mass, geometric_centre)
  endif
!
! Redo if we are going to splice the layer in half.
!
  if (splice_twice) then
    repeat = 2
  else if (splice) then 
    repeat = 1
  else 
    repeat = 0
  endif

  do i=1,repeat
!    if (verbose) then
      write(*,*)
      write(*,*)'SPLICING UNIT CELL IN HALF IN THE Z-DIRECTION'
!    endif

    allocate(temp_xyzq(nlength,nwidth))
    call splice_coords_in_half(nlength, nwidth, 3, xyzq, unit_cell, geometric_centre(3), temp_xyzq, new_nlength, splice_top)
    deallocate(xyzq)

    nlength = new_nlength
    allocate(xyzq(nlength,nwidth))

    do j=1,nlength
      xyzq(j,1:nwidth) = temp_xyzq(j,1:nwidth)
    enddo

    deallocate(temp_xyzq)
!
! We need to repeat this for charge if it is allocated
!
    if (allocated(charge)) then
      if (verbose) then
        write(*,*)
        write(*,*)'SPLICING CHARGE IN HALF IN THE Z-DIRECTION'
      endif
      allocate(temp_xyzq(chg_length,nwidth))
      call splice_coords_in_half(chg_length, nwidth, 3, charge, unit_cell, geometric_centre(3), temp_xyzq, new_nlength, & 
                                 splice_top)
      deallocate(charge)

      chg_length = new_nlength
      allocate(charge(chg_length,nwidth))

      do j=1,chg_length
        charge(j,1:nwidth) = temp_xyzq(j,1:nwidth)
      enddo

      deallocate(temp_xyzq)
    endif

!    splice_top = .not.splice_top
!
! Debug: Print To Screen
!
    if (verbose) then
      call print_coords('SPLICED UNIT CELL', 3, 3, unit_cell)
      call print_coords('SPLICED COORDINATES', nlength, nwidth, xyzq)
    endif

    call calculate_centre_of_mass(nlength, nwidth, xyzq, minimums, maximums, centre_of_mass, geometric_centre)
    call print_min_max_com(nwidth, minimums, maximums, centre_of_mass, geometric_centre)
  enddo

! We now need to combine in the CHGCAR density around the atomic centres, if allocated
!
  if (allocated(charge) .and. .not. compensate) then
!
! We can just add the charge density in with the other information on atomic centres
!
    allocate(temp_xyzq(nlength,nwidth))
    do i=1,nlength
      temp_xyzq(i,1:nwidth) = xyzq(i,1:nwidth)
    enddo
    deallocate(xyzq)

    allocate(xyzq(nlength+chg_length,nwidth))
    do i=1,nlength
      xyzq(i,1:nwidth) = temp_xyzq(i,1:nwidth)
    enddo
    do i=1,chg_length
      xyzq(nlength+i,1:nwidth) = charge(i,1:nwidth)
    enddo

    nlength = nlength + chg_length
!
! Tidy up
!
    deallocate(temp_xyzq)
    deallocate(charge)
!
  endif
!
! Calculate all the values now
!
  call print_dipoles(nlength, nwidth, xyzq, geometric_centre, the_dipole)
  call print_dipole_voltage_offset(unit_cell, the_dipole(3), the_dipole_voltage_offset)

  if (calculate) then
!    call print_second_moment(nlength, nwidth, xyzq, centre_of_mass, the_second_moment)
!    call print_quadrupole(nlength, nwidth, xyzq, centre_of_mass, the_quadrupole)
!    call print_quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole, the_quadrupole_decomposed)
!    call print_voltage_offset(nwidth, minimums, maximums, centre_of_mass, unit_cell, the_quadrupole_decomposed, the_voltage_offset)
!    call print_diagonalised_quadrupole(the_quadrupole)
!    call print_third_moment(nlength, nwidth, xyzq, centre_of_mass, the_third_moment)
!    call print_octupole(nlength, nwidth, xyzq, centre_of_mass, the_octupole)
!    if (splice) call print_dipole_voltage_offset(unit_cell, the_dipole(3), the_dipole_voltage_offset)
    call print_second_moment(nlength, nwidth, xyzq, geometric_centre, the_second_moment)
    call print_quadrupole(nlength, nwidth, xyzq, geometric_centre, the_quadrupole)
    call print_quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole, the_quadrupole_decomposed)
    call print_quadrupole_voltage_offset(nwidth, minimums, maximums, geometric_centre, unit_cell, &
                              the_quadrupole_decomposed, the_quadrupole_voltage_offset)
    call print_diagonalised_quadrupole(the_quadrupole)
    call print_third_moment(nlength, nwidth, xyzq, geometric_centre, the_third_moment)
    call print_octupole(nlength, nwidth, xyzq, geometric_centre, the_octupole)
  endif

  if (compensate) then
!
! Editing just z- direction: 3
!
    write(*,*)
    write(*,*)
!
! Calculate distance from COM for dipole neutralising layer
! Currently dist(COM->Max/Min)*5
!
!    if (abs(centre_of_mass(3)-maximums(3)) > abs(centre_of_mass(3)-minimums(3))) then  
!      the_dipole_neutralising_layer = abs(centre_of_mass(3)-maximums(3))*5
!    else
!      the_dipole_neutralising_layer = abs(centre_of_mass(3)-minimums(3))*5
!    endif

    the_dipole_neutralising_layer = abs(geometric_centre(3)-maximums(3))*5 

!
! Print position of dipole neutralising layer
!
    write(*,*) 'PLANES FOR CHARGE NEUTRALISING LAYERS (Z-AXIS)'
    write(*,*) '==================='
!    write(*,*) centre_of_mass(3)+the_dipole_neutralising_layer, ' AND ', &
!    centre_of_mass(3)-the_dipole_neutralising_layer
    write(*,*) geometric_centre(3)+the_dipole_neutralising_layer, ' AND ', &
    geometric_centre(3)-the_dipole_neutralising_layer
    write(*,*)
    
!
! Find out necessary charge at this distance to remove dipole
!
!    the_dipole_neutralising_q = remove_dipole(nlength, nwidth, 3, xyzq, centre_of_mass, the_dipole_neutralising_layer)
    the_dipole_neutralising_q = remove_dipole(nlength, nwidth, 3, xyzq, geometric_centre, the_dipole_neutralising_layer)
!
! It would be better if we passed the minimums and maximums for the unit cell in
! here, but this will do for now. The add_new_q routine will need editing when
! updated for the unit cell. Otherwise atoms will be distributed with respect to the current atoms
!
! We need to update size of the arrays first.
! Firstly, we back up the old data
!
    allocate(temp_xyzq(nlength,nwidth))
    temp_xyzq(:,:) = xyzq(:,:)
!
! Not the prettiest but move_alloc was proving troublesome
!
    new_nlength = nlength+(2*(the_dipole_neutralising_density**2))
    deallocate(xyzq) 
    allocate(xyzq(new_nlength,nwidth))
    xyzq=0.0d0

    do i=1,nlength
        xyzq(i,1:nwidth) = temp_xyzq(i,1:nwidth)
    enddo

!
! Debug: Print To Screen
!
    if (verbose) then
      call print_coords('COPIED COORDINATES', nlength, nwidth, xyzq)
    endif
!
    deallocate(temp_xyzq)
    allocate(temp_xyzq(new_nlength-nlength,nwidth))
!
!
!    call get_new_q_values_to_coordinates(new_nlength-nlength, nwidth, temp_xyzq, minimums, maximums, centre_of_mass, &
    call get_new_q_values_to_coordinates(new_nlength-nlength, nwidth, temp_xyzq, minimums, maximums, geometric_centre, &
    the_dipole_neutralising_layer, the_dipole_neutralising_q, the_dipole_neutralising_density, unit_cell)
!
! Combine with previous xyzq values and update nlength
!
    if (verbose) then
      call print_coords('TO BE ADDED COORDINATES', new_nlength-nlength, nwidth, temp_xyzq)
    endif

    do i=1,new_nlength-nlength
      xyzq(nlength+i,1:nwidth) = temp_xyzq(i,1:nwidth)
    enddo

    deallocate(temp_xyzq)
    nlength=new_nlength
!
! Debug: Print To Screen 
!
    if (verbose) then 
      call print_coords('FINAL UNIT CELL', 3, 3, unit_cell)
      call print_coords('FINAL COORDINATES', nlength, nwidth, xyzq)
    endif
!
! Calculate COM
!
    call calculate_centre_of_mass(nlength, nwidth, xyzq, minimums, maximums, centre_of_mass, geometric_centre)
!
! Debug: Print To Screen
!
    call print_min_max_com(nwidth, minimums, maximums, centre_of_mass, geometric_centre)

!    call print_dipoles(nlength, nwidth, xyzq, centre_of_mass, the_dipole)
    call print_dipoles(nlength, nwidth, xyzq, geometric_centre, the_dipole)

    if (calculate) then
!      call print_second_moment(nlength, nwidth, xyzq, centre_of_mass, the_second_moment)
!      call print_quadrupole(nlength, nwidth, xyzq, centre_of_mass, the_quadrupole)
!      call print_third_moment(nlength, nwidth, xyzq, centre_of_mass, the_third_moment)
!      call print_octupole(nlength, nwidth, xyzq, centre_of_mass, the_octupole)
      call print_second_moment(nlength, nwidth, xyzq, geometric_centre, the_second_moment)
      call print_quadrupole(nlength, nwidth, xyzq, geometric_centre, the_quadrupole)
      call print_third_moment(nlength, nwidth, xyzq, geometric_centre, the_third_moment)
      call print_octupole(nlength, nwidth, xyzq, geometric_centre, the_octupole)
    endif
!
! Are the x- and y- coordinates originally fractional?
! If so convert back to fractional coordinates
!
    if (fractional(1) .and. fractional(2)) then
!      call print_coords('UNIT CELL', 3, 3, unit_cell)
      call convert_to_fractional(nlength,nwidth,xyzq,fractional,unit_cell)
    !
    ! Debug: Print To Screen
    !
      if (verbose) then 
        write(*,*)
        write(*,*) "Converting Coordinates to Fractional"
        call print_coords('FRACTIONAL COORDINATES', nlength, nwidth, xyzq)
      endif
    endif
  
  endif
!
  write(*,*)
  deallocate(xyzq) 
!
end program surface_dipole_sling
