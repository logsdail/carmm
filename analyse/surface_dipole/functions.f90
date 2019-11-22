! AJL 02/2013
!
! I'm putting all the functions in here to make my life "easier"
!
  real function dipole(nlength, nwidth, direction, xyzq, centre_of_mass) 
  integer :: i
  integer :: nlength, nwidth, direction
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth)
!
! Set initial value
!
  dipole = 0.0d0
!
! Loop and calculate
!
  do i=1,nlength
    dipole = dipole + ((xyzq(i,direction) - centre_of_mass(direction))*xyzq(i,4))
!    dipole = dipole + (xyzq(i,direction)*xyzq(i,4))
  enddo
!
  return
  end

  real function second_moment(nlength, nwidth, direction1, direction2, xyzq, centre_of_mass)
  integer :: i
  integer :: nlength, nwidth, direction1, direction2
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth)
!
! Set initial value
!
  second_moment = 0.0d0
!
! Loop and calculate
!
  do i=1,nlength
    second_moment = second_moment + xyzq(i,4) * (xyzq(i,direction1) - centre_of_mass(direction1)) * &
                 (xyzq(i,direction2) - centre_of_mass(direction2)) 
  enddo
!
  return
  end

  real function third_moment(nlength, nwidth, direction1, direction2, direction3, xyzq, centre_of_mass)
  integer :: i
  integer :: nlength, nwidth, direction1, direction2, direction3
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth)
!
! Set initial value
!
  third_moment = 0.0d0
!
! Loop and calculate
!
  do i=1,nlength
    third_moment = third_moment + xyzq(i,4) * &
                 (xyzq(i,direction1) - centre_of_mass(direction1)) * &
                 (xyzq(i,direction2) - centre_of_mass(direction2)) * &
                 (xyzq(i,direction3) - centre_of_mass(direction3))
  enddo
!
  return
  end

  real function quadrupole(nlength, nwidth, direction1, direction2, xyzq, centre_of_mass)
  integer :: i
  integer :: nlength, nwidth, direction1, direction2
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), distance, kronecker
!
! Set initial value
!
  quadrupole = 0.0d0
  distance = 0.0d0
!
! Define Delta Function
!
  kronecker = 0.0d0
  if (direction1 .eq. direction2) then
    kronecker = 1.0d0
  endif
!
! Loop and calculate
!
  do i=1,nlength
    if (kronecker .eq. 1.0d0) then
      distance = ((xyzq(i,1) - centre_of_mass(1))**2) + ((xyzq(i,2) - centre_of_mass(2))**2) + ((xyzq(i,3) - centre_of_mass(3))**2)
! 
! We use the square of the distance in this so square root is redundant
!      distance = sqrt(distance)
!
    end if

    quadrupole = quadrupole + xyzq(i,4) * &
                 (( 3.0d0 * (xyzq(i,direction1) - centre_of_mass(direction1)) * &
                  (xyzq(i,direction2) - centre_of_mass(direction2))) - distance*kronecker )
  enddo
!
  return
  end

  real function quadrupole_decomposed(nwidth, minimums, maximums, the_quadrupole)
!
! This function decomposes the quadrupole as a linear chain of dipoles in the normal (Z-)
! Q_zz should be equal to 2qa^2, where q is the charge to be decided and a is the half depth of the slab
! Taken from "Molecular Modeling and Simulation of Hydrogen Bonding Pure Fluids and Mixtures", Thorsten Schnabel
!
  integer :: nwidth
  real    :: minimums(nwidth), maximums(nwidth), the_quadrupole(3,3)
  real    :: half_breadth
!
! Set initial value
!
  quadrupole_decomposed = 0.0d0
  half_breadth = 0.0d0
!
  half_breadth = (abs(minimums(3) - maximums(3))) / 2.0d0

  quadrupole_decomposed = the_quadrupole(3,3) / (2.0d0 * (half_breadth**2))
!
  end

  function voltage_offset(unit_cell, the_dipole)
!
! Parameters
!
  real, parameter :: coloumb = 1.602176565d-19
  real, parameter :: angstrom = 1.d-10
  real, parameter :: epsilon0 = 8.854187817d-12
  real, parameter :: pi = 4.d0*datan(1.d0)
!
  real    :: unit_cell(3,3)
  real    :: the_dipole, voltage_offset, surface_area
!
  voltage_offset = 0.0d0
  surface_area   = 0.0d0
!
  surface_area = ( abs(unit_cell(1,1) * unit_cell(2,2)) - abs(unit_cell(1,2) * unit_cell(2,1)) )
!  write(*,*) the_dipole, surface_area
  voltage_offset = the_dipole / surface_area
!
! This is in e.Angs/Angs^2 ==> e/Angs
!
  voltage_offset = voltage_offset * coloumb / angstrom
!
! This is now SI, of C/m. Finally multiply by Coloumb constant.
!
  voltage_offset = voltage_offset / (4 * pi * epsilon0)
  end  

  function quadrupole_offset(nwidth, minimums, maximums, centre_of_mass, unit_cell, the_quadrupole_decomposed)
!
  integer :: nwidth, i
  real    :: minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth), unit_cell(3,3)
  real    :: the_quadrupole_decomposed 
  real    :: half_breadth
  real    :: temp_xyzq(2,nwidth)
!
  real    :: voltage_offset
!
  temp_xyzq = 0.0d0
  quadrupole_offset = 0.0d0
!  surface_area   = 0.0d0
!
  half_breadth = (abs(minimums(3) - maximums(3))) / 2.0d0

  do i=1,2
    temp_xyzq(i,1) = centre_of_mass(1)
    temp_xyzq(i,2) = centre_of_mass(2)
    temp_xyzq(i,3) = ( centre_of_mass(3) + ( half_breadth / 2 ) ) - ( (i-1) * (half_breadth) )
    temp_xyzq(i,4) = (1 - (i-1) * 2 ) * the_quadrupole_decomposed
!    write(*,*) temp_xyzq(i,1), temp_xyzq(i,2), temp_xyzq(i,3), temp_xyzq(i,4)
  enddo 

  quadrupole_offset = voltage_offset(unit_cell,dipole(2, nwidth, 3, temp_xyzq, centre_of_mass))
!
  end

  real function octupole(nlength, nwidth, direction1, direction2, direction3, xyzq, centre_of_mass)
  integer :: i
  integer :: nlength, nwidth, direction1, direction2, direction3
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth), distance
  real    :: kronecker_12, kronecker_13, kronecker_23
!
! Set initial value
!
  octupole = 0.0d0
  distance = 0.0d0
  present = 0.0d0
!
! As taken from "The Theory of Intermolecular Forces", Anthony Stone (Page 17)
!
! Define Delta Function
!
  kronecker_12 = 0.0d0
  kronecker_13 = 0.0d0
  kronecker_23 = 0.0d0
  if (direction1 .eq. direction2) kronecker_12 = 1.0d0
  if (direction1 .eq. direction3) kronecker_13 = 1.0d0
  if (direction2 .eq. direction3) kronecker_23 = 1.0d0
!
! Loop and calculate
!
  do i=1,nlength
    distance = ((xyzq(i,1) - centre_of_mass(1))**2) + ((xyzq(i,2) - centre_of_mass(2))**2) + ((xyzq(i,3) - centre_of_mass(3))**2) 
!
! We use the distance^(2) in this, so no point taking the square root
!
    octupole = octupole + xyzq(i,4) * &
                         ( ( 5.0d0 *   (xyzq(i,direction1) - centre_of_mass(direction1)) * &
                                       (xyzq(i,direction2) - centre_of_mass(direction2)) * &
                                       (xyzq(i,direction3) - centre_of_mass(direction3)) &
                           ) - distance * &
                           ( (xyzq(i,direction1)-centre_of_mass(direction1))*kronecker_23 + &
                             (xyzq(i,direction2)-centre_of_mass(direction2))*kronecker_13 + &
                             (xyzq(i,direction3)-centre_of_mass(direction3))*kronecker_12 &
                         ) ) 
  enddo
!
  return
  end

  real function remove_dipole(nlength, nwidth, direction, xyzq, centre_of_mass, the_dipole_neutralising_layer)
!
  integer :: i 
  integer :: nlength, nwidth, direction
  real    :: xyzq(nlength,nwidth), centre_of_mass(nwidth)
  real    :: temp_xyzq(nlength+2, nwidth)
  real    :: the_dipole_neutralising_layer
  real    :: current_dipole, last_dipole
  real    :: factor
!
! External functions
!
  real    :: dipole
!
! Calculate initial value
!
  i = 1
  factor = 1.0d0
  last_dipole = dipole(nlength, nwidth, direction, xyzq, centre_of_mass)
!
! Turns out I can calculate this with algebra
!
  remove_dipole = -last_dipole/the_dipole_neutralising_layer
  remove_dipole = remove_dipole/2

  write(*,*) 'STARTING DIPOLE = ', last_dipole
!
! Set values to recheck dipole
!
  do i=1,nlength
    temp_xyzq(i,1:nwidth) = xyzq(i,1:nwidth)
  enddo

  temp_xyzq(nlength+1,direction) = centre_of_mass(direction) + the_dipole_neutralising_layer
  temp_xyzq(nlength+1,4)         = remove_dipole 
  temp_xyzq(nlength+2,direction) = centre_of_mass(direction) - the_dipole_neutralising_layer
  temp_xyzq(nlength+2,4)         = -remove_dipole
!
! Calculate new dipole
!
  current_dipole = dipole(nlength+2, nwidth, direction, temp_xyzq, centre_of_mass)
  write(*,*) 'FINAL DIPOLE    = ', current_dipole, 'FINAL CHARGE    = ', remove_dipole

  return
  end
