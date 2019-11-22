! AJL 04/2014
!
! I'm putting all the functions in here to make my life "easier"
!
  subroutine convert_to_cartesian(nlength, nwidth, xyzq, fractional, unit_cell)
!
  integer               :: nlength, nwidth
  real, dimension(:)    :: xyzq(nlength,nwidth), unit_cell(3,3)
  logical, dimension(:) :: fractional(3)
!
  real    :: temp_coord
  integer :: i
!
  do i=1,nlength
    if (fractional(1) .and. fractional(2)) then
      temp_coord=xyzq(i,1)*unit_cell(1,2)+xyzq(i,2)*unit_cell(2,2)
      xyzq(i,1)=xyzq(i,1)*unit_cell(1,1)+xyzq(i,2)*unit_cell(2,1)
! This makes sure we don't contaminate the components
      xyzq(i,2)=temp_coord
    endif
    if (fractional(3)) then
      xyzq(i,3)=xyzq(i,3)*unit_cell(3,3)
    endif
  enddo
!
  return
  end

  subroutine convert_to_fractional(nlength, nwidth, xyzq, fractional, unit_cell)
!
  integer               :: nlength, nwidth
  real, dimension(:)    :: xyzq(nlength,nwidth), unit_cell(3,3)
  logical, dimension(:) :: fractional(3)
!
  real    :: temp_coord
  integer :: i
!
      do i=1,nlength
        if (fractional(1) .and. fractional(2)) then
          temp_coord=(unit_cell(1,1)*xyzq(i,2) - unit_cell(1,2)*xyzq(i,1)) / &
                     (unit_cell(1,1)*unit_cell(2,2)-unit_cell(1,2)*unit_cell(2,1))
          xyzq(i,1)=(xyzq(i,1)-unit_cell(1,2)*temp_coord)/unit_cell(1,1)
          xyzq(i,2)=temp_coord
        endif
        if (fractional(3)) then
          xyzq(i,3)=xyzq(i,3)/unit_cell(3,3)
        endif
      enddo
!
  return
  end

  subroutine calculate_centre_of_mass(nlength, nwidth, xyzq, minimums, maximums, centre_of_mass, geometric_centre) 
  integer :: i, j
  integer :: nlength, nwidth
  real    :: xyzq(nlength,nwidth), minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth)
  real    :: geometric_centre(nwidth)
!
! Clean arrays
!
  minimums = 10000.0d0
  maximums = -10000.0d0
  centre_of_mass = 0.0d0
  geometric_centre = 0.0d0
!
! Loop over all values
!
  do i=1,nlength
    do j=1,nwidth
!
! Sum averages 
!
      centre_of_mass(j) = centre_of_mass(j) + ( xyzq(i,j) / nlength )
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
! Calculate mid-point between maxima and minima
!
  do j=1,nwidth
    geometric_centre(j) = (minimums(j) + maximums(j)) * 0.5
  enddo

  return
  end

  subroutine restructure_charge(chg_length, nwidth, charge, unit_cell, direction, midpoint)
!
  integer :: chg_length, nwidth, direction
  real    :: charge(chg_length, nwidth), unit_cell(3,3), midpoint
!
  real    :: threshold 
  integer :: i
!
  threshold = midpoint + (unit_cell(direction,direction)/2)
  do i=1,chg_length
    if (charge(i,direction).gt.threshold) then
      charge(i,direction) =  charge(i,direction) - unit_cell(direction,direction)
    endif
  enddo
  
  return
  end

  subroutine splice_coords_in_half(nlength, nwidth, direction, xyzq, unit_cell, threshold, temp_xyzq, new_nlength, splice_top)
!
  integer :: i, counter
  integer :: nlength, nwidth, direction
  real    :: xyzq(nlength, nwidth), temp_xyzq(nlength, nwidth), unit_cell(3,3) 
  real    :: threshold, sumQ
  logical :: splice_top
!
  counter = 1 
  sumQ = 0

!  write(*,*) threshold

  do i=1,nlength
    if (splice_top) then
      if (xyzq(i,direction) .gt. threshold) then
        temp_xyzq(counter,1:nwidth) = xyzq(i,1:nwidth)
        sumQ = sumQ + xyzq(i,4)
        counter = counter + 1
      endif
    else
      if (xyzq(i,direction) .lt. threshold) then
        temp_xyzq(counter,1:nwidth) = xyzq(i,1:nwidth)
        sumQ = sumQ + xyzq(i,4)
        counter = counter + 1
      endif
    endif
  enddo

  write(*,*) 
  write(*,*) 'INTEGRATED CHARGE OF SPLIT CELL = ', sumQ
  write(*,*)
!
! Set new length. We have overcounted by one so accommodate that
! 
  new_nlength = counter - 1
!
! Print Warning if there is some inconsistency in the sizes...
!
  if (new_nlength*2 .ne. nlength) then
    write(*,*)'WARNING: SPLICE SIZE IS NOT HALF THE ORIGINAL SIZE...'
    write(*,*)'BE WARY OF THESE RESULTS AS THEY MAY NOT HOLD TRUE' 
    write(*,*)'ORIGINAL LENGTH = ', nlength
    write(*,*)'NEW LENGTH      = ', new_nlength
!    
! Lets see if we can find the problematic cases (essentially the same Z-coordinate as the cutoff)
! And then we'll add these to the "spliced" xyzq array, but with a halved contribution to the     
!
  endif

  return
  end

  subroutine print_coords(title, nlength, nwidth, xyzq)
  character(*) :: title
  integer        :: i, nlength, nwidth
  real           :: xyzq(nlength,nwidth)
!
  write(*,*)
  write(*,*)
  write(*,*) title, nlength, nwidth
  write(*,*) '==================='
  do i=1,nlength
    write(*,*) xyzq(i,1:nwidth)
  enddo

  return
  end

  subroutine print_min_max_com(nwidth, minimums, maximums, centre_of_mass, geometric_centre)
  integer :: nwidth
  real    :: minimums(nwidth), maximums(nwidth), centre_of_mass(nwidth), geometric_centre(nwidth)
!
  write(*,*)
  write(*,*)'MINIMA'
  write(*,*)'==================='
  write(*,*) minimums(1:nwidth)
  write(*,*)
  write(*,*)'MAXIMA'
  write(*,*)'==================='
  write(*,*) maximums(1:nwidth)
  write(*,*)
  write(*,*)'GEOMETRIC CENTRE'
  write(*,*)'==================='
  write(*,*) geometric_centre(1:nwidth)
  write(*,*)
  write(*,*)'CENTRE OF MASS'
  write(*,*)'==================='
  write(*,*) centre_of_mass(1:nwidth)
  
!
  return
  end

