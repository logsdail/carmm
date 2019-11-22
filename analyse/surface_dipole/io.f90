! AJL 02/2013
!
! I'm putting all the functions in here to make my life "easier"
!
  integer function atomslength(filename)
!
  character(*)       :: filename
  character(255)     :: line
  integer            :: stat, i
  integer            :: atomcount
  character(2)       :: atoms(10)
  integer            :: quantities(10)
!
  quantities = 0
  atoms = ""
!
  open(unit = 2, file = filename, action = 'READ')

  if (filename.ne."CHGCAR") then
    atomslength = 1
!
    do
      read(2,*,iostat=stat) line
      if (stat /= 0) exit
      atomslength = atomslength + 1
    enddo
!
    atomslength = atomslength - 1
  else
!
    do i=1,5
      read(2,*) line
    enddo
!
! Read in the atomic species
! 
    read(2,'(A255)') line ! atom names
    call tokenise_words(line, atoms)
!
! Count the different types of atoms
!
    atomcount = 0
    do i=1,size(atoms)
      if (atoms(i).ne."") then
        atomcount = atomcount + 1
     endif
    enddo
!
! Read in quantities of each species and sum to get total atoms
!
    read(2,*) quantities(1:atomcount)
    atomslength = 0
    do i=1,atomcount
      atomslength = atomslength + quantities(i)
    enddo
  endif

  close(2)
  return
  end

  integer function chglength(filename,nlength)
!
  character(*)       :: filename
  character(255)     :: line
  integer            :: i
  integer            :: nlength
!
  open(unit = 2, file = filename, action = 'READ')
!
  do i=1,8+nlength
    read(2,*) line
  enddo

  read(2,*) maxx,maxy,maxz

  chglength = maxx*maxy*maxz

  close(2)
  return
  end

  integer function first_line_length(filename)
!
  character(*)       :: filename
  character(255)     :: line
  real               :: values(10)
  integer            :: i
!
  values = 0.0d0
  open(unit = 2, file = filename, action = 'READ')
!
  read(2,'(A255)') line
!
  read (line, *, end=888) (values(i), i = 1, 10)
888 continue
!
  close(2)
!
  if (values(5) .ne. 0.0d0) then
    first_line_length = 5
  else
    first_line_length = 4
  endif
!
  return
  end

  subroutine read_xyzq(nlength, nwidth, xyzq, fractional, fractional_charge)
!
  real               :: xyzq(nlength,nwidth)
  real               :: fractional_charge(nlength)
  logical            :: fractional(3)
  integer            :: nlength, nwidth
  integer            :: i, j, line_length, stat
  integer            :: first_line_length
!
! Check line lengths
!
  line_length = first_line_length("xyzq.txt")
!
  open(unit = 2, file = "xyzq.txt", action = 'READ')

  do i=1,nlength
!    read(2,*) xyzq(i,1:nwidth)
    if (line_length .eq. 5) then 
      read(2,*) xyzq(i,1:nwidth), fractional_charge(i)
    else
      read(2,*) xyzq(i,1:nwidth)
    endif
      
!
! Lets automate the check for if this is fractional coordinates
! For x-, y- and z- coordinates check if they are less than 1.0000001.
! If so, fractional remains true.
!
    do j=1,3
      if (abs(xyzq(i,j)) .gt. 1.0000001d0) then
        fractional(j) = .false.
      endif
    enddo
  enddo

  close(2)
  return
  end

  subroutine read_uc(unit_cell)
!
  real  :: unit_cell(3,3)
!
  open(unit = 2, file = "uc.txt", action = 'READ')
!  read(2,*) unit_cell(1,1), unit_cell(1,2), unit_cell(1,3)
!  read(2,*) unit_cell(2,1), unit_cell(2,2), unit_cell(2,3)
!  read(2,*) unit_cell(3,1), unit_cell(3,2), unit_cell(3,3)
  read(2,*) unit_cell(1,1:3)
  read(2,*) unit_cell(2,1:3)
  read(2,*) unit_cell(3,1:3)
  close(2)

  return
  end

  subroutine read_chgcar(nlength, nwidth, xyzq, chg_length, charge, unit_cell, fractional) 
!
  integer :: nlength, nwidth, chg_length
  integer :: ix, iy, iz, maxx, maxy, maxz  
  real    :: scale 
  real    :: xyzq(nlength, nwidth)
  real    :: charge(chg_length, nwidth)
  real    :: unit_cell(3,3)
  character(255) :: line 
  character(2)   :: atoms(10)
  integer        :: quantities(10)
  integer        :: atomcount, i
  logical        :: fractional(3)
!
  real, allocatable ::CHgrid(:,:,:) ! dynamic 3D matrix 
!
  quantities = 0
  atoms = ""
!
  open(unit = 2, file = "CHGCAR", action = 'READ')
  read(2,*) line !reads title string 
  read(2,*) scale  
  read(2,*) unit_cell(1,1:3)
  read(2,*) unit_cell(2,1:3)
  read(2,*) unit_cell(3,1:3)
!
! Read in the atomic species
! 
  read(2,'(A255)') line ! atom names
  call tokenise_words(line, atoms)
!
! Count the different types of atoms
!
  atomcount = 0

  do i=1,size(atoms)
    if (atoms(i).ne."") then
      atomcount = atomcount + 1
!      write(*,*) i, atoms(i)
   endif 
  enddo 
!
! Read in quantities of each species and make cumulative 
!
  read(2,*) quantities(1:atomcount)
  do i=2,atomcount
    quantities(i) = quantities(i) + quantities(i-1)
  enddo
!
! Read in Direct/Cartesian line
!
  read(2,*) line
  if (trim(line).ne."Direct") then
    fractional = .false.
  endif
!
! Read in atom coordinates
!
  atomcount = 1
  do i=1,nlength 
    read(2,*) xyzq(i,1:3)
!
! This needs some hack to figure out our electronic charge
!
    if (i.gt.quantities(atomcount)) then
      atomcount = atomcount + 1
    endif
 
    xyzq(i,4) = atom_charge(atoms(atomcount))
  enddo
!
!  call print_coords('CHGCAR ATOMIC COORDINATES', nlength, nwidth, xyzq) 
!  
! Read in the charge on the grid
!
  read(2,*) maxx,maxy,maxz 
  allocate(CHGrid(1:maxx,1:maxy,1:maxz)) 
  read(2,*) (((CHGrid(ix,iy,iz),ix=1,maxx),iy=1,maxy),iz=1,maxz) !this reads the entire dataset of charges. 
  CHGrid = - CHGrid / (maxx*maxy*maxz)
 
!  totQ = 0.0d0 
  atomcount = 1
  do ix=1, maxx
    do iy=1,maxy 
       do iz=1,maxz 
!       totQ=totQ+CHgrid(ix,iy,iz)  
        charge(atomcount,1) = real(ix-1)/real(maxx)
        charge(atomcount,2) = real(iy-1)/real(maxy)
        charge(atomcount,3) = real(iz-1)/real(maxz)
        charge(atomcount,4) = CHgrid(ix,iy,iz)
!        write(*,*) charge(atomcount,:)
        atomcount = atomcount + 1
      end do 
    end do
  end do 
  
!  write(*,*) 'The total chargedensity*Volume=', totC 
!  write(*,*) 'The total charge=', totC/(maxx*maxy*maxz) 

!  write(*,*) 'The total charge=', totQ

  deallocate(CHGrid)
  close(2)

  return
  end

  subroutine tokenise_words(line,tokens)
  character(*) :: line
  character(*) :: tokens(10)
  integer      :: pos1, pos2, n
!
  n = 1 
  pos1 = 1
  pos2 = 0
!
!  write(*,*) line
  do
    pos2 = INDEX(line(pos1:), " ")
    if (pos2.eq.0) then
       n = n + 1
       tokens(n) = line(pos1:)
       exit
    endif
    tokens(n) = line(pos1:pos1+pos2-2)
    if (tokens(n).ne." ") then
      n = n + 1
    endif
    pos1 = pos2+pos1
  enddo

  return
  end

!  subroutine get_atomnr(line,atomnr)
!  character(*),intent(IN):: line
!  integer,intent(OUT)      ::atomnr
!
!  integer :: ios,nratoms,atnr,atnrpos
!  character(255):: line2,atl
!
!  !a little bit more complex...read it as a line and sum till it crashes
!  ios=0
!  nratoms=0
!  line2=line
!  do while (ios==0)
!    read(line2,*,IOSTAT=ios) atnr
!    if (ios==0) then
!      nratoms=nratoms+atnr
!      write(atl,*) atnr
!      atl=ADJUSTL(atl)!remove leading blanks
!      atnrpos=Index(line2,trim(atl))+len_trim(atl)
!      line2=" "
!      line2(atnrpos:255)=line(atnrpos:255)
!      if (len_trim(line2)<=0) then
!        ios=10
!      endif
!    endif
!  enddo
!
!  atomnr=nratoms
!
!  return
!  end

  real function atom_charge(atom_label)
!
  character(*)       :: atom_label
!
!  write(*,*) atom_label
  select case(trim(atom_label))
    case("O")
      atom_charge = 6.0d0
    case("Mg")
      atom_charge = 8.0d0
    case("Ca")
      atom_charge = 10.0d0
    case("Sr")
      atom_charge = 10.0d0
    case("Ba")
      atom_charge = 10.0d0
    case default
      write(*,*) atom_label, " IS NOT YET DEFINED. PLEASE EDIT SOURCE CODE"
      stop
  endselect  
!
  return
  end
