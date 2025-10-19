!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C) 2017 Jorge Cebola Borbinha, 
! Faculty of Sciences and Technology, NOVA University of Lisbon
!
! @author:  Jorge Cebola Borbinha
! @github:  jorge-borbinha
! @website: jorge-borbinha.github.io
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the MIT License.
!
! For more information, see https://choosealicense.com/licenses/mit/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! ==============================================================================
!
! MODULE: PhantomGlobalVars
!   This module encapsulates all global data and configuration parameters for 
!   the ReDimLungs program. Using a module allows for easy sharing of variables 
!   between the main program and subroutines without passing them as arguments 
!   repeatedly.
!
! ==============================================================================
MODULE PhantomGlobalVars
  implicit none
  save

  ! --- Configuration parameters loaded from input file (reDimLungs.in) ---
  integer            :: organ_choice
  character(len=256) :: phantom_file
  character(len=256) :: tasks_file
  integer            :: tasks_file_lines
  character(len=256) :: organlist_file
  integer            :: organlist_file_lines
  integer :: n_slices, n_rows, n_cols
  real    :: x_res, y_res, z_res
  integer :: right_lung_id, right_lung_blood_id
  integer :: left_lung_id, left_lung_blood_id
  integer :: stomach_id, esophagous_id
  integer :: heart_blood_id
  integer :: min_bone_id, max_bone_id
  integer :: slice_tolerance

  ! --- Global Arrays and Variables ---
  ! Main 3D array to hold the entire phantom data
  integer, allocatable, dimension(:,:,:) :: phantom_vector
  ! 3D array to hold the smaller, isolated section of the organ for modification
  integer, allocatable, dimension(:,:,:) :: organ_vector
  ! 2D arrays to store the min/max coordinates of the organ along each axis
  integer, allocatable, dimension(:,:)   :: x_vector, y_vector, z_vector

  ! --- File and UI related variables ---
  character(len=256) :: phantom_out_file
  character(len=70)  :: fwd_line = '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  character(len=11)  :: fwd_title_prompt = '>>>>>>>>>> '
  character(len=6)   :: section_prompt = '>>>>> '
  character(len=6)   :: line_prompt = ' >    '

END MODULE PhantomGlobalVars


! ==============================================================================
!
! PROGRAM: ReDimLungs
!
! PURPOSE:
!   Main driver for the ReDimLungs program. This program semi-automatically
!   modifies a 3D reference voxel phantom to more accurately match the anatomy 
!   of an individual patient, mainly focusing on scaling the lungs.
!
! WORKFLOW:
!   - Loads configuration parameters from an input file (reDimLungs.in).
!   - Reads the input phantom data and identifies the initial dimensions
!     of the target organ (left or right lung).
!   - Applies a series of modifications (adding/deleting voxel layers) based
!     on a user-defined task list.
!   - Writes the newly modified phantom to output files suitable for
!     visualization and Monte Carlo simulation.
!
! AUTHOR:
!   Jorge Cebola Borbinha
!
! DATE:
!   March 2020
! ==============================================================================
PROGRAM ReDimLungs
  use PhantomGlobalVars
  implicit none

  ! --- Local Variables ---
  integer :: organ_id, other_lung_id, lung_blood_id
  integer :: n_voxels, ios
  logical :: ex
  character(len=256) :: organ_choice_name

  ! Declare a namelist to read parameters using a namelist. 
  ! More reliable than formatted read for this purpose.
  namelist /PHANTOM_CONFIG/ organ_choice, phantom_file, &
                            tasks_file, tasks_file_lines, &
                            organlist_file, organlist_file_lines, &
                            n_slices, n_rows, n_cols, &
                            x_res, y_res, z_res, &
                            right_lung_id, left_lung_id, & 
                            right_lung_blood_id, left_lung_blood_id, &
                            stomach_id, esophagous_id, &
                            heart_blood_id, min_bone_id, max_bone_id, &
                            slice_tolerance

  ! Initialize report file and console output
  open (unit=10, file='reDimLungs.out', status='replace')
  call WriteHeader()

  inquire(file='reDimLungs.in', exist=ex)
    if (.not. ex) then
      call HandleError("Input file 'reDimLungs.in' not found.")
    end if

    open(unit=99, file='reDimLungs.in', status='old', iostat=ios)
    if (ios /= 0) call HandleError("Could not open input file 'reDimLungs.in'.")

    read(99, nml=PHANTOM_CONFIG, iostat=ios)
    if (ios /= 0) then
        close(99)
        call HandleError("Error reading or parsing 'reDimLungs.in'. Check the format.")
    end if
    close(99)

    write(*,*) " Input file loaded successfully."
    write(10,*) section_prompt, "Input file 'reDimLungs.in' read successfully."
    write(10,*) " "

  if (organ_choice==1) then
      organ_choice_name = 'Right Lung'
      organ_id = right_lung_id
      other_lung_id = left_lung_id
      lung_blood_id = right_lung_blood_id
  else if (organ_choice==2) then
      organ_choice_name = 'Left Lung'
      organ_id = left_lung_id
      other_lung_id = right_lung_id
      lung_blood_id = left_lung_blood_id
  end if
  write(*,*) ' The organ to redimension is: ', trim(organ_choice_name)
  write(10,*)' The organ to redimension is: ', trim(organ_choice_name)
  write(*,*)' '

  n_voxels = n_slices * n_rows * n_cols

  ! Allocate global arrays now that dimensions are known
  allocate(phantom_vector(n_cols, n_rows, n_slices), stat=ios)
  if (ios /= 0) call HandleError("Failed to allocate memory for the main phantom_vector.")

  ! - Read phantom data and get initial organ dimensions
  call ReadPhantomDataAndGetDimensions(organ_id, organ_choice_name, n_voxels)

  ! - Apply modifications based on the tasks file
  call ApplyModifications(organ_id, other_lung_id, lung_blood_id, tasks_file_lines)

  ! - Write the final modified phantom to output files
  call WriteOutputFiles(organ_id, lung_blood_id, n_voxels)

  ! --- Finalization ---
  write(10,*)fwd_line
  write(10,*)'Program finished successfully.'
  close(10)
  write(*,*) 'Program finished successfully. Check reDimLungs.out for details.'

CONTAINS

  ! ==============================================================================
  ! SUBROUTINE: WriteHeader
  ! PURPOSE: Writes the initial header to the console and report file.
  ! ==============================================================================
  subroutine WriteHeader()
    ! Write header to report file
    write(10,*) fwd_line
    write(10,*) fwd_line
    write(10,*) fwd_title_prompt,'      ReDimLungs Fortran Routine       ',fwd_title_prompt
    write(10,*) fwd_line
    write(10,*) fwd_line
    write(10,*) '  '
    write(10,*) '  '
    ! Write header to console
    write(*,*)' >>> ReDimLungs Fortran Routine <<< '
    write(*,*)'  '
  end subroutine WriteHeader

  ! ==============================================================================
  ! SUBROUTINE: HandleError
  ! PURPOSE: 
  !   Prints an error message and terminates the program with informative
  !   error message.
  !
  ! ARGUMENTS:
  !   message (input, character): The error message to display.
  ! ==============================================================================
  subroutine HandleError(message)
    character(len=*), intent(in) :: message
    write(*,*) " "
    write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*) "! ERROR:"
    write(*,*) "! ", trim(message)
    write(*,*) "! Program will now terminate."
    write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    stop 1
  end subroutine HandleError

  ! ==============================================================================
  ! SUBROUTINE: ReadProgress
  ! PURPOSE: Displays a progress indicator to the console.
  ! ARGUMENTS: z (input, integer): slice index
  ! ==============================================================================
  subroutine ReadProgress(z)
    use PhantomGlobalVars, only: n_slices
    integer, intent(in) :: z
    real :: percentage
    integer, save :: last_percentage_shown = 0

    percentage = real(z) / real(n_slices) * 100.0

    ! print header for progress bar, use advance='no' to stay on the same line. Reset counter
    if (z == 1) then
      last_percentage_shown = 0
      write(*, fmt='(a)', advance='no') '  Progress: '
    end if
    ! Check if the percentage has reached the next 10% milestone and hasn't been displayed before.
    if (int(percentage) / 10 > last_percentage_shown / 10) then
      write(*, fmt='(I3,A)', advance='no') int(percentage), '%, '
      last_percentage_shown = int(percentage) !update last_percantage_shown to prevent repeated printing
    end if
    ! handles final slice, ensures progess bar is completes and reports n_slices
    if (z == n_slices) then
      write(*, fmt='(A,I3)', advance='yes') 'Complete. Total slices processed: ', z
    end if

    ! if (z==1)  then
    !   write(*,fmt='(a)', advance='no') '  No of slices processed:'
    ! else if (mod(z, 50) == 0) then
    !   write(*,fmt='(i3,a)', advance='no') z, ','
    ! else if (z==n_slices) then
    !   write(*,fmt='(i3)') z
    ! end if
  end subroutine ReadProgress

! ==============================================================================
  ! SUBROUTINE: SkipBOMAndTabs
  ! PURPOSE:
  !   Skips the UTF-8 Byte Order Mark (BOM) and any leading tab characters at the
  !   beginning of a file. This ensures that file reads occur without errors.
  !   The UTF-8 BOM is a specific sequence of three bytes. The BOM bytes are 
  !   EF BB BF in hexadecimal.
  !
  ! ARGUMENTS:
  !   unit (input, integer): The file unit number to process.
  ! ==============================================================================
  subroutine SkipBOMAndTabs(unit)
    implicit none
    integer, intent(in) :: unit
    character(len=3) :: bom_buffer
    character(len=1) :: c
    integer :: ios

    ! Check for bom character sequence at the beginning of the file.
    ! Note: ' ' is the unicode representation of the BOM character.
    character(len=3) :: utf8_bom
    utf8_bom = ' '  ! This is the UTF-8 BOM character sequence.

    ! Read the first 3 characters into a buffer.
    read(unit, fmt='(a3)', iostat=ios, advance='no') bom_buffer
    
    ! Check if the buffer equals the UTF-8 BOM.
    if (ios == 0 .and. bom_buffer == utf8_bom) then
      ! If a BOM was detected, rewind the file pointer to the beginning of the file.
      rewind(unit=unit)
      ! Now, read and discard the 3 BOM characters.
      read(unit, '(a3)')
    else
      ! If no BOM was detected, just rewind to the beginning.
      rewind(unit=unit)
    end if

    ! Skip any leading tab characters.
    do
      read(unit, fmt='(a)', iostat=ios, advance='no') c
      if (ios /= 0) then
        ! Error or end of file reached.
        exit
      end if
      if (c == char(9)) then
        ! Backspace to stay on the tab character and then read again to consume it.
        backspace(unit)
        read(unit, '(a1)')
      else
        ! Not a tab, this is the first real character. Backspace to reposition the file pointer for the next read statement.
        backspace(unit)
        exit
      end if
    end do
  end subroutine SkipBOMAndTabs

    ! ==============================================================================
  ! SUBROUTINE: CalculateLungDimensions
  ! PURPOSE:
  !   Calculates and reports the dimensions of the currently selected lung
  !   in both voxels and centimeters, using a specific output format.
  !
  ! ARGUMENTS:
  !   min_x (input, integer) : Minimum x-coordinate of the lung.
  !   max_x (input, integer) : Maximum x-coordinate of the lung.
  !   min_y (input, integer) : Minimum y-coordinate of the lung.
  !   max_y (input, integer) : Maximum y-coordinate of the lung.
  !   min_z (input, integer) : Minimum z-coordinate of the lung.
  !   max_z (input, integer) : Maximum z-coordinate of the lung.
  ! ==============================================================================
  subroutine CalculateLungDimensions(min_slice, max_slice, organ_id_in)
    use PhantomGlobalVars

    implicit none
    integer, intent(in) :: min_slice, max_slice, organ_id_in

    integer :: x, y, z
    integer :: min_x, max_x, abs_min_x, abs_max_x 
    integer :: min_y, max_y, abs_min_y, abs_max_y
    integer :: min_z, max_z, abs_min_z, abs_max_z
    integer :: counter_x, counter_y, counter_z


        ! --- X-Axis Analysis ---
    do z=min_slice,max_slice
      do y=1,n_rows
        counter_x = y + n_rows*(z-1)
        max_x = 0
        min_x = n_cols + 1
        do x=1,n_cols
          if (organ_vector(x,y,z) == organ_id_in) then
            if (x < min_x) min_x = x
            if (x > max_x) max_x = x
          end if
        end do
        x_vector(counter_x,1) = z
        x_vector(counter_x,2) = y
        x_vector(counter_x,3) = min_x
        x_vector(counter_x,4) = max_x
      end do
    end do
    abs_min_x=minval(x_vector(:,3))
    abs_max_x=maxval(x_vector(:,4))
    !write (*,*)'  '
    write (*,*)' > X AXIS LUNG COORDINATES:'
    write (10,*)' > X AXIS LUNG COORDINATES:'
    write (*,'(A,I3,2X,A,I3)') '  Min x :', abs_min_x, '  Max x:', abs_max_x
    write (10,'(A,I3,2X,A,I3)')'  Min x :', abs_min_x, '  Max x:', abs_max_x
    write (*,'(A,I3,A,F7.4,A)') '  Size:', abs_max_x-abs_min_x+1,' slices; ',(abs_max_x-abs_min_x+1)*x_res,' cm.'
    write (10,'(A,I3,A,F7.4,A)')'  Size:', abs_max_x-abs_min_x+1,' slices; ',(abs_max_x-abs_min_x+1)*x_res,' cm.'

    ! --- Y-Axis Analysis ---
    do z=min_slice,max_slice
      do x=1,n_cols
        counter_y=x+n_cols*(z-1)
        max_y=0
        min_y=n_rows+1
        do y=1,n_rows
          if (organ_vector(x,y,z) == organ_id_in) then
            if (y < min_y) min_y = y
            if (y > max_y) max_y = y
          end if
        end do
        y_vector(counter_y,1) = z
        y_vector(counter_y,2) = x
        y_vector(counter_y,3) = min_y
        y_vector(counter_y,4) = max_y
      end do
    end do
    abs_min_y=minval(y_vector(:,3))
    abs_max_y=maxval(y_vector(:,4))
    write(*,*)' > Y AXIS LUNG COORDINATES:'
    write(10,*)' > Y AXIS LUNG COORDINATES:'
    write(*,'(A,I3,2X,A,I3)') '  Min y :', abs_min_y, '  Max y:', abs_max_y
    write(10,'(A,I3,2X,A,I3)')'  Min y :', abs_min_y, '  Max y:', abs_max_y
    write(*,'(A,I3,A,F7.4,A)') '  Size:', abs_max_y-abs_min_y+1,' slices; ',(abs_max_y-abs_min_y+1)*y_res,' cm.'
    write(10,'(A,I3,A,F7.4,A)')'  Size:', abs_max_y-abs_min_y+1,' slices; ',(abs_max_y-abs_min_y+1)*y_res,' cm.'

    ! --- Z-Axis Analysis ---
    do y=1,n_rows
      do x=1,n_cols
        counter_z=x+n_cols*(y-1)
        max_z=0
        min_z=n_slices+1
        do z=min_slice,max_slice
          if (organ_vector(x,y,z) == organ_id_in) then
            if (z < min_z) min_z = z
            if (z > max_z) max_z = z
          end if
        end do
        z_vector(counter_z,1) = y
        z_vector(counter_z,2) = x
        z_vector(counter_z,3) = min_z
        z_vector(counter_z,4) = max_z
      end do
    end do
    abs_min_z = minval(z_vector(:,3))
    abs_max_z = maxval(z_vector(:,4))
    write(*,*)' > Z AXIS LUNG COORDINATES:'
    write(10,*)' > Z AXIS LUNG COORDINATES:'
    write(*,'(A,I3,2X,A,I3)') '  Min z :', abs_min_z, '  Max z :', abs_max_z
    write(10,'(A,I3,2X,A,I3)')'  Min z :', abs_min_z, '  Max z :', abs_max_z
    write(*,'(A,I3,A,F7.4,A)') '  Size:', abs_max_z-abs_min_z+1,' slices; ',(abs_max_z-abs_min_z+1)*z_res,' cm.'
    write(10,'(A,I3,A,F7.4,A)')'  Size:', abs_max_z-abs_min_z+1,' slices; ',(abs_max_z-abs_min_z+1)*z_res,' cm.'

  end subroutine CalculateLungDimensions

  ! =================================================================================
  ! SUBROUTINE: ReadPhantomDataAndGetDimensions
  !
  ! PURPOSE:
  !   Handles initial data loading and analysis. This subroutine reads the phantom 
  !   data into a 3D array, from a specified ct-den-mat type file. Despite taking 
  !   ct-den-mat files taking more time to read, it's necessary for this routine 
  !   to input and output ct-den-mat files, because: 
  !     1) '.vox' files do not contain organ IDs, only material IDs and densities; 
  !     2) '.dat' raw data files are not standardized, may vary in format and do 
  !        not contain voxel coordinates; 
  !     3) ct-den-mat type files are used to visualize the phantom (employing gnuplot)
  !        in an iterative cycle with this routine until the target size for the
  !        lungs is achieved in all directions (x,y and z). Therefore, the ct-den-mat 
  !        file must be used as input and output format.
  !   Furthermore, this subroutine isolates the relevant organ (i.e. lungs) section
  !   (adding some tolerance for possible layer additions), calculates and reports the 
  !   initial 3D dimensions (bounding box) of the target organ before any modifications
  !   are made.
  !
  ! ARGUMENTS:
  !   organ_id_in      (input, integer): The numerical ID of the target organ.
  !   organ_name_in    (input, character): The name of the organ for reporting.
  !   x_res, y_res, z_res (input, real): Voxel resolutions in cm.
  !   n_voxels_in      (input, integer): Total number of voxels in the phantom.
  !
  ! =================================================================================
  subroutine ReadPhantomDataAndGetDimensions(organ_id_in, organ_name_in, n_voxels_in)
    use PhantomGlobalVars
    implicit none

    ! --- Arguments ---
    integer, intent(in)           :: organ_id_in
    character(len=*), intent(in)  :: organ_name_in
    integer, intent(in)           :: n_voxels_in

    ! --- Local Variables ---
    integer :: x, y, z, x_read, y_read, z_read, organ, material
    real    :: density
    integer :: min_z, max_z, min_slice, max_slice
    integer :: i, ios
    logical :: ex
    character(len=256) :: error_msg

    ! --- SECTION: Load Phantom File ---
    write(*,*)' Reading the phantom ct-den-mat type file...'
    inquire(file=trim(phantom_file), exist=ex)
    if (.not. ex) call HandleError("Phantom file not found: " // trim(phantom_file))
    open(unit=11, file=trim(phantom_file), status='old', iostat=ios)
    if (ios /= 0) call HandleError("Error opening phantom file: " // trim(phantom_file))
    
    ! Patch to skip BOM and leading tabs before reading the file content
    call SkipBOMAndTabs(11)

    do i=1,11 ! Read and discard header lines
      read (11,*,iostat=ios) 
      if (ios /= 0) call HandleError("Error reading header from phantom file.")
    end do

    do z=1,n_slices
      do y=1,n_rows
        do x=1,n_cols
            read (11,*,iostat=ios) x_read,y_read,z_read,density,material,organ
            if (ios > 0) then
              write(error_msg,'(A,I3,A,I3,A,I3,A)')"Formatting error while reading phantom file&
                               & at voxel (x,y,z) = (",x_read,",",y_read,",",z_read,")."
              call HandleError(trim(error_msg))
            else if (ios < 0) then
              write(error_msg,'(A,I3,A,I3,A,I3,A)')"Unexpected end of file reached while reading phantom data&
                               & at voxel (x,y,z) = (",x_read,",",y_read,",",z_read,")."
              call HandleError(trim(error_msg))
            end if
            phantom_vector(x,y,z) = organ
        end do
        read(11,*,iostat=ios) ! Skip blank line
      end do
      read(11,*,iostat=ios) ! Skip blank line
      call ReadProgress(z)
    end do
    close(11)

    write(*,*)' Phantom written in a 3D vector.'
    write(*,*)' Detecting organ section...'

    ! Write details to report file
    write(10,*)fwd_line
    write(10,*)section_prompt,'SECTION: Load phantom file and detect organ section.'
    write(10,*)'  '
    write(10,*)section_prompt,'Load Phantom File'
    write(10,*)'Opening phantom file:'
    write(10,*)line_prompt,trim(phantom_file)
    write(10,*)'Organ considered for redimension and organ ID:'
    write(10,'(A,A,I3)')line_prompt,trim(organ_name_in), organ_id_in
    write(10,*)'Voxel dimensions in x,y,z (cm):'
    write(10,'(A, 3(F6.4,1x))')line_prompt,x_res,y_res,z_res
    write(10,*)'Voxels volume (cm^3):'
    write(10,'(A,F8.6)')line_prompt,x_res*y_res*z_res
    write(10,*)'No. of voxels in x,y,z and total:'
    write(10,'(A,3(I4,1X),2X,I9)')line_prompt,n_cols,n_rows,n_slices,n_voxels_in
    write(10,*)'Reading the phantom file...'
    write(10,*)'Phantom writen in a 3D vector.'
    write(10,*)'  '

    ! --- SECTION: Isolate Organ Section ---
    max_z = 0
    min_z = n_slices + 1
    do x=1,n_cols
      do y=1,n_rows
        do z=1,n_slices
          if (phantom_vector(x,y,z) == organ_id_in) then
            if (z < min_z) min_z = z
            if (z > max_z) max_z = z
          end if
        end do
      end do
    end do

    max_slice = max_z + slice_tolerance
    min_slice = min_z - slice_tolerance
    if (min_slice < 1) min_slice = 1
    if (max_slice > n_slices) max_slice = n_slices

    write(*,*)' > Organ Section Detected:'
    write(*,*)' Min and Max z values:'
    write(*,'(A,I4,2X,I4)') line_prompt, min_z, max_z
    write(*,*)' Section considered for redimension:'
    write(*,'(A,I4,2X,I4)') line_prompt,min_slice,max_slice
    write(*,*)'  '

    allocate(organ_vector(n_cols, n_rows, min_slice:max_slice), stat=ios)
    if (ios /= 0) call HandleError("Failed to allocate memory for the organ_vector.")
    organ_vector = phantom_vector(:,:,min_slice:max_slice)

    write(10,*)section_prompt,'Detecting Organ Section'
    write(10,*)'Min and Max z values:'
    write(10,'(A,I4,2X,I4)') line_prompt,min_z, max_z
    write(10,*)'Section size: (slices, cm)'
    write(10,'(A,I4,2X,F7.4)') line_prompt,max_z-min_z+1, (max_z-min_z+1)*z_res
    write(10,*)'Section considered for redimension:'
    write(10,'(A,I4,2X,I4)') line_prompt,min_slice,max_slice
    write(10,*)'  '

    write(*,*)' Listing the initial coordinates on all 3 axis.'
    write(10,*)section_prompt,'Listing the initial coordinates on all 3 axis'

          ! --- SECTION: Get Initial Dimensions ---
    allocate( x_vector( (1+n_rows*(min_slice-1)):(n_rows+n_rows*(max_slice-1)), 1:4), stat=ios)
    if (ios /= 0) call HandleError("Failed to allocate memory for x_vector.")
    allocate( y_vector( (1+n_cols*(min_slice-1)):(n_cols+n_cols*(max_slice-1)), 1:4), stat=ios)
    if (ios /= 0) call HandleError("Failed to allocate memory for y_vector.")
    allocate( z_vector( (1+(n_cols*(1-1))):(n_cols+(n_cols*(n_rows-1))), 1:4), stat=ios)
    if (ios /= 0) call HandleError("Failed to allocate memory for z_vector.")

    call CalculateLungDimensions(min_slice, max_slice, organ_id_in)
    write(*,*)'  '
    write(10,*)'  '

  end subroutine ReadPhantomDataAndGetDimensions


  ! ==============================================================================
  ! SUBROUTINE: ApplyModifications
  !
  ! PURPOSE:
  !   Reads a list of modification tasks from specified tasks file and applies
  !   them sequentially to the target organ, updating the 'organVector'. Each task
  !   involves adding or deleting a single voxel layer from the target organ along
  !   a specified axis and side.
  !
  ! ARGUMENTS:
  !   organ_id_in, other_lung_id_in, lung_blood_id_in (input, integer): IDs for organs.
  !   tasks_file_lines_in (input, integer): Number of tasks to process.
  !   x_res, y_res, z_res (input, real): Voxel resolutions for reporting.
  !
  ! ==============================================================================
  subroutine ApplyModifications(organ_id_in, other_lung_id_in, lung_blood_id_in, tasks_file_lines_in)
    use PhantomGlobalVars
    implicit none

    ! --- Arguments ---
    integer, intent(in) :: organ_id_in, other_lung_id_in, lung_blood_id_in
    integer, intent(in) :: tasks_file_lines_in

    ! --- Local Variables ---
    integer :: i, ios, min_slice, max_slice
    character(len=1) :: choiceAxis, choiceAddDel, choicePosNeg
    integer :: z, y, x
    integer :: xaxis_z, xaxis_y, xaxis_min_x, xaxis_max_x, counter_x
    integer :: yaxis_z, yaxis_x, yaxis_min_y, yaxis_max_y, counter_y
    integer :: zaxis_y, zaxis_x, zaxis_min_z, zaxis_max_z, counter_z
    logical :: ex
    character(len=256) :: error_msg

    tasks_file_lines = tasks_file_lines_in
    min_slice = lbound(organ_vector, 3)
    max_slice = ubound(organ_vector, 3)

    inquire(file=trim(tasks_file), exist=ex)
    if (.not. ex) call HandleError("Tasks file not found: " // trim(tasks_file))
    open(unit=12, file=trim(tasks_file), status='old', iostat=ios)
    if (ios /= 0) call HandleError("Error opening tasks file: " // trim(tasks_file))

    ! Patch to skip BOM and leading tabs before reading the file content
    call SkipBOMAndTabs(12)

    do i=1,4! Read and discard header line
      read (12,*,iostat=ios)
      if (ios /= 0) call HandleError("Error reading header from tasks file.")
    end do

    ! Write to report file
    write(10,*)fwd_line
    write(10,*)section_prompt,'SECTION: Organ Redimension Structure.'
    write(10,*)'  '
    write(10,*)section_prompt,'Load Tasks File'
    write(10,*)'Opening tasks file:'
    write(10,*)line_prompt,trim(tasks_file)
    write(10,*)'No of lines to read:'
    write(10,'(A,I3)')line_prompt,tasks_file_lines
    write(10,*)'Starting organ redimension...'
    write(10,*)'  '

    ! --- Main Modification Loop ---
    do i = 1, tasks_file_lines
      read (12,*,iostat=ios) choiceAxis, choiceAddDel, choicePosNeg
      if (ios /= 0) then
        write(error_msg,'(A,I1,A)')"Error reading task line ",i," from tasks file."
        call HandleError(trim(error_msg))
      endif


      write(10,*)section_prompt,'Reading line from tasks file. Task no:',i
      if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'n') ) then
              write(10,*)' Adding a slice... On the negative side of the ',choiceAxis, ' axis...'
              write(*,*)' Adding a slice... On the negative side of the ',choiceAxis, ' axis...'
      else if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'p') ) then
              write(10,*)' Adding a slice... On the positive side of the ',choiceAxis, ' axis...'
              write(*,*)' Adding a slice... On the positive side of the ',choiceAxis, ' axis...'
      else if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'n') ) then
              write(10,*)' Deleting a slice... On the negative side of the ',choiceAxis, ' axis...'
              write(*,*)' Deleting a slice... On the negative side of the ',choiceAxis, ' axis...'
      else if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'p') ) then
              write(10,*)' Deleting a slice... On the positive side of the ',choiceAxis, ' axis...'
              write(*,*)' Deleting a slice... On the positive side of the ',choiceAxis, ' axis...'
      end if

      ! --- Resize organ in the X direction ---
      if (choiceAxis.eq.'x') then
        do z=min_slice,max_slice
          X_LOOP: do y=1,n_rows
                      counter_x=y+n_rows*(z-1)
                      xaxis_z = x_vector(counter_x,1)
                      xaxis_y = x_vector(counter_x,2)
                      xaxis_min_x = x_vector(counter_x,3)
                      xaxis_max_x = x_vector(counter_x,4)

                      if (xaxis_min_x > n_cols .or. xaxis_max_x == 0) cycle X_LOOP ! Skip if no organ in this row

                      if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'n') ) then
                        if( (organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)==other_lung_id_in) .or.&
                            ((organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)<=max_bone_id) .and.&
                            (organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)>=min_bone_id)) ) then
                            cycle X_LOOP ! Skip this modification
                        end if
                        organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)=organ_id_in
                      end if

                      if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'p') ) then
                        if( (organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)==other_lung_id_in) .or.&
                            ((organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)<=max_bone_id) .and.&
                            (organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)>=min_bone_id)) ) then
                            cycle X_LOOP ! Skip this modification
                        end if
                        organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)=organ_id_in
                      end if

                      if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'n') ) then
                        if( (organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)==other_lung_id_in) .or.&
                            ((organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)<=max_bone_id) .and.&
                            (organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)>=min_bone_id) ) .or.&
                            (organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)==0) ) then
                            organ_vector(xaxis_y,xaxis_min_x,xaxis_z)=heart_blood_id
                        else if(organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)==lung_blood_id_in) then
                            organ_vector(xaxis_y,xaxis_min_x,xaxis_z)=heart_blood_id
                            organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)=heart_blood_id
                        else
                            organ_vector(xaxis_y,xaxis_min_x,xaxis_z) = organ_vector(xaxis_y,xaxis_min_x-1,xaxis_z)
                        end if
                      end if

                      if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'p') ) then
                        if( (organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)==other_lung_id_in) .or.&
                            ((organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)<=max_bone_id) .and.&
                            (organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)>=min_bone_id) ) .or.&
                            (organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)==0) ) then
                            organ_vector(xaxis_y,xaxis_max_x,xaxis_z)=heart_blood_id
                        else if((organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)==lung_blood_id_in)) then
                            organ_vector(xaxis_y,xaxis_max_x,xaxis_z)=heart_blood_id
                            organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)=heart_blood_id
                        else
                            organ_vector(xaxis_y,xaxis_max_x,xaxis_z) = organ_vector(xaxis_y,xaxis_max_x+1,xaxis_z)
                        end if
                      end if
                  end do X_LOOP
        end do
      end if

      ! --- Resize organ in the Y direction ---
      if (choiceAxis.eq.'y') then
        do z=min_slice,max_slice
          Y_LOOP: do x=1,n_cols
                        counter_y = x+n_cols*(z-1)
                        yaxis_z = y_vector(counter_y,1)
                        yaxis_x = y_vector(counter_y,2)
                        yaxis_min_y = y_vector(counter_y,3)
                        yaxis_max_y = y_vector(counter_y,4)
                        if (yaxis_min_y > n_rows .or. yaxis_max_y == 0) cycle Y_LOOP ! Skip this modification

                        if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'n') ) then
                          if( (organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)==other_lung_id_in) .or.&
                              ((organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)<=max_bone_id) .and.&
                              (organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)>=min_bone_id)) ) then
                              cycle Y_LOOP ! Skip this modification
                          end if
                          organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)=organ_id_in
                        end if

                        if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'p') ) then
                          if( (organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)==other_lung_id_in) .or.&
                              ((organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)<=max_bone_id) .and.&
                              (organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)>=min_bone_id)) ) then
                              cycle Y_LOOP ! Skip this modification
                          end if
                          organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)=organ_id_in
                        end if

                        if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'n') ) then
                          if( (organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)==other_lung_id_in) .or.&
                              ((organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)<=max_bone_id) .and.&
                              (organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)>=min_bone_id)) .or.&
                              (organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)==0) ) then
                              organ_vector(yaxis_min_y,yaxis_x,yaxis_z)=heart_blood_id
                          else if (organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)==lung_blood_id_in) then
                              organ_vector(yaxis_min_y,yaxis_x,yaxis_z)=heart_blood_id
                              organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)=heart_blood_id
                          else
                              organ_vector(yaxis_min_y,yaxis_x,yaxis_z) = organ_vector(yaxis_min_y-1,yaxis_x,yaxis_z)
                          end if
                        end if

                        if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'p') ) then
                          if( (organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)==other_lung_id_in) .or.&
                              ((organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)<=max_bone_id) .and.&
                              (organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)>=min_bone_id)) .or.&
                              (organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)==0) ) then
                              organ_vector(yaxis_max_y,yaxis_x,yaxis_z)=heart_blood_id
                          else if (organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)==lung_blood_id_in) then
                              organ_vector(yaxis_max_y,yaxis_x,yaxis_z) = heart_blood_id
                              organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)=heart_blood_id
                          else
                              organ_vector(yaxis_max_y,yaxis_x,yaxis_z) = organ_vector(yaxis_max_y+1,yaxis_x,yaxis_z)
                          end if
                        end if
                  end do Y_LOOP
        end do
      end if

      ! --- Resize organ in the Z direction ---
      if (choiceAxis.eq.'z') then
        do y=1,n_rows
          Z_LOOP: do x=1,n_cols
                        counter_z = x+n_cols*(y-1)
                        zaxis_y = z_vector(counter_z,1)
                        zaxis_x = z_vector(counter_z,2)
                        zaxis_min_z = z_vector(counter_z,3)
                        zaxis_max_z = z_vector(counter_z,4)
                        if (zaxis_min_z > n_slices .or. zaxis_max_z == 0) cycle Z_LOOP ! Skip this modification

                        if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'n') ) then
                          if( (organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)==other_lung_id_in) .or.&
                              ((organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)<=max_bone_id) .and.&
                              (organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)>=min_bone_id)) ) then
                              cycle Z_LOOP ! Skip this modification
                          end if
                          organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)=organ_id_in
                        end if

                        if ( (choiceAddDel.eq.'a').and.(choicePosNeg.eq.'p') ) then
                          if( (organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)==other_lung_id_in) .or.&
                              ((organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)<=max_bone_id) .and.&
                              (organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)>=min_bone_id)) ) then
                              cycle Z_LOOP ! Skip this modification
                          end if
                          organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)=organ_id_in
                        end if

                        if ( (choiceAddDel.eq.'d').and.(choicePosNeg.eq.'n') ) then
                          if( (organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)==other_lung_id_in) .or.&
                              ((organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)<=max_bone_id) .and.&
                              (organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)>=min_bone_id)) .or.&
                              (organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)==0) ) then
                              organ_vector(zaxis_y,zaxis_x,zaxis_min_z)=heart_blood_id
                          else if (organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)==lung_blood_id_in) then
                              organ_vector(zaxis_y,zaxis_x,zaxis_min_z)=heart_blood_id
                              organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)=heart_blood_id
                          else
                              organ_vector(zaxis_y,zaxis_x,zaxis_min_z) = organ_vector(zaxis_y,zaxis_x,zaxis_min_z-1)
                          end if
                        end if

                        if ( (choiceAddDel.eq.'d' ).and.(choicePosNeg.eq.'p') ) then
                          if( (organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)==other_lung_id_in) .or.&
                              ((organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)<=max_bone_id) .and.&
                              (organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)>=min_bone_id)) .or.&
                              (organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)==0) ) then
                              organ_vector(zaxis_y,zaxis_x,zaxis_max_z)=heart_blood_id
                          else if (organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)==lung_blood_id_in) then
                              organ_vector(zaxis_y,zaxis_x,zaxis_max_z)=heart_blood_id
                              organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)=heart_blood_id
                          else
                              organ_vector(zaxis_y,zaxis_x,zaxis_max_z) = organ_vector(zaxis_y,zaxis_x,zaxis_max_z+1)
                          end if
                        end if
                  end do Z_LOOP
        end do
      end if

      ! --- Recalculate dimensions after each modification ---
      write(10,*)' Current dimensions in all 3 axis after this modification:'
      call CalculateLungDimensions(min_slice, max_slice, organ_id_in)

      ! Write empty line(s) to report file
      write(10,*)' '
      write(10,*)' '
      write(*,*)' '
      write(*,*)' '

    end do
    close(12)

  end subroutine ApplyModifications

  ! ==============================================================================
  ! SUBROUTINE: WriteOutputFiles
  !
  ! PURPOSE:
  !   Output file generation and final operations. Cleanup process is performed to
  !   check for any lung blood voxels outside lung. Then, a final count the number 
  !   of voxels for key organs is performed and the complete, modified phantom is
  !   written to output files: a simulation file ('.vox') and a visualization file
  !   (type 'ct-den-mat-XY.dat').
  !
  ! ARGUMENTS:
  !   organ_id_in (input, integer): ID of the target organ.
  !   lung_blood_id_in (input, integer): ID for the target lung's blood.
  !   organ_name_in    (input, character): The name of the organ for reporting.
  !   x_res, y_res, z_res (input, real): Voxel resolutions for the .vox file header.
  !   n_voxels_in      (input, integer): Total number of voxels.
  !
  ! ==============================================================================
  subroutine WriteOutputFiles(organ_id_in, lung_blood_id_in, n_voxels_in)
    use PhantomGlobalVars
    implicit none

    ! --- Arguments ---
    integer, intent(in)           :: organ_id_in
    integer, intent(in)           :: lung_blood_id_in
    integer, intent(in)           :: n_voxels_in

    ! --- Local Variables ---
    integer :: blood_out_lungs_counter
    integer :: x, y, z, i, counter, min_slice, max_slice, ios
    integer :: abs_min_x, abs_max_x, abs_min_y, abs_max_y, abs_min_z, abs_max_z
    character(len=1) :: choice_end_script
    character(len=256) :: error_msg
    integer, dimension(:), allocatable :: organlist_org, organlist_mat
    real, dimension(:), allocatable    :: organlist_den
    real, allocatable, dimension(:,:) :: mat_den_vector
    integer :: n_right_voxels, n_right_blood, n_right_tissue
    integer :: n_left_voxels, n_left_blood, n_left_tissue
    integer :: n_lung_voxels, n_esophagus_voxels, n_stomach_voxels
    integer :: organ, material
    real    :: density

    min_slice = lbound(organ_vector, 3)
    max_slice = ubound(organ_vector, 3)

    ! --- Final Dimension Analysis ---
    ! Recalculate dimensions after all modifications have been applied.
    write(10,*)
    write(*,*)' Recalculate dimensions after all modifications have been applied.'
    write(10,*)fwd_line
    write(10,*)section_prompt,'SECTION: Final Organ Dimensions Analysis.'
    write(10,*)' Recalculate dimensions after all modifications have been applied.'
    call CalculateLungDimensions(min_slice, max_slice, organ_id_in)

    ! --- SECTION: Cleanup Lung Blood ---
    blood_out_lungs_counter = 0
    write(10,*)fwd_line
    write(10,*)section_prompt,'SECTION: Lung Blood Redimension.'
    write(10,*)'  '
    write(10,*)section_prompt,'Check if there are any Lung Blood Voxels outside the Lung'
    do z=1,n_slices
      do y=1,n_rows
        do x=1,n_cols
          if ( (x<abs_min_x .or. x>abs_max_x .or. y<abs_min_y .or. y>abs_max_y .or. &
               z<abs_min_z .or. z>abs_max_z) .and. phantom_vector(x,y,z)==lung_blood_id_in ) then
            phantom_vector(x,y,z)=heart_blood_id
            blood_out_lungs_counter=blood_out_lungs_counter+1
          end if
        end do
      end do
    end do
    write(*,*)' Analysis completed.'
    write(10,*)' Number of Lung Blood Voxels found outside the Lung:'
    write(10,'(A,I6)')line_prompt, blood_out_lungs_counter
    write(10,*)' These Lung Blood Voxels were eliminated and replaced by the Organ ID ', heart_blood_id
    write(*,*)' '
    write(10,*)' '

    ! --- SECTION: Write new phantom file ---
    write(*,*)' Do you wish to proceed to writing new phantom in file?'
    read(*,*) choice_end_script
    if ((choice_end_script.eq.'n') .or. (choice_end_script.eq.'N')) stop 'The script will stop at your request.'

    write(10,*)' '
    write(10,*)fwd_line
    write(10,*)section_prompt,'SECTION: Writing new phantom in file.'
    write(10,*)'  '
    phantom_vector(:,:,min_slice:max_slice)=organ_vector
    write(*,*)' What is the name of the file to write the new phantom?'
    read (*,*) phantom_out_file
    write(10,*)' Writing the new file:'
    write(10,*)line_prompt,trim(phantom_out_file)
    open(unit=13, file=trim(phantom_out_file), status='replace', iostat=ios)
    if (ios /= 0) call HandleError("Could not create output phantom file: " // trim(phantom_out_file))

    ! Load organ list to map organ IDs to materials and densities
    write(*,*)' Reading organlist file:', trim(organlist_file)
    write(10,*)' Reading organlist file: ', trim(organlist_file)
    allocate(organlist_org(organlist_file_lines), organlist_mat(organlist_file_lines), organlist_den(organlist_file_lines))

    ! Inquire if organlist file exists. Try to open.
    inquire(file=trim(organlist_file), exist=ex)
    if (.not. ex) call HandleError("Organlist file not found: " // trim(organlist_file))
    open(unit=14, file=trim(organlist_file), status='old', iostat=ios)
    if (ios /= 0) call HandleError("Error opening organlist file: " // trim(organlist_file))

    ! Patch to skip BOM and leading tabs before reading the file content
    call SkipBOMAndTabs(14)
 
    read (14,*,iostat=ios) ! Read and discard header line from organlist file
    if (ios /= 0) call HandleError("Error reading header from organlist file.")

    ! Read organlist file
    do i=1,organlist_file_lines
        read (14,*,iostat=ios) organ, material, density
        if (ios > 0) then
            write(error_msg,'(A,I3,A)')"Error reading line ",i," from organlist file."
            call HandleError(trim(error_msg)) ! Call error message from error reading line
        else if (ios < 0) then
            write(error_msg,'(A,I3,A)')"Unexpected end of file on line ",i," from organlist file."
            call HandleError(trim(error_msg)) ! Call error message from unexpected end of file
        end if
        organlist_org(i) = organ
        organlist_mat(i) = material
        organlist_den(i) = density
    end do
    close(14)
    write(10,*) ' Finished analysing organlist file.'
    ! Write new ct-den-matXY phantom file
    write(*,*) ' Writing the new file: ',trim(phantom_out_file)
    allocate(mat_den_vector(n_voxels_in,1:2))
    do z=1,n_slices
      do y=1,n_rows
        do x=1,n_cols
            counter = x+n_cols*((y-1)+n_rows*(z-1))
            do i=1, organlist_file_lines
                if (phantom_vector(x,y,z) == organlist_org(i)) then
                    material = organlist_mat(i)
                    density = organlist_den(i)
                    organ = organlist_org(i)
                    exit
                end if
            end do
            write(13,'(1x,i4,1x,i4,1x,i4,1x,f5.3,1x,i3,3x,i3)')x,y,z, &
                                                        & density,material,organ
            mat_den_vector(counter,1)=material
            mat_den_vector(counter,2)=density
        end do
        write(13,*)'  '
      end do
      write(13,*)'  '
      call ReadProgress(z)
    end do
    close(13)
    deallocate(organlist_org, organlist_mat, organlist_den)

    ! --- SECTION: Count voxels in final organs ---
    write(*,*)' Checking number of voxels in the phantom...'
    n_right_voxels=0; n_right_blood=0; n_right_tissue=0
    n_left_voxels=0; n_left_blood=0; n_left_tissue=0
    n_lung_voxels=0; n_esophagus_voxels=0; n_stomach_voxels=0
    do z=1,n_slices
      do y=1,n_rows
        do x=1,n_cols
            if (phantom_vector(x,y,z) == right_lung_id) then
                n_right_voxels=n_right_voxels+1; n_right_tissue=n_right_tissue+1
            else if (phantom_vector(x,y,z) == right_lung_blood_id) then
                n_right_voxels=n_right_voxels+1; n_right_blood=n_right_blood+1
            else if (phantom_vector(x,y,z) == left_lung_id) then
                n_left_voxels=n_left_voxels+1; n_left_tissue=n_left_tissue+1
            else if (phantom_vector(x,y,z) == left_lung_blood_id) then
                n_left_voxels=n_left_voxels+1; n_left_blood=n_left_blood+1
            else if (phantom_vector(x,y,z) == esophagous_id) then
                n_esophagus_voxels=n_esophagus_voxels+1
            else if (phantom_vector(x,y,z) == stomach_id) then
                n_stomach_voxels=n_stomach_voxels+1
            end if
        end do
      end do
    end do
    n_lung_voxels = n_right_voxels + n_left_voxels
    write(10,*)'Number of voxels in the final phantom:'
    write(10,*)line_prompt,'Right Lung       :', n_right_tissue
    write(10,*)line_prompt,'Left Lung.       :', n_left_tissue
    write(10,*)line_prompt,'Total Lung.      :', n_lung_voxels
    write(10,*)line_prompt,'Right Lung Blood :', n_right_blood
    write(10,*)line_prompt,'Left Lung Blood  :', n_left_blood
    write(10,*)line_prompt,'Stomach          :', n_stomach_voxels
    write(10,*)line_prompt,'Esophagus        :', n_esophagus_voxels

    ! --- SECTION: Write .vox file ---
    write(*,*)' Do you wish to proceed to writing the .vox file?'
    read(*,*) choice_end_script
    if ((choice_end_script.eq.'n') .or. (choice_end_script.eq.'N')) stop 'The script will stop at your request.'

    open(unit=16, file='phantom.vox', status='replace', iostat=ios)
    if (ios /= 0) call HandleError("Could not create output file 'phantom.vox'.")
    write(*,*)' Writing the .vox file ...'
    write(16,*)'[SECTION VOXELS HEADER v.2008-04-13]'
    write(16,'(3(2x,I4))') n_cols,n_rows,n_slices
    write(16,'(3(2X,F7.4))') x_res,y_res,z_res
    write(16,*)' 1'
    write(16,*)' 2'
    write(16,*)' 0'
    write(16,*)'[END OF VXH SECTION]'
    do z=1,n_slices
      do y=1,n_rows
        do x=1,n_cols
            counter = x+n_cols*((y-1)+n_rows*(z-1))
            write(16,'(1X,I3,2X,F5.3)') int(mat_den_vector(counter,1)), mat_den_vector(counter,2)
        end do
      end do
      call ReadProgress(z)
    end do
    close(16)
    deallocate(mat_den_vector)
    write(*,*)' Finished writing .vox file.'

  end subroutine WriteOutputFiles


end program ReDimLungs
