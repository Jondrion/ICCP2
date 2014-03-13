program md

use gas
use plot

implicit none

integer, parameter :: LinCell = 6, CorrPts = 250
type(gasType) :: LJSys

real(8), allocatable :: LastDiffPos(:,:,:), AvDist(:), NewDiff(:,:)

integer :: InitStep, SimulStep, rescaleStep, CorrUpdates, &
            DiffStatNum, DiffStepMax, DisplayStep, UpdateStep

real(8) :: density, temperature, TimeStep, PI, BoxSize
integer :: CorrArray(0:CorrPts)

  call init_simul(DisplayStep, RescaleStep, UpdateStep)
  call simulate(.false., InitStep)
  call simulate(.true., SimulStep)
  call analyze_diff
#ifdef Plot
  call end_plot
#endif

contains
  
  subroutine init_simul(DisplayStep, RescaleStep, UpdateStep)
    integer, intent(out) :: DisplayStep, RescaleStep, UpdateStep 
    
    open (12, file="md.params")
    read (12, *) density
    read (12, *) temperature
    read (12, *) InitStep
    read (12, *) SimulStep
    read (12, *) TimeStep
    read (12, *) DiffStatNum
    read (12, *) DisplayStep, RescaleStep, UpdateStep 
    close(12)

    CorrArray = 0
    
    PI = 4._8*atan(1._8)
    call LJSys%create(LinCell, Density, Temperature)
    call LJSys%set_temp(Temperature)
    call LJSys%get_boxsize(BoxSize)
    
    DiffStepMax = int(SimulStep/(DiffStatNum+1))
    ALLOCATE (LastDiffPos(1:DiffStepMax, 3, LJSys%n))
    ALLOCATE(AvDist(DiffStepMax))
    ALLOCATE (NewDiff(3,LJSys%n))
    AvDist = 0._8

#ifdef Plot
    call init_plot(BoxSize)
#endif
  end subroutine init_simul

!end program md

  subroutine simulate(produce, StepNum)
  logical, intent(in) :: produce
  integer, intent(in) :: StepNum
  integer :: Step
  real(8) :: Potential, Kinetic
  
  call LJSys%make_PairList
  do Step = 1, StepNum
    call LJSys%move_part(TimeStep)
    if (produce) then
      if (step == 1) then
        NewDiff = 0._8
      else
        newdiff = NewDiff + LJSys%mom*TimeStep
        write (31,*) sum(newdiff*newdiff)
      end if
      if (mod(step,updatestep)==0) call update_corr(LJSys%Pos,CorrArray)
      call update_diff_pos(step)
      call LJSys%get_kin(Kinetic)
      write (17,*) Kinetic
      call flush(31)
    end if
    if (mod(step,DisplayStep)==0) then
       call LJSys%get_pot(Potential)
       call LJSys%get_kin(Kinetic)
       print *, Potential, Kinetic, Potential + Kinetic
#ifdef Plot       
       call plw3d(1d0, 1d0, 1d0, 0._8, BoxSize, 0._8, BoxSize, 0._8, BoxSize, 45d0, 45*sin(0.001_8*Step))
       CALL plot_points(modulo(LJSys%Pos,BoxSize), BoxSize)
#endif
    end if
    if (.not.(produce).and.(mod(step,RescaleStep)==0)) then
      call LJSys%set_temp(Temperature)
    end if
    if (mod(step, 100) == 0) call LJSys%make_PairList
  end do
  if (produce) call write_corr
  end subroutine simulate 

  subroutine update_corr(Pos, CorrArray)
  integer :: i, j, CorrIndex
  real(8) :: relpos(3), r2
  real(8), intent(in) :: Pos(:,:)
  integer, intent(inout) :: CorrArray(:)

  do i = 1, LJSys%n-1
    do j = i+1, LJSys%n
      relpos = pos(:,i) - pos(:,j)
      relpos = relpos - nint(relpos/BoxSize)*BoxSize
      r2 = sum(relpos*relpos)
      CorrIndex = int(2._8*sqrt(r2)/BoxSize*CorrPts)
      if (CorrIndex<=CorrPts) CorrArray(CorrIndex) = CorrArray(CorrIndex)+1
    end do
  end do
  CorrUpdates = CorrUpdates + 1
  end subroutine update_corr
  
  
  subroutine write_corr
  integer :: i
  real(8) :: DeltaR, Radius
   
  DeltaR = BoxSize/(2._8*CorrPts)
  do i = 0, CorrPts
    Radius = (i+0.5_8)*DeltaR
    WRITE(11,'(F12.5)') 2._8*CorrArray(i)/(4*PI*Radius**2*DeltaR)/Density/CorrUpdates/LJSys%n
  end do
  end subroutine write_corr



  subroutine update_diff_pos(step)
  integer, intent(in) :: step
  real(8) :: RelPos(3,LJSys%n)
  integer :: k,l

  do k=1, DiffStepMax
    if (mod(step, k).eq.0) then 
      l = step/k
      if (l<=DiffStatNum+1) then
        if (l>1) then
          RelPos = LJSys%Pos - LastDiffPos(k,:,:)
          RelPos = RelPos - nint(RelPos/BoxSize)*BoxSize
          AvDist(k) = AvDist(k) + sum(RelPos*RelPos)/LJSys%n
        end if
        LastDiffPos(k,:,:) = LJSys%Pos
      end if
    end if
  end do
  end subroutine update_diff_pos

  

  subroutine analyze_diff
  integer :: m
  
  do m = 1, DiffStepMax
    AvDist(m) = AvDist(m)/DiffStatNum
    write (13,*) m, AvDist(m)
  end do
  end subroutine analyze_diff 
      
      


end program md
