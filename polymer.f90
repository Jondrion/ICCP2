module polymer
  use plplot

  implicit none
  
  real(8), parameter :: PI = 4 * atan(1d0)

  private

  public polymerType
  type polymerType
    private
   
    integer, private :: Length
    integer, private :: NumberAngles
    real(8), private :: PolWeight
    real(8), allocatable, private :: Position(:,:)

   contains
    private
     
    procedure, public :: init, destroy
    procedure, private :: create, calc_Angles, calc_Energy, calc_Weights
    procedure, private :: choose_Angle, cumsum
    procedure, public :: get_Position, get_Length
    procedure, public :: plot
    procedure, private :: plot_polymer
     
   end type

contains
  
  ! -- initialize the polymer: allocate position array of certain length and create the first two beads
  subroutine init(this,Length)
  
    class(polymerType) :: this
    integer, intent(in) :: Length
    
    this%Length=Length
    allocate (this%Position(2,this%Length))

    ! -- create inital two beads
    this%Position(1,1)=0
    this%Position(2,1)=0
    this%Position(1,2)=0
    this%Position(2,2)=1

    ! -- set polymer weight to 1
    this%PolWeight=1

    this%NumberAngles=10

    call this%create(this%PolWeight, 3)

  end subroutine


  ! -- creating the polymer by adding enough beads until length is reached, this routine is recursive
  recursive subroutine create(this, PolWeight, Number)

    class(polymerType) :: this
    real(8) :: Angle(this%NumberAngles)
    real(8) :: Energy(this%NumberAngles)
    integer, intent(in) :: Number
    real(8), intent(in) :: PolWeight
    real(8) :: Weights(this%NumberAngles)
    real(8) :: Norm_Weights(this%NumberAngles)
    real(8) :: SumWeights

    ! -- Calculate weights w_j^L and their product W^L
    call this%calc_Angles(Angle,this%NumberAngles)
    print *, "Angles"
    print *, Angle(:)/PI
    call this%calc_Energy(Angle,Number,Energy)
    call this%calc_Weights(Energy,Weights)

    SumWeights = SUM(Weights)

    Norm_Weights = Weights/SumWeights
    print *, "Weights"
    print *, Norm_Weights(:)

    ! -- Add bead number L
    this%Position(1,Number)=0
    this%Position(2,Number)=Number-1

    ! -- Update PolWeight

    ! -- recursive part
    if ( Number < this%Length ) then
      call this%create(PolWeight, Number+1)
    end if

  end subroutine


  subroutine destroy(this)
    class(polymerType) :: this

    deallocate(this%Position)
  end subroutine


  subroutine calc_Angles(this, Angle, N)

    class(polymerType) :: this
    integer, intent(in) :: N
    integer :: i
    real(8), intent(out) :: Angle(N)
    real(8) :: Interval
    real(8) :: Offset

    Interval = 2 * PI / N
    call RANDOM_NUMBER(Offset)
    Offset = Offset * Interval
    Angle(:) = [(Offset + ((i-1)*Interval), i=1,N)]

  end subroutine


  subroutine calc_Energy(this, Angle, N, E)
    
    class(polymerType) :: this
    real(8), intent(in) :: Angle(this%NumberAngles)
    integer, intent(in) :: N
    real(8) :: Ri(2)
    real(8) :: Rsqi
    real(8) :: rm2,rm6,rm12
    real(8) :: new_pos(2,this%NumberAngles)
    real(8), intent(out) :: E(this%NumberAngles)

    integer :: i,j

    E=0

    ! -- calculate position given the angle
    new_pos(1,:) = this%Position(1,N-1) + SIN(Angle(:))
    new_pos(2,:) = this%Position(2,N-1) + COS(Angle(:))

    ! -- loop over all angles
    do j = 1, this%NumberAngles, 1
      ! -- loop over all existing elements
      do i = 1, N-1, 1
        Ri = new_pos(:,j) - this%Position(:,i)
        Rsqi = dot_product(Ri,Ri)
        rm2 = 1.d0/Rsqi
        rm6 = rm2**3
        rm12 = rm6**2
        E(j) = E(j) + 4.d0 * ( rm12 - rm6 )
      end do
    end do

  end subroutine


  subroutine calc_Weights(this, Energy, Weights)
    
    class(polymerType) :: this
    real(8), intent(out) :: Weights(this%NumberAngles)
    real(8), intent(in) :: Energy(this%NumberAngles)

    Weights(:) = EXP(-Energy(:))

  end subroutine


  subroutine choose_Angle(this, Angles, Norm_Weights, output_Angle)

    class(polymerType) :: this
    real(8), intent(inout) :: Angles(this%NumberAngles), Norm_Weights(this%NumberAngles)
    real(8), intent(out) :: output_Angle
    real(8) :: Roulette
    integer :: angle_number

    call RANDOM_NUMBER(Roulette)
    call this%cumsum(Norm_Weights)

    where (Norm_Weights<Roulette) Norm_Weights=1._8
    where (Norm_Weights>=Roulette) Norm_Weights=0._8

    angle_number = sum(nint(Norm_Weights))+1

    output_Angle = Angles(angle_number)

  end subroutine


  subroutine cumsum(this, Array)
    class(polymerType) :: this
    real(8), intent(inout) :: Array(this%NumberAngles)
    integer :: i
    
    do i=1, this%NumberAngles-1
      Array(i+1)=Array(i+1)+Array(i)
    end do

  end subroutine


  subroutine get_Position(this, Position)

    class(polymerType) :: this
    real(8), allocatable, intent(out) :: Position(:,:)

    Position=this%Position

  end subroutine

  subroutine get_Length(this, Length)

    class(polymerType) :: this
    real(8), intent(out) :: Length

    Length=this%Length

  end subroutine


  subroutine plot(this)

    class(polymerType) :: this

    call plsdev("xcairo")

    ! gnuplot color scheme
    call plscol0(0, 255, 255, 255)  ! white
    call plscol0(1, 255, 0, 0)      ! red
    call plscol0(2, 0, 255, 0)      ! green
    call plscol0(3, 0, 0, 255)      ! blue
    call plscol0(4, 255, 0, 255)    ! magenta
    call plscol0(5, 0, 255, 255)    ! cyan
    call plscol0(6, 255, 255, 0)    ! yellow
    call plscol0(7, 0, 0, 0)        ! black
    call plscol0(8, 255, 76, 0)     ! orange
    call plscol0(9, 128, 128, 128)  ! gray

    call plinit()

    call this%plot_polymer

    call plend()

  end subroutine
  

  subroutine plot_polymer(this)

    class(polymerType) :: this
    real(8) :: min, max,border
  
    min = minval(this%Position)
    max = maxval(this%Position)
  
    border=(max-min)/20
    min = min - border
    max = max + border
    
    call plcol0(7)
    call plenv(min, max, min, max, 0, 0)
    call pllab("x", "y", "Polymer")
  
  
    call plcol0(1)
    call plline(this%Position(1,:),this%Position(2,:))
    
    call plcol0(2)
    call plpoin(this%Position(1,:),this%Position(2,:),4)

  end subroutine


end module polymer