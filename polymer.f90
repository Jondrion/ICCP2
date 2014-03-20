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
    procedure, private :: create, calc_Angles
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
    this%Position(1,2)=1
    this%Position(2,2)=1

    ! -- set polymer weight to 1
    this%PolWeight=1

    this%NumberAngles=10

    call this%create(this%PolWeight, 3)

  end subroutine


  ! -- creating the polymer by adding enough beads until length is reached, this routine is recursive
  recursive subroutine create(this, PolWeight, Number)

	class(polymerType) :: this
    real(8) :: Angle(1,this%NumberAngles)
    integer, intent(in) :: Number
    real(8), intent(in) :: PolWeight

    ! -- Calculate weights w_j^L and their product W^L

    ! -- Add bead number L
    call this%calc_Angles(Angle,this%NumberAngles)

    ! 
    this%Position(1,Number)=Number-1
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
  	real(8), intent(out) :: Angle(1,N)
  	real(8) :: Interval
  	real(8) :: Offset

  	Interval = 2 * PI / N
  	call RANDOM_NUMBER(Offset)
  	Offset = Offset * Interval
  	Angle(1,:) = [(Offset + ((i-1)*Interval), i=1,N)]

  end subroutine


  subroutine calc_Energy(this, Angle, N, E)
  	
  	class(polymerType) :: this
  	real(8), intent(in) :: Angle(1,this%NumberAngles)
  	integer, intent(in) :: N
  	real(8) :: Ri(2)
  	real(8) :: Rsqi
  	real(8) :: rm2,rm6,rm12
  	real(8), intent(out) :: E

  	integer :: i

  	E=0

  	do i = 1, N-1, 1
  		Ri = this%Position(:,N) - this%Position(:,i)
  		Rsqi = dot_product(Ri,Ri)
  		rm2 = 1.d0/Rsqi
  		rm6 = rm2**3
  		rm12 = rm6**2
  		E = E + 4.d0 * ( rm12 - rm6 )
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
    call pllab("x", "y", "polymer")
  
  
    call plcol0(1)
    call plline(this%Position(1,:),this%Position(2,:))
    
    call plcol0(2)
    call plpoin(this%Position(1,:),this%Position(2,:),4)

  end subroutine



end module polymer