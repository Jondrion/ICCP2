module polymer
	implicit none
  private

  public polymerType
  type polymerType
    private
   
    integer, private :: Length
    real(8), private :: PolWeight
    real(8), allocatable, private :: Position(:,:)

   contains
    private
     
    procedure, public :: init, destroy
    procedure, private :: create
    procedure, public :: get_Position, get_Length
     
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

    call this%create(this%PolWeight, 3)

  end subroutine

  ! -- creating the polymer by adding enough beads until length is reached, this routine is recursive
  subroutine create(this, PolWeight, Number)

    class(polymerType) :: this
    integer, intent(in) :: Number
    real(8), intent(in) :: PolWeight

    ! -- Calculate weights w_j^L and their product W^L

    ! -- Add bead number L
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


end module polymer