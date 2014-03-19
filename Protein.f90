module Protein
  use plplot

  implicit none
  private



  public proteinType
  type proteinType
    private
   
    integer, public :: Length
    real(8), allocatable, private :: Position(:,:)
    real(8), public :: Energy
    real(8), public :: EtoEdis

   contains
    private
     
    procedure, public :: init, destroy
    procedure, private :: create
    procedure, public :: clone
    procedure, public :: get_Position, get_Length, get_Energy, get_EtoEdis
    procedure, public :: plot
    procedure, private :: plot_protein
     
   end type

contains
  
  subroutine init(this,Length)
  
    class(proteinType) :: this
    integer, intent(in) :: Length
    
    this%Length=Length
    allocate (this%Position(2,this%Length))

    call this%create(Length)

  end subroutine

  subroutine destroy(this)
    class(proteinType) :: this

    deallocate(this%Position)
  end subroutine

  subroutine create(this, Number)

    class(proteinType) :: this
    integer :: i
    integer, intent(in) :: Number


    do i=1, Number
      this%Position(1,i)=1+i
      this%Position(2,i)=10+i
    end do

  end subroutine

  subroutine clone(this,cloned)

    class(proteinType) :: this
    class(proteinType), intent(out), allocatable :: cloned

    allocate (cloned, source=this)

  end subroutine

  subroutine get_Position(this, Position)

    class(proteinType) :: this
    real(8), allocatable, intent(out) :: Position(:,:)

    Position=this%Position

  end subroutine

  subroutine get_Length(this, Length)

    class(proteinType) :: this
    real(8), intent(out) :: Length

    Length=this%Length

  end subroutine

  subroutine get_Energy(this, Energy)

    class(proteinType) :: this
    real(8), intent(out) :: Energy

    Energy=this%Energy

  end subroutine

  subroutine get_EtoEdis(this, EtoEdis)

    class(proteinType) :: this
    real(8), intent(out) :: EtoEdis

    EtoEdis=this%EtoEdis

  end subroutine

  subroutine plot(this)

    class(proteinType) :: this


    

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

    call this%plot_protein

    call plend()

  end subroutine
  

  subroutine plot_protein(this)

  class(proteinType) :: this
  real(8) :: min, max,border

  min = minval(this%Position)
  max = maxval(this%Position)

  border=(max-min)/20
  print *, min, max, border
  min = min - border
  max = max + border
 print *, min, max
  
  call plcol0(7)
  call plenv(min, max, min, max, 0, 0)
  call pllab("x", "y", "Protein")


  call plcol0(1)
  call plline(this%Position(1,:),this%Position(2,:))
  
  call plcol0(2)
  call plpoin(this%Position(1,:),this%Position(2,:),4)

   end subroutine

  


end module Protein