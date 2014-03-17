module Protein
  
  implicit none
  private

  public proteinType
  type proteinType
     private
   
     integer, public :: Length
     real(8), allocatable, public :: Position(2,:)
     real(8), public :: Energy
     real(8), public :: EtoEdis

   contains
     public
     
     procedure :: init, destroy
     procedure :: create
     procedure :: clone
     procedure :: get_Position, get_Length, get_Energy, get_EtoEdis
     
   end type

contains
  
  subroutine init(this,Length)
  
    class(proteinType) :: this
    integer, intent(in) :: Length
    
    this%Length=Length
  end subroutine

  subroutine destroy(this)
    class(proteinType) :: this

    deallocate(this%pos)
  end subroutine

  subroutine create(this, Number)

    class(proteinType) :: this
    integer, intent(in) :: Number

    do i=1, Number
      this%Position(1,i)=1+i
      this%Position(2,i)=1+i
    end do

  end subroutine

  subroutine clone(this,clone)

    class(proteinType) :: this
    class(proteinType), intent(out) :: clone

    clone=this

  end subroutine

  subroutine get_Position(this, Position)

    class(proteinType) :: this
    integer, intent(out) :: Position

    Position=this%Position

  end subroutine

  subroutine get_Length(this, Length)

    class(proteinType) :: this
    integer, intent(out) :: Length

    Length=this%Length

  end subroutine

  subroutine get_Energy(this, Energy)

    class(proteinType) :: this
    integer, intent(out) :: Energy

    Energy=this%Energy

  end subroutine

  subroutine get_EtoEdis(this, EtoEdis)

    class(proteinType) :: this
    integer, intent(out) :: EtoEdis

    EtoEdis=this%EtoEdis

  end subroutine




