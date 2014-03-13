module Protein
  
  implicit none
  private

  public Prot
  type Prot
     private
   
     integer, public :: N
     real(8), allocatable, public :: x(2,:)

   contains
     public
     
     procedure :: create, destroy
     procedure :: clone
     procedure :: get_N, get_Ener, get_endtoend
     
   end type

contains
  
  subroutine create(this,N)
  
    class(Prot) :: this
    integer, intent(in) :: N
    
    this%N=N
  end subroutine
