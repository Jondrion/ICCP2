program Poly
  use Protein

  implicit none
  real(8) :: density,temperature
  real(8), allocatable :: Position(:,:)
  type(proteinType) :: polymer

  call parameters(density,temperature)
  print *, "Density=",density, "Temperature=",temperature

  call polymer%init(10)
  call polymer%get_Position(Position)
  print *, "Position"
  print *, Position(:,:)
  ! call plot! 


end program
