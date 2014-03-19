program Poly
  use Protein
  use polymer

  implicit none
  real(8) :: density,temperature
  real(8), allocatable :: Position(:,:)
  type(proteinType) :: prot
  type(polymerType) :: pol

  call parameters(density,temperature)
  print *, "Density=",density, "Temperature=",temperature

  call prot%init(10)
  call prot%get_Position(Position)
  print *, "Position"
  print *, Position(0,:)
  call polymer%plot 

  call pol%init(10)
  call pol%get_Position(Position)
  print *, "Position"
  print *, Position(:,:)

end program
