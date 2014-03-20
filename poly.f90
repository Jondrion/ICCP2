program Poly

  use polymer

  implicit none
  real(8) :: density,temperature
  real(8), allocatable :: Position(:,:)

  type(polymerType) :: pol

  call parameters(density,temperature)
  print *, "Density=",density, "Temperature=",temperature

  call pol%init(10)
  call pol%get_Position(Position)
  print *, "Position"
  print *, Position(:,:)
  call pol%plot

end program