program Poly

  use polymer

  implicit none
  real(8) :: density,temperature
  real(8), allocatable :: Position(:,:)

  type(polymerType) :: pol

  call init_random_seed
  call parameters(density,temperature)
  print *, "Density=",density, "Temperature=",temperature

  call pol%init(40)
  call pol%get_Position(Position)

  call pol%plot
  
end program
