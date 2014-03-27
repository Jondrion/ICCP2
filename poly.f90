program Poly

  use polymer

  implicit none
  real(8) :: temperature, epsilon, sigma
  real(8), allocatable :: Position(:,:)


  type(polymerType) :: pol

  call init_random_seed
  call parameters(temperature, epsilon, sigma)
  print *, "temperature=",temperature, "epsilon=",epsilon, "sigma", sigma

 
  call pol%init(30,temperature, epsilon, sigma)
  
 

  
  call pol%get_Position(Position)

  call pol%plot
  
end program


