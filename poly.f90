program Poly

  use polymer

  implicit none
  real(8) :: temperature, epsilon, sigma
  real(8), allocatable :: Position(:,:)
  integer :: i


  type(polymerType) :: pol

  call init_random_seed
  call parameters(temperature, epsilon, sigma)
  print *, "temperature=",temperature, "epsilon=",epsilon, "sigma", sigma

  do i=1,999
    call pol%init(50,temperature, epsilon, sigma)
    call pol%destroy
    print *, "iteration: ", i
  end do
   
  call pol%init(50,temperature, epsilon, sigma)
  call pol%get_Position(Position)

  call pol%plot
  
end program


