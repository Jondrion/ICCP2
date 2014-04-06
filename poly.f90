program Poly

  use polymer

  implicit none
  real(8) :: temperature, epsilon, sigma, Weight
  integer :: i
  type(polymerType) :: pol

  open (5, file="initialpolymer.txt", status="unknown")
  close (5, status="DELETE")
  open (5, file="polymerdata.txt", status="unknown")
  close (5, status="DELETE")

  call init_random_seed
  call parameters(temperature, epsilon, sigma)
  print *, "temperature=",temperature, "epsilon=",epsilon, "sigma", sigma

  Weight=1._8
  call pol%init(350,temperature, epsilon, sigma)
  do i=1,10000
    
    call pol%create(Weight,3)    
    call pol%reset
    Weight=1._8

    print *, "iteration: ", i
  end do
  
  call pol%destroy
  
end program


