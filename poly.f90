program Poly

  use polymer

  implicit none
  real(8) :: temperature, epsilon, sigma
  real(8), allocatable :: Position(:,:)
  integer :: i
  type(polymerType) :: pol

  open (5, file="initialpolymer.txt", status="unknown")
  close (5, status="DELETE")
  open (5, file="polymerdata.txt", status="unknown")
  close (5, status="DELETE")

  call init_random_seed
  call parameters(temperature, epsilon, sigma)
  print *, "temperature=",temperature, "epsilon=",epsilon, "sigma", sigma

  do i=1,4999
    call pol%init(100,temperature, epsilon, sigma)
    call pol%destroy
    print *, "iteration: ", i
  end do
  
end program


