subroutine parameters(temperature, epsilon, sigma)

  real(8), intent(out) :: temperature, epsilon, sigma
  open(10,file="parameters.txt")
  read(10,*) temperature
  read(10,*) epsilon
  read(10,*) sigma
  close(10)

end subroutine
