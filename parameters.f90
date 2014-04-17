subroutine parameters(temperature, epsilon, sigma, numBeads, numAngles, numPopStart, usePerm)

  real(8), intent(out) :: temperature, epsilon, sigma
  integer, intent(out) :: numBeads, numAngles, numPopStart
  logical, intent(out) :: usePerm

  open(10,file="parameters.txt")
  read(10,*) temperature
  read(10,*) epsilon
  read(10,*) sigma
  read(10,*) numBeads
  read(10,*) numAngles
  read(10,*) numPopStart
  read(10,*) usePerm
  close(10)

end subroutine
