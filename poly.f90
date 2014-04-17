program Poly
  use polymer
  implicit none

  real(8) :: temperature, epsilon, sigma, Weight
  integer :: i, numBeads, numAngles, numPopStart
  logical :: usePerm
  type(polymerType) :: pol

  !-- Delete date from previous simulations
  open (5, file="initialpolymer.txt", status="unknown")
  close (5, status="DELETE")
  open (5, file="polymerdata.txt", status="unknown")
  close (5, status="DELETE")

  call init_random_seed
  call parameters(temperature, epsilon, sigma, numBeads, numAngles, numPopStart, usePerm)

  Weight=1._8
  call pol%init(numBeads, temperature, epsilon, sigma, numAngles, numPopStart, usePerm)

  ! -- grow the polymers
  do i=1,numPopStart

    call pol%addBead(Weight,3)    
    call pol%reset
    Weight=1._8

  end do
  
  call pol%destroy
  
end program


