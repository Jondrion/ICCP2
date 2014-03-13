subroutine parameters(density,temperature)

  real(8), intent(out) :: density, temperature
  open(10,file="parameters.txt")
  read(10,*) density
  read(10,*) temperature
  close(10)

end subroutine
