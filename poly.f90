program Poly
  implicit none
  real(8) :: density,temperature

  call parameters(density,temperature)
  print *, "Density=",density, "Temperature=",temperature

  call plot


end program
