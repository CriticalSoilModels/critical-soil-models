program main
  use stdlib_kinds, only: dp
  use mod_lin_elastic
  use critical_soil_models, only: say_hello
  implicit none

  real(dp) :: test = 0.0_dp
  call say_hello()
end program main
