#=
Computes the Moens-Korteweg wave speed
=#
function _compute_wavespeed(modulus,thickness,radius,blood_density)
    return sqrt(modulus*thickness/(2*blood_density*radius))
end


"""
This routine calculates the characteristic impedance Zc and the complex wavenumber kappa,
given the tree physiological data and the angular frequency `omega` (rad/sec).

The calculation here is based on the approach of Avolio (1980), which itself is adapted from
Westerhof and Noodergraaf (1970).
"""
function _characteristic_impedance_and_wavenumber(omega,blood_density,blood_viscosity,wavespeed,radius,poisson_ratio,ve_parameter)
  alpha = radius*sqrt(omega*blood_density/blood_viscosity)
  area = π*radius^2
  beta = max(alpha,1.0e-3)*im^(3/2)
  F10 = 2.0*besselj(1,beta)/beta/besselj(0,beta)
  F10fact = 1.0/sqrt(1.0-F10)
  phi = ve_parameter*(1.0 - exp(-2.0*omega))
  phi_fact = exp(im*0.5*phi)

  Zc = blood_density*wavespeed/area/sqrt(1.0-poisson_ratio^2)
  Zc *= (omega > 0.0) ? F10fact*phi_fact : real(F10fact)

  kappa = omega/wavespeed*F10fact/phi_fact

  return Zc, kappa
end
