module Cardio1D

  using FFTW
  using SpecialFunctions
  using UnPack

  const DEFAULT_BLOOD_DENSITY_KGM3 = 1050.0
  const DEFAULT_BLOOD_VISCOSITY_KGMS = 0.004
  const DEFAULT_POISSON_RATIO = 0.5
  const DEFAULT_VE_PARAMETER = 15π/180
  const DEFAULT_TERMINAL_RESISTANCE = 6.0e9

  export CVTree

  include("physical_equations.jl")
  include("cvtree.jl")
  include("network.jl")




end # module Cardio1D
