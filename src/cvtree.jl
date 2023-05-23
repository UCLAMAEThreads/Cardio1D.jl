struct CVTree
  n :: Int
  blood_density :: Float64
  blood_viscosity :: Float64
  poisson_ratio :: Float64
  ve_parameter :: Float64
  names :: Vector{String}
  lengths :: Vector{Float64}
  radii :: Vector{Float64}
  moduli :: Vector{Float64}
  thicknesses :: Vector{Float64}
  wavespeeds :: Vector{Float64}
  Rp :: Vector{Float64}
  parents :: Vector{Int}
  daughters :: Vector{Vector{Int}}
end

function CVTree(names::Vector{String},lengths::Vector{Float64},radii::Vector{Float64},moduli::Vector{Float64},
                  thicknesses::Vector{Float64},parents::Vector{Int},daughters::Vector{Vector{Int}};
                  blood_density = DEFAULT_BLOOD_DENSITY_KGM3,
                  blood_viscosity = DEFAULT_BLOOD_VISCOSITY_KGMS,
                  ve_parameter = DEFAULT_VE_PARAMETER,
                  poisson_ratio = DEFAULT_POISSON_RATIO,
                  terminal_resistance = DEFAULT_TERMINAL_RESISTANCE)
  n = length(names)
  length(lengths) == n || error("Inconsistent lengths vector")
  length(radii) == n || error("Inconsistent radii vector")
  length(moduli) == n || error("Inconsistent moduli vector")
  length(thicknesses) == n || error("Inconsistent thicknesses vector")
  length(parents) == n || error("Inconsistent parents vector")
  length(daughters) == n || error("Inconsistent daughters vector")

  # Transform units to SI
  lengths_m = lengths./100.0
  radii_m = radii./100.0
  thicknesses_m = thicknesses./100.0
  moduli_Pa = moduli*1.0e5

  wavespeeds = _compute_wavespeed.(moduli_Pa,thicknesses_m,radii_m,blood_density)

  Rp = zero(wavespeeds)
  for j in 1:n
    Rp[j] = _set_terminal_resistance(terminal_resistance,Val(isterminal(daughters[j])))
  end

  return CVTree(n,blood_density,blood_viscosity,poisson_ratio,ve_parameter,
                names,lengths_m,radii_m,moduli_Pa,thicknesses_m,wavespeeds,Rp,parents,daughters)


end

function Base.show(io::IO, m::MIME"text/plain", cv::CVTree)
  println(io, "Cardiovascular tree data with $(cv.n) segments")
end


isterminal(segment_daughters::Vector{Int}) = length(segment_daughters) > 0 ? false : true


_set_terminal_resistance(terminal_resistance,::Val{true}) = terminal_resistance
_set_terminal_resistance(terminal_resistance,::Val{false}) = 0.0

### Functions on individual segments

characteristic_impedance_and_wavenumber(omega::Float64,id::Int,cv::CVTree) = _characteristic_impedance_and_wavenumber(omega,cv.blood_density,cv.blood_viscosity,cv.wavespeeds[id],cv.radii[id],
                                            cv.poisson_ratio,cv.ve_parameter)

isterminal(id::Int,cv::CVTree) = isterminal(cv.daughters[id])
ndaughters(id::Int,cv::CVTree) = length(cv.daughters[id])
daughters(id::Int,cv::CVTree) = cv.daughters[id]
