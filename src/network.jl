function build_lower_network_impedances(id::Int,freq::AbstractVector,cv::CVTree)
  omega = 2π*freq
  nfreq = length(freq)

  Z0data = zeros(ComplexF64,cv.n,nfreq)
  Ztdata = zeros(ComplexF64,cv.n,nfreq)
  ptransdata = zeros(ComplexF64,cv.n,nfreq)

  for (j,wj) in enumerate(omega)
    Z0, Zt, p_trans, check_sum = calculate_impedances_of_segment(wj,id,cv)
    Z0data[:,j] = Z0
    Ztdata[:,j] = Zt
    ptransdata[:,j] = p_trans
  end

  return Z0data, Ztdata, ptransdata

end

build_network_impedances(freq::AbstractVector,cv::CVTree) = build_lower_network_impedances(1,freq,cv)


function calculate_impedances_of_segment(omega::Float64,id::Int,cv::CVTree)
  @unpack n, lengths, radii, Rp, blood_density, blood_viscosity = cv

  Z0 = zeros(ComplexF64,n)
  Zt = zeros(ComplexF64,n)
  p_trans = zeros(ComplexF64,n)

  Zc, kappa = characteristic_impedance_and_wavenumber(omega,id,cv)

  l_k, rad_k = lengths[id], radii[id]

  exp_minus = exp(-im*kappa*l_k)
  exp_plus = 1/exp_minus
  exp_2minus = exp_minus*exp_minus

  check_sum = 1
  if isterminal(id,cv)
    Zt_k = Rp[id]
  else
    Ztinv = complex(0.0)
    for daught_id in daughter_ids(id,cv)
      Z0_d, Zt_d, p_trans_d, check_sum_d = calculate_impedances_of_segment(omega,daught_id,cv)

      Z0 .+= Z0_d
      Zt .+= Zt_d
      p_trans .+= p_trans_d
      Ztinv += 1.0/Z0_d[daught_id]
      check_sum += check_sum_d
    end
    Zt_k = 1.0/Ztinv
  end
  Refl = (Zt_k-Zc)/(Zt_k+Zc)

  if (omega > 0.0)
    Z0_k = Zc*(1.0+Refl*exp_2minus)/(1.0-Refl*exp_2minus)
  else
    Rp = 8.0*blood_viscosity*l_k/(π*rad_k^4)
    Z0_k = Rp + Zc*(1.0+Refl)/(1.0-Refl)
  end

  Z0[id] = Z0_k
  Zt[id] = Zt_k

  # Calculate transmission ratio p(0)/p(-l) for segment k
  p_trans[id] = (1.0+Refl)/(exp_plus + Refl*exp_minus)


  return Z0, Zt, p_trans, check_sum
end


function get_signals_in_segment(id::Int,Q0_mlpersec,Z0data,ptransdata,cv::CVTree)
    @unpack parents = cv

    mmHg_to_Pa = 133.3224
    m3sec_to_mlsec = 100^3


    Ptrans_k = ones(ComplexF64,size(Z0data,2))
    parent_id = parents[id]
    while (parent_id > 0)
      Ptrans_k = Ptrans_k.*ptransdata[parent_id,:]
      parent_id = parents[parent_id]
    end

    Z0_0 = Z0data[1,:]
    Z0_k = Z0data[id,:]

    #=
    # Using pressure input
    p0_Pa = mmHg_to_Pa*p0_mmHg
    p0hat = fft(p0_Pa)
    pkhat = apply_transfer_function(p0hat,Ptrans_k)
    Qkhat = apply_transfer_function(pkhat,1.0./Z0_k)
    =#
    
    # Using flow rate input
    Q0_m3sec = Q0_mlpersec/m3sec_to_mlsec
    Q0hat = fft(Q0_m3sec)
    p0hat = apply_transfer_function(Q0hat,Z0_0)
    pkhat = apply_transfer_function(p0hat,Ptrans_k)
    Qkhat = apply_transfer_function(pkhat,1.0./Z0_k)

    pk_mmHg = real(ifft(pkhat))/mmHg_to_Pa
    Qk_mlpersec = m3sec_to_mlsec*real(ifft(Qkhat))

    return pk_mmHg, Qk_mlpersec

end

function apply_transfer_function(Fhat,That)
    Ghat = zero(Fhat)
    Ghat[1] = Fhat[1]*That[1]
    n = length(Fhat)
    nfreq = n÷2+1
    for k in 2:nfreq-1
        Ghat[k] = Fhat[k]*That[k]
        Ghat[n-k+2] = conj(Ghat[k]) # ensure that the complex conjugate pairs are intact
    end
    return Ghat
end
