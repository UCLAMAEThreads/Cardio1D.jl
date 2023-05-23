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
    for daught_id in daughters(id,cv)
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


function get_pressure_sig(id::Int,freqbpm,P0_mmHg,Z0data,ptransdata,cv::CVTree)
    @unpack parents = cv

    mmHg_to_Pa = 133.3224
    m3sec_to_mlsec = 100^3

    P0_Pa = mmHg_to_Pa*P0_mmHg
    p0hat = fft(P0_Pa)

    Z0_k = Z0data[id,:]
    Ptrans_k = ones(ComplexF64,size(Z0_k))
    parent_id = parents[id]
    while (parent_id > 0)
      Ptrans_k = Ptrans_k.*ptransdata[parent_id,:]
      parent_id = parents[parent_id]
    end
    pkhat = apply_transfer_function(p0hat,Ptrans_k)
    Qkhat = apply_transfer_function(pkhat,1.0./Z0_k)

    pk_mmHg = real(ifft(pkhat))/mmHg_to_Pa
    Qk_mlpersec = m3sec_to_mlsec*real(ifft(Qkhat))

    t = collect(range(0.0,60.0/freqbpm,length(pk_mmHg)))

    return t, pk_mmHg, Qk_mlpersec

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
