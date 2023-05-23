function get_frequencies(tsec::AbstractVector;heartrate_bpm=Inf)
  Δt = tsec[2] - tsec[1]
  heartperiod_sec = tsec[end]-tsec[1]+Δt
  nt = length(tsec)

  tsec_renorm = copy(tsec)
  freq0_Hz = 1.0/heartperiod_sec
  if !isinf(heartrate_bpm)
    freq0_Hz = heartrate_bpm/60
    heartperiod_sec_renorm = 1.0/freq0_Hz
    tsec_renorm .*= heartperiod_sec_renorm/heartperiod_sec
  end

  nfreq = nt÷2+1
  fmax = (nt/2)*freq0_Hz
  freq = range(0.0,fmax,length=nfreq)
  return freq, tsec_renorm

end
