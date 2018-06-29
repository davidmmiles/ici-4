function findPeakFreq, bx=bx, by=by, bz=bz, samplingCadence=samplingCadence, lowFreq=lowFreq, highFreq=highFreq

  maxLen = (N_ELEMENTS(bx)gt 250000)?250000:N_ELEMENTS(bx)-1

; FFT_POWERSPECTRUM gets incredibly slow with large data sets
  bx = bx(0:maxLen)
  by = by(0:maxLen)
  bz = bz(0:maxLen)

  ;; mean/median filter to suppress the spin tone.
  Px = FFT_POWERSPECTRUM(bx, 1.0/samplingCadence,FREQ=xfreq, /amp)
  Py = FFT_POWERSPECTRUM(by, 1.0/samplingCadence,FREQ=yfreq, /amp)
  Pz = FFT_POWERSPECTRUM(bz, 1.0/samplingCadence,FREQ=zfreq, /amp)

  PxMax = max(Px[where(xfreq gt lowFreq and xfreq lt highFreq)],/nan)
  PyMax = max(Py[where(yfreq gt lowFreq and yfreq lt highFreq)],/nan)
  PzMax = max(Pz[where(zfreq gt lowFreq and zfreq lt highFreq)],/nan)

  xFreq = xfreq( where(Px eq PxMax))
  yFreq = yfreq( where(Py eq PyMax))
  zFreq = zfreq( where(Pz eq PzMax))
  avgFreq = (xFreq + yFreq + zFreq) / 3.0

  print, "Frequency Peak at X: " + string(xFreq) + " Y: " + string(yFreq) + " Z: " + string(zFreq) + " Avg: " + string(avgFreq)
  
  return, avgFreq
  
  end