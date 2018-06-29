function checkSpinLevels, bx=bx, by=by, bz=bz, samplerate=samplerate, spinFreq=spinFreq, coneFreq=coneFreq


  

  bmag = sqrt( bx^2 + by^2 + bz^2 )
  
  ; FFT_POWERSPECTRUM gets really slow for long datasets
  
  bmag = (N_ELEMENTS(bmag) gt 250000)?bmag(0:250000):bmag
  
  Px = FFT_POWERSPECTRUM(bmag, 1.0/samplerate,FREQ=xfreq, /amp)
  
  spinLevel = Px(value_locate(xfreq,spinFreq))
  coneLevel = Px(value_locate(xfreq,coneFreq))
  
  print, 'Spin level: ' + string(spinLevel) + ' Coning Level: ' + string(coneLevel)

  return, [spinLevel, coneLevel]

end