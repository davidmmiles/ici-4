function ici4_fft_bandstop, data=data, lowStopFreq=lowStopFreq, highStopFreq=highStopFreq, cadence=cadence

  
  data_fft = fft( data*hanning(n_elements(data)), -1)
  nfft = n_elements(data_fft)
  
  data_fft[ round((nfft/2)*(lowStopFreq*2.0/cadence)) : round((nfft/2)*(highStopFreq*2.0/cadence)) ] = 1e-4
  data_fft[ round(nfft-(nfft/2)*(highStopFreq*2.0/cadence)-1) : round(nfft-(nfft/2)*(lowStopFreq*2.0/cadence)-1) ] = 1e-4
  
  data_filtered = fft( data_fft, 1 )

	return, data_filtered
	end
