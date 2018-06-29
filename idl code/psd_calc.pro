function psd_calc, xin=xin, han_win=han_win, res=res, detrend=detrend

  if keyword_set(detrend) then detrend = 1 else detrend = 0
  
  if keyword_set(han_win) then begin
    data_window = hanning(n_elements(xin))
    data_corr=2.67038
  endif else begin
    data_window = dblarr(n_elements(xin))
    data_window[*]=1.0
    data_corr = 1.0
  endelse
  
  if keyword_set(res) then res=res else res=1.0

  
  if detrend then begin
    ; Remove a linear trend from each run of data
    index = dindgen(n_elements(xin))
    coeffs = poly_fit(index, xin, 1)
    xin = xin - poly(index, coeffs)
  endif

  xwin = xin*data_window
  xfft = FFT(xwin)
  xpsd = abs(xfft[0:n_elements(xfft)/2])
  xpsd = xpsd*xpsd*2*data_corr
  xpsd = xpsd/(1./(n_elements(xin)*res))
  
  return, xpsd
end