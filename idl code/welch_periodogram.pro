function welch_periodogram, xin=xin, tin=tin, nfft=nfft, cadence=cadence, overlap=overlap, detrend=detrend

  if keyword_set(detrend) then detrend = 1 else detrend = 0
  
  if ~keyword_set(cadence) then cadence = 1

  step = round((1-overlap)*nfft)

;  Coeff = DIGITAL_FILTER(4/nfft, 1, 50., 400, /DOUBLE)
;  xfiltered = CONVOL(xin, Coeff, /NAN, /NORMALIZE, /EDGE_TRUNCATE)

  xfiltered = xin

  px = MAKE_ARRAY( nfft/2+1, /double)
  count = 0

  for i = 0ul, floor(double(n_elements(xfiltered)-nfft)/step)-1 do begin
    px += psd_calc(xin=xfiltered[i*step:i*step+nfft], /han_win, res=1.0/cadence, detrend=detrend)
    count++
  endfor

  return, px / count
end