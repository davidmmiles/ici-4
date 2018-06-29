function sonogram, xin=xin, tin=tin, nfft=nfft, cadence=cadence, overlap=overlap, detrend=detrend

if  keyword_set(detrend) then detrend = 1

if ~keyword_set(cadence) then cadence = 1

step = round((1-overlap)*nfft)

Coeff = DIGITAL_FILTER(4/nfft, 1, 50., 400, /DOUBLE)
xfiltered = CONVOL(xin, Coeff, /NAN, /NORMALIZE, /EDGE_TRUNCATE)
 
px = MAKE_ARRAY( ceil(double(n_elements(xfiltered)-nfft)/step), nfft/2+1, /double)

for i = 0ul, floor(double(n_elements(xfiltered)-nfft)/step)-1 do begin 
  px[i,*] = psd_calc(xin=xfiltered[i*step:i*step+nfft], /han_win, res=1.0/cadence, detrend=detrend)
endfor

return, px
end