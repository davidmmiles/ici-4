function rotateMagDataTimeVarying, magx=magx, magy=magy, magz=magz, theta=theta, roe=roe, phi=phi

  numSamples = N_ELEMENTS(magx)

  magArray = MAKE_ARRAY(3,numSamples, /DOUBLE)

  magArray(0, *) = magx
  magArray(1, *) = magy
  magArray(2, *) = magz

  for j = 0UL, numSamples - 1 do begin
    o = theta(j)*!DtoR
    rotMatrix = [ [ 1, 0, 0], [ 0, COS(o), -SIN(o)], [ 0, SIN(o), cos(o) ] ]
    magArray(*, j) = rotMatrix # magArray(*,j)
  endfor
  
  
;  for j = 0UL, numSamples - 1 do begin
;    o = roe(j)*!DtoR
;    rotMatrix = [ [ COS(o), 0, SIN(o)], [ 0, 1, 0], [ -SIN(o), 0, COS(o)] ]
;    magArray(*, j) = rotMatrix # magArray(*,j)
;  endfor
;
;  for j = 0UL, numSamples - 1 do begin
;    o = phi(j)*!DtoR
;    rotMatrix = [ [ COS(o), -SIN(o), 0], [ SIN(o), COS(o), 0], [ 0, 0, 1] ]
;    magArray(*, j) = rotMatrix # magArray(*,j)
;  endfor
; 



  
  return, magArray
  
end
