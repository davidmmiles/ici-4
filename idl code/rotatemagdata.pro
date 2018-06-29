function rotateMagData, magx=magx, magy=magy, magz=magz, anglex=anglex, angley=angley, anglez=anglez

  magArray = MAKE_ARRAY(3,N_ELEMENTS(magx), /DOUBLE)

  magArray(0, *) = magx
  magArray(1, *) = magy
  magArray(2, *) = magz
  
  o = anglex*!DtoR
  rotMatrix = [ [ 1, 0, 0], [ 0, COS(o), -SIN(o)], [ 0, SIN(o), cos(o) ] ]

  magArray = rotMatrix # magArray

  o = angley*!DtoR
  rotMatrix = [ [ COS(o), 0, SIN(o)], [ 0, 1, 0], [ -SIN(o), 0, COS(o)] ]

  magArray = rotMatrix # magArray

  o = anglez*!DtoR
  rotMatrix = [ [ COS(o), -SIN(o), 0], [ SIN(o), COS(o), 0], [ 0, 0, 1] ]

  magArray = rotMatrix # magArray
  
  return, magArray

end
