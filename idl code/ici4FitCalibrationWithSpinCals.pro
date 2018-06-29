function ici4FitCalibrationWithSpinCals, bx=bx, by=by, bz=bz, igrfbx=igrfbx, igrfby=igrfby, igrfbz=igrfbz, spinCals=spinCals

  XZCoeffA  = -6.604E-5
  YCoeffA   = -6.604E-5 
  YCoeffB   = 1

  bx = spinCals(0) * (bx-spinCals(3))
  by = spinCals(1) * (by-spinCals(4))
  bz = spinCals(2) * (bz-spinCals(5))


  igrfmag = sqrt(igrfbx^2 + igrfby^2 + igrfbz^2)
  
  ;; Fit first order using spincals.
  print, "Fitting first order function using spin cals: "
  XY = [ [bx], [by], [bz] ]
  Z = (igrfmag)^2
  A = [ XZCoeffA, YCoeffA, YCoeffB ]

  yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='magFitFromSpinFirstOrder', ITER=iterations, ITMAX=maxIterations, STATUS=fitStatus, TOL=fitTolerance, YERROR=yerror, /DOUBLE)

  print, "FitStatus was " + STRING(fitStatus) + " after " + STRING(iterations) + " of " + STRING(maxIterations) + "."
  print, "YError was " + STRING(yerror)
  print, A
  
  XCoeffA = spinCals(0)*A(0)
  XCoeffB = XCoeffA*spinCals(3)
  YCoeffA = spinCals(1)*A(1)
  YCoeffB = YCoeffA*A(2)
  ZCoeffA = spinCals(2)*A(0)
  ZCoeffB = ZCoeffA*spinCals(5)
  
  print, 'Fitted Calibrations are:"
  print, 'X = ' + string(XCoeffA) + ' ' + string(XcoeffB)
  print, 'Y = ' + string(YCoeffA) + ' ' + string(YcoeffB)
  print, 'Z = ' + string(ZCoeffA) + ' ' + string(ZcoeffB)
  
  return, [XCoeffA, YCoeffA, ZCoeffA, XCoeffB, YCoeffB, ZCoeffB]

end

PRO magFitFromSpinFirstOrder, XY, A, F, PDER
  X = XY[*,0]
  Y = XY[*,1]
  Z = XY[*,2]
  F = ( A[0]*X)^2 + (A[1]*(Y - A[2]))^2 + (A[0]*Z)^2
  PDER = [ [2*X*A[0]*X+2*Z*A[0]], [2*Y*(A[1]*Y+A[2])], [2*(A[1]*Y+A[2])]  ]
END
