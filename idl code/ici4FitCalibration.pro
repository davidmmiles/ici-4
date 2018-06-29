pro ici4FitCalibration, sampleArray=sampleArray

XCoeffA = 0
XCoeffB = 0
XCoeffC = -6.604E-5
XCoeffD = -69.43
YCoeffA = 0
YCoeffB = 0
YCoeffC = -6.604E-5
YCoeffD = 72.37
ZCoeffA = 0
ZCoeffB = 0
ZCoeffC = -6.604E-5
ZCoeffD = 51.3

xrange = [MIN(sampleArray.time), MAX(sampleArray.time) ]
p_bmagfit = plot(sampleArray.time, sampleArray.B, XRANGE=xrange, COLOR='blue', NAME='Bmag')
p_bmagfit = plot(sampleArray.time, sampleArray.igrfmag, XRANGE=xrange, COLOR='red', NAME='BIGRFmag', /OVERPLOT)


;; Fit first order
print, "Fitting first order function: "
XY = [ [SampleArray.Bxmean], [SampleArray.Bymean], [SampleArray.Bzmean] ]
Z = (sampleArray.igrfmag)^2
A = [ XCoeffC, XCoeffD, YCoeffC, YCoeffD, ZCoeffC, ZCoeffD ]

yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='magFitFirstOrder', ITER=iterations, ITMAX=maxIterations, STATUS=fitStatus, TOL=fitTolerance, YERROR=yerror, /DOUBLE)

print, "FitStatus was " + STRING(fitStatus) + " after " + STRING(iterations) + " of " + STRING(maxIterations) + "."
print, "YError was " + STRING(yerror)

print, "XCoeffA = 0"
print, "XCoeffB = 0"
print, "XCoeffC = " + STRING( A[0] )
print, "XCoeffD = " + STRING( A[1] )

print, "YCoeffA = 0"
print, "YCoeffB = 0"
print, "YCoeffC = " + STRING( A[2] )
print, "YCoeffD = " + STRING( A[3] )

print, "ZCoeffA = 0"
print, "ZCoeffB = 0"
print, "ZCoeffC = " + STRING( A[4] )
print, "ZCoeffD = " + STRING( A[5] )

BxOrder1 = A[0]*(sampleArray.Bxmean) + A[1]
ByOrder1 = A[2]*(sampleArray.Bymean) + A[3]
BzOrder1 = A[4]*(sampleArray.Bzmean) + A[5]

;; Fit second order
print, "Fitting second order function: "
XY = [ [SampleArray.Bxmean], [SampleArray.Bymean], [SampleArray.Bzmean] ]
Z = (sampleArray.igrfmag)^2
A = [ XCoeffB, XCoeffC, XCoeffD, YCoeffB, YCoeffC, YCoeffD, ZCoeffB, ZCoeffC, ZCoeffD ]

yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='magFitSecondOrder', ITER=iterations, ITMAX=maxIterations, STATUS=fitStatus, TOL=fitTolerance, YERROR=yerror, /DOUBLE)

print, "FitStatus was " + STRING(fitStatus)+ " after " + STRING(iterations) + " of " + STRING(maxIterations) + "."
print, "YError was " + STRING(yerror)

print, "XCoeffA = 0"
print, "XCoeffB = " + STRING( A[0] )
print, "XCoeffC = " + STRING( A[1] )
print, "XCoeffD = " + STRING( A[2] )

print, "YCoeffA = 0"
print, "YCoeffB = " + STRING( A[3] )
print, "YCoeffC = " + STRING( A[4] )
print, "YCoeffD = " + STRING( A[5] )

print, "ZCoeffA = 0"
print, "ZCoeffB = " + STRING( A[6] )
print, "ZCoeffC = " + STRING( A[7] )
print, "ZCoeffD = " + STRING( A[8] )

BxOrder2 = A[0]*(sampleArray.Bxmean)^2 + A[1]*(sampleArray.Bxmean) + A[2]
ByOrder2 = A[3]*(sampleArray.Bymean)^2 + A[4]*(sampleArray.Bymean) + A[5]
BzOrder2 = A[6]*(sampleArray.Bzmean)^2 + A[7]*(sampleArray.Bzmean) + A[8]

;; Fit third order
print, "Fitting third order function: "
XY = [ [SampleArray.Bxmean], [SampleArray.Bymean], [SampleArray.Bzmean] ]
Z = (sampleArray.igrfmag)^2
A = [ XCoeffA, XCoeffB, XCoeffC, XCoeffD, YCoeffA, YCoeffB, YCoeffC, YCoeffD, ZCoeffA, ZCoeffB, ZCoeffC, ZCoeffD ]

yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='magFitThirdOrder', ITER=iterations, ITMAX=maxIterations, STATUS=fitStatus, TOL=fitTolerance, YERROR=yerror, /DOUBLE)

print, "FitStatus was " + STRING(fitStatus)+ " after " + STRING(iterations) + " of " + STRING(maxIterations) + "."
print, "YError was " + STRING(yerror)

print, "XCoeffA = " + STRING( A[0] )
print, "XCoeffB = " + STRING( A[1] )
print, "XCoeffC = " + STRING( A[2] )
print, "XCoeffD = " + STRING( A[3] )

print, "YCoeffA = " + STRING( A[4] )
print, "YCoeffB = " + STRING( A[5] )
print, "YCoeffC = " + STRING( A[6] )
print, "YCoeffD = " + STRING( A[7] )

print, "ZCoeffA = " + STRING( A[8] )
print, "ZCoeffB = " + STRING( A[9] )
print, "ZCoeffC = " + STRING( A[10] )
print, "ZCoeffD = " + STRING( A[11] )

BxOrder3 = A[0]*(sampleArray.Bxmean)^3 + A[1]*(sampleArray.Bxmean)^2 + A[2]*(sampleArray.Bxmean) + A[3]
ByOrder3 = A[4]*(sampleArray.Bymean)^3 + A[5]*(sampleArray.Bymean)^2 + A[6]*(sampleArray.Bymean) + A[7]
BzOrder3 = A[8]*(sampleArray.Bzmean)^3 + A[9]*(sampleArray.Bzmean)^2 + A[10]*(sampleArray.Bzmean) + A[11]

BMagOrder1 = sqrt(BxOrder1^2+ByOrder1^2+BzOrder1^2)
BMagOrder2 = sqrt(BxOrder2^2+ByOrder2^2+BzOrder2^2)
BMagOrder3 = sqrt(BxOrder3^2+ByOrder3^2+BzOrder3^2)


p_bmagfit = plot(sampleArray.time, sampleArray.B, XRANGE=xrange, COLOR='black', NAME='Original')
p_bmagfit = plot(sampleArray.time, BmagOrder1, XRANGE=xrange, COLOR='red', NAME='Order1', /OVERPLOT)
p_bmagfit = plot(sampleArray.time, BmagOrder2, XRANGE=xrange, COLOR='blue', NAME='Order2', /OVERPLOT)
p_bmagfit = plot(sampleArray.time, BmagOrder3, XRANGE=xrange, COLOR='green', NAME='Order3', /OVERPLOT)
p_bmagfit = plot(sampleArray.time, sampleArray.igrfmag, XRANGE=xrange, COLOR='pink', NAME='BIGRFmag', /OVERPLOT)


end

PRO magFitFirstOrder, XY, A, F, PDER
  X = XY[*,0]
  Y = XY[*,1]
  Z = XY[*,2]
  F = ( A[0]*X + A[1])^2 + (A[2]*Y + A[3])^2 + (A[4]*Z + A[5] )^2
  PDER = [ [2*X*(A[0]*X+A[1])], [2*(A[0]*X+A[1])], $
    [2*Y*(A[2]*Y+A[3])], [2*(A[2]*Y+A[3])], $
    [2*Z*(A[4]*Z+A[5])], [2*(A[4]*Z+A[5])] ]
END

PRO magFitSecondOrder, XY, A, F, PDER
  X = XY[*,0]
  Y = XY[*,1]
  Z = XY[*,2]
  F = ( A[0]*X^2 + A[1]*X + A[2])^2 + (A[3]*Y^2 + A[4]*Y + A[5])^2 + (A[6]*Z^2 + A[7]*Z + A[8] )^2
  PDER = [ [2*X^2*(A[0]*X^2+A[1]*X+A[2])], [2*X*(A[0]*X^2+A[1]*X+A[2])], [2*(A[0]*X^2+A[1]*X+A[2])], $
    [2*Y^2*(A[3]*Y^2+A[4]*Y+A[5])], [2*Y*(A[3]*Y^2+A[4]*Y+A[5])], [2*(A[3]*Y^2+A[4]*Y+A[5])], $
    [2*Z^2*(A[6]*Z^2+A[7]*Z+A[8])], [2*Z*(A[6]*Z^2+A[7]*Z+A[8])], [2*(A[6]*Z^2+A[7]*Z+A[8])] ]
END

PRO magFitThirdOrder, XY, A, F, PDER
  X = XY[*,0]
  Y = XY[*,1]
  Z = XY[*,2]
  F = ( A[0]*X^3 + A[1]*X^2 + A[2]*X + A[3] )^2 + (A[4]*Y^3 + A[5]*Y^2 + A[6]*Y + A[7] )^2 + (A[8]*Z^3 + A[9]*Z^2 + A[10]*Z + A[11] )^2
  PDER = [ [2*X^3*(A[0]*X^3 + A[1]*X^2 + A[2]*X + A[3])], [2*X^2*(A[0]*X^3 + A[1]*X^2 + A[2]*X + A[3])], [2*X*(A[0]*X^3 + A[1]*X^2 + A[2]*X + A[3])], [2*(A[0]*X^3 + A[1]*X^2 + A[2]*X + A[3])], $
    [2*Y^3*(A[4]*Y^3 + A[5]*Y^2 + A[6]*Y + A[7])], [2*Y^2*(A[4]*Y^3 + A[5]*Y^2 + A[6]*Y + A[7])], [2*Y*(A[4]*Y^3 + A[5]*Y^2 + A[6]*Y + A[7])], [2*(A[4]*Y^3 + A[5]*Y^2 + A[6]*Y + A[7])], $
    [2*Z^3*(A[8]*Z^3 + A[9]*Z^2 + A[10]*Z + A[11])], [2*Z^2*(A[8]*Z^3 + A[9]*Z^2 + A[10]*Z + A[11])], [2*Z*(A[8]*Z^3 + A[9]*Z^2 + A[10]*Z + A[11])], [2*(A[8]*Z^3 + A[9]*Z^2 + A[10]*Z + A[11])] ]



END