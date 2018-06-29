pro fitSpinCoeffs, b1a=b1a, b1b=b1b, offset1=offset1, offset2=offset2, spininterval=spininterval

;offset = 8*round(spininterval/8)
;
;b1a = b1a(offset:-1)
;b1b = b1b(0:-offset-1)


XY = [[b1a], [b1b]]
A = [0.000001, 0.000001, -6.6e5, 0.000001, 0.000001, -6.6e5, 0.000001]
Z = dblarr(N_ELEMENTS(b1a))
Z(*) = offset2-offset1

yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='mymulticurvefitthird')
PRINT, 'Function parameters: ', A

b2a = A[0]*b1a^3 + A[1]*b1a^2 + A[2]*b1a
b2b = A[3]*b1b^3 + A[4]*b1b^2 + A[5]*b1b

spinHodogram, b1a=b1a, b1b=b1b, b2a=b2a, b2b=b2b

A = [-6.6e5, -6.6e5 ]

yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='mymulticurvefitfirst')
PRINT, 'Function parameters: ', A

b2a = A[0]*b1a
b2a = A[1]*b2a

spinHodogram, b1a=b1a, b1b=b1b, b2a=b2a, b2b=b2b

end

PRO mymulticurvefitfirst, XY, A, F, PDER
  X = XY[*,0]
  Y = XY[*,1]
  F = A[0]*X - A[1]*Y
  PDER = [ [X], [-Y] ]
END



PRO mymulticurvefitthird, XY, A, F, PDER
  X = XY[*,0]
  Y = XY[*,1]
  F = A[0]*X^3 + A[1]*X^2 + A[2]*X - A[3]*Y^3 - A[4]*Y^2 - A[5]*Y + A[6]
  PDER = [ [X^3], [X^2], [X], [-Y^3], [-Y^2], [-Y], [0] ]
END
