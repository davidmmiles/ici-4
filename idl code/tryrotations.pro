;; Angles from the boom to the rocket (degrees) FIXME: Unconfirmed and approximate
xAngleBoom  = ( 0 )
yAngleBoom  = ( 0 )
zAngleBoom  = ( 85 )


;; Calculate |B|
sampleArray.B = sqrt( sampleArray.Bx^2 + sampleArray.By^2 + sampleArray.Bz^2 )

;temp = sampleArray.Bx
;sampleArray.Bx = -sampleArray.By
;sampleArray.By = -temp


;; Create data rotated into the rocket (SC) frame
magArray = rotateMagData(magx=sampleArray.By, magy=sampleArray.Bx, magz=sampleArray.Bz, anglex=xAngleBoom, angley=yAngleBoom, anglez=zAngleBoom)
sampleArray.Bscx = reform(magArray(0,*))
sampleArray.Bscy = reform(magArray(1,*))
sampleArray.Bscz = reform(magArray(2,*))

dummyRoeE = MAKE_ARRAY(N_ELEMENTS(sampleArray.thetaE), /DOUBLE, VALUE=0)

rotatedB = rotateMagDataTimeVarying(magx=sampleArray.Bscx, magy=sampleArray.Bscy, magz=sampleArray.Bscz, theta=sampleArray.phiE, roe=dummyRoeE, phi=sampleArray.thetaE)

ici4timeseries, time=sampleArray.time, bx1=rotatedB(0,*), by1=rotatedB(1,*), bz1=rotatedB(2,*), bx2=sampleArray.bscx, by2=sampleArray.bscy, bz2=sampleArray.bscz, spikex=sampleArray.spikex, spikey=sampleArray.spikeY, spikeZ=sampleArray.SpikeZ, coordinates='Rotated', ScriptName=ScriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer

ici4timeseries, time=sampleArray.time, bx1=rotatedB(0,*), by1=rotatedB(1,*), bz1=rotatedB(2,*), bx2=sampleArray.IGRFgeix, by2=sampleArray.IGRFgeiy, bz2=sampleArray.IGRFgeiz, spikex=sampleArray.spikex, spikey=sampleArray.spikeY, spikeZ=sampleArray.SpikeZ, coordinates='Rotated', ScriptName=ScriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer
