function calSpinPlane, sampleArray=sampleArray, magSamplingCadence=magSamplingCadence

    sampleArray.Bxraw = -sampleArray.Bxraw
    
    print, 'Looking for spin freqency'
    spinFreq = findPeakFreq(bx=sampleArray.bxraw, by=sampleArray.byraw, bz=sampleArray.bzraw, samplingCadence=magsamplingCadence, lowFreq=1, highFreq=5)
    
    print, 'Looking for coning frequency'
    coneFreq = findPeakFreq(bx=sampleArray.bxraw, by=sampleArray.byraw, bz=sampleArray.bzraw, samplingCadence=magsamplingCadence, lowFreq=1.0/20, highFreq=1.0/5)
    
    
    ; B_1 uses the generic calibrations
    
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
    
    Bx_1 = XCoeffA*(sampleArray.Bxraw)^3 + XCoeffB*(sampleArray.Bxraw)^2 + XCoeffC*(sampleArray.Bxraw) + XCoeffD
    By_1 = YCoeffA*(sampleArray.Byraw)^3 + YCoeffB*(sampleArray.Byraw)^2 + YCoeffC*(sampleArray.Byraw) + YCoeffD
    Bz_1 = ZCoeffA*(sampleArray.Bzraw)^3 + ZCoeffB*(sampleArray.Bzraw)^2 + ZCoeffC*(sampleArray.Bzraw) + ZCoeffD
    
    
    levels1 = checkspinLevels(bx=Bx_1, by=By_1, bz=Bz_1, samplerate=magsamplingCadence, spinFreq=spinFreq, coneFreq=coneFreq)
    
    ; Work with spin plane axis to calculate zeros and scalings
    
    print, 'Work with spin plane axis to calculate zeros and scalings'
    

    smoothWidth = round(magSamplingCadence / spinFreq)

    print, "Mean filtering to remove spin with a window of " + string(smoothWidth)

    sampleArray.Bxmean = smooth(sampleArray.Bxraw, smoothwidth, /edge_truncate, /nan)
    sampleArray.Bymean = smooth(sampleArray.Byraw, smoothwidth, /edge_truncate, /nan)
    sampleArray.Bzmean = smooth(sampleArray.Bzraw, smoothwidth, /edge_truncate, /nan)
    
    sampleArray = sampleArray[smoothWidth/2:-smoothWidth/2]
    
    bxOffset = mean(sampleArray.Bxmean)
    bzOffset = mean(sampleArray.Bzmean)
    
    p1 = plot(sampleArray.time, sampleArray.bxmean, 'red')
    p1 = plot(sampleArray.time, sampleArray.bzmean, 'blue', /overplot)
    
    print, "Spin averaging gives Bx offset: " + string(bxOffset) + " and Bz offset: " + string(bzoffset)
    
    p2 = plot(sampleArray.time, sampleArray.bxmean-bxOffset, 'red')
    p2 = plot(sampleArray.time, sampleArray.bzmean-bzOffset, 'blue', /overplot)
    
    b1a = sampleArray.bxraw
    b1b = sampleArray.bzraw
    b2a = sampleArray.bxraw - bxOffset
    b2b = sampleArray.bzraw - bzOffset
    
    print, 'Hodogram with X and Z offsets removed'
    
    hodo1 = spinHodogramCompare(b1a=b1a, b1b=b1b, b2a=b2a, b2b=b2b)
    
    print, 'Comparing X to Z at 90 degree offset to find gain scaling.'
    
    b1a = sampleArray.bxraw - bxOffset
    b1b = sampleArray.bzraw - bzOffset
    b2a = b1a(1:-round(smoothwidth/4))
    b2b = b1b(round(smoothwidth/4):-1)
    
    ratio = b2a/b2b
    index = where(abs(ratio) gt 10)
    
    ratio(index) = !Values.D_NaN
    
    p3 = plot(ratio)
    
    scale = mean(ratio,/nan)
    
    print, 'Scale of X/Z is: ' + string(scale)
    
    b2a = b1a
    b2b = b1b*scale
    
    print, 'Hodogram with X and Z offsets removed, Z scaling adjusted'
    
    hodo2 = spinHodogramCompare(b1a=b1a, b1b=b1b, b2a=b2a, b2b=b2b)
    
    Bx_2 = XCoeffC*(sampleArray.Bxraw-bxOffset)*1
    By_2 = YCoeffC*(sampleArray.Byraw)
    Bz_2 = ZCoeffC*(sampleArray.Bzraw-bzOffset)*scale
    
    levels2 = checkspinLevels(bx=Bx_2, by=By_2, bz=Bz_2, samplerate=magsamplingCadence, spinFreq=spinFreq, coneFreq=coneFreq)
    
    Bx_3 = XCoeffC*(sampleArray.Bxraw-bxOffset-hodo2(4))*1
    By_3 = YCoeffC*(sampleArray.Byraw)
    Bz_3 = ZCoeffC*(sampleArray.Bzraw-bzOffset-hodo2(5))*scale
    
    hodo2 = spinHodogramCompare(b1a=Bx_2, b1b=Bz_2, b2a=Bx_3, b2b=Bz_3)
    
    levels3 = checkspinLevels(bx=Bx_3, by=By_3, bz=Bz_3, samplerate=magsamplingCadence, spinFreq=spinFreq, coneFreq=coneFreq)
    
    
    return, [ 1, 1, scale, -bxOffset-hodo2(4), 0, -bzOffset-hodo2(5)]
    
    stop

;    spinperiod = 4*round(spinperiod/4)
;
;temp = sampleArray[spinperiod/4:-1].bxmean - sampleArray[1:-spinperiod/4].bzmean
;
;p = plot(sampleArray.time, sampleArray.bxmean, 'red')
;p = plot(sampleArray.time, sampleArray.bzmean, 'green', /overplot)
;p = plot(sampleArray.time, temp, 'black', /overplot)




end