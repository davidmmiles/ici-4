pro ici4parser, InputMagFileName=InputMagFileName, InputHKFileName=InputHKFileName, startsecond=startsecond, stopsecond=stopsecond, $
  disabledespike=disabledespike, disableinterpol=disableinterpol, disablespinmean=disablespinmean, disableslowplots=disableslowplots, $
  screen=screen, cadence=cadence, fitcalibrations=fitcalibrations

;; User Configurable Variables

; Software version number
VersionMajor  = 0
VersionMinor  = 3
VersionRev    = 0

;; Liftoff time (YYYYY MO DD HH MM SS) - 19.02.2015:22:06:41 UT
loYear    = 2015
loMonth   = 02
loDay     = 19
loHour    = 22
loMinute  = 06
loSecond  = 41.0D

;; Boom deployed time (seconds) - FIXME: Need definitive value
boomDeployedTime = 70

;; Time before liftoff to preserve (seconds)
tMinus = 10

;; Parameters for fitting the magnetic calibration
fitstart = 125
fitstop = 450
maxIterations = 50
fitTolerance = 1e-3

; Earth radius (m)
Re = 6371000

;; Angles from the boom to the rocket (degrees) FIXME: Unconfirmed and approximate
xAngleBoom  = ( 0 ) 
yAngleBoom  = ( 0 ) 
zAngleBoom  = ( 85 )

; Date format string
dateFormat = '(C(CDwA, X, CMoA, X, CDI2.2, X, CHI2.2, ":", CMI2.2, ":", CSF9.6, CYI5))'

; Output file prefix used for all plots
FilePrefix = ""

; The name of this script
ScriptName = "ici4parser"

; AA Filter cutoff (percent of nyquist)
aaCutoff = 0.95

; Resolution to create pngs of plots
PlotResolution = 120

; Spike removal
magDespikeInterval = 100          ; Interval over which to calculate the mean for spike detection for magnetic data
magdespikeThreshold = 10000000      ; Threshold (in counts) above which to declare a spike for magnetic data
hkDespikeInterval = 10           ; Interval over which to calculate the mean for spike detection for housekeeping data
hkdespikeThreshold = 100         ; Threshold (in counts) above which to declare a spike for housekeeping data

; Sampling rates
magSamplingCadence = 2879
hkSamplingCadence = 2879/64

; Sensor Transfer Function (nT) from (counts) = Ax^3 + Bx^2 + Cx + d

;XCoeffA = 0
;XCoeffB = 0
;XCoeffC = -6.604E-5
;XCoeffD = -69.43
;YCoeffA = 0
;YCoeffB = 0
;YCoeffC = -6.604E-5
;YCoeffD = 72.37
;ZCoeffA = 0
;ZCoeffB = 0
;ZCoeffC = -6.604E-5
;ZCoeffD = 51.3

XCoeffA = 0
XCoeffB = 0
XCoeffC =  -6.1079159e-005
XCoeffD =       4094.7108
YCoeffA = 0
YCoeffB = 0
YCoeffC =  -5.8944181e-005
YCoeffD =        0.093014254
ZCoeffA = 0
ZCoeffB = 0
ZCoeffC =  -8.2365302e-005
ZCoeffD =       774.37546


;XCoeffA = 0
;XCoeffB = 0
;XCoeffC =  -6.1714967e-005
;XCoeffD =       -1828.1928
;YCoeffA = 0
;YCoeffB = 0
;YCoeffC =  -6.0746401e-005
;YCoeffD =        1287.5491
;ZCoeffA = 0
;ZCoeffB = 0
;ZCoeffC =  -8.3865916e-005
;ZCoeffD =       -409.31033

;XCoeffA = 0
;XCoeffB = 0
;XCoeffC =  -6.1543522e-005
;XCoeffD =       -1822.1193
;YCoeffA = 0
;YCoeffB = 0
;YCoeffC =  -6.0385900e-005
;YCoeffD =        1005.5572
;ZCoeffA = 0
;ZCoeffB = 0
;ZCoeffC =  -8.3633321e-005
;ZCoeffD =       -408.55314

;XCoeffA = 0
;XCoeffB =   5.1781678e-016
;XCoeffC =  -6.1490307e-005
;XCoeffD =       -1872.6825
;YCoeffA = 0
;YCoeffB =   3.0574379e-015
;YCoeffC =  -6.4699347e-005
;YCoeffD =        2525.1808
;ZCoeffA = 0
;ZCoeffB =   2.2917662e-016
;ZCoeffC =  -8.3679586e-005
;ZCoeffD =       -423.49710


; Dynamic Spectra parameters
dynSpectraNFFT    = 512
dynSpectraOverlap = 0.50
welchNFFT         = 4096
welchOverlap      = 0.50 

; Housekeeping Transfer Functions (Various Units) from (counts) = Cx + d

HK0CoeffC = 610.501E-6
Hk0CoeffD = 0

HK1CoeffC = 122.1E-3
Hk1CoeffD = -273.15

HK2CoeffC = 31.3077E-3
Hk2CoeffD = -20.51

HK3CoeffC = 610.501E-6
Hk3CoeffD = 0

HK4CoeffC = 1.8315E-3
Hk4CoeffD = 0

HK5CoeffC = 9.6432E-3
Hk5CoeffD = 1.0136

HK6CoeffC = 31.3077E-3
Hk6CoeffD = -20.51

HK7CoeffC = 64.8466E-3
Hk7CoeffD = 2.02814

PlotFontSize = 12

VersionString = "V" + STRING(VersionMajor, FORMAT='(I0)') + "." + STRING(VersionMinor, FORMAT='(I0)') + "." + STRING(VersionRev, FORMAT='(I0)')


if keyword_set(startsecond) then begin
  print, "Keyword startsecond set to " + STRING(startSecond, FORMAT='(I0)') + "."
endif else begin
  print, "Keyword startsecond not set. Will default to start of data"
endelse

if keyword_set(stopsecond) then begin
  print, "Keyword stopsecond set to " + STRING(stopSecond, FORMAT='(I0)') + "."
endif else begin
  print, "Keyword stopsecond not set. Will default to end of data"
endelse

if(not keyword_set(screen)) then buffer = 1 else buffer = 0
  
if not keyword_set(cadence) then cadence = magSamplingCadence

if cadence gt magSamplingCadence then begin
  print, "ERROR: Won't upsample magnetic data from " + STRING(magSamplingCadence) + " to " + STRING(cadence) + "."
endif else if cadence lt magSamplingCadence then begin
  if cadence ne round(cadence) then begin
    print, "Rounding cadence to " + STRING( round(cadence),FORMAT='(I0)') + " sps."
    cadence = round(cadence)
  endif
  
  print, "Magnetic data will be downsampled from " + STRING(magSamplingCadence,FORMAT='(I0)') + " to " + STRING(cadence,FORMAT='(I0)') + " sps."
endif else begin
  print, "Magdnetic data left at native cadence of " + STRING(cadence,FORMAT='(I0)') + " sps."
endelse

print, "Running " + ScriptName + " " + VersionString
tic

if keyword_set(InputMagFileName) then begin
  print, "Parsing file: " + InputMagFileName + " for magnetic data." 
  
  ;; Parse the raw magnetic telemtry file
  worstCaseMagLength = file_lines(InputMagFileName)
  
  print, InputMagFileName + " contains approximately " + STRING(worstCaseMagLength) + " records or " + STRING(worstCaseMagLength/magSamplingCadence) + " seconds of data."
  
  a = {magSample, index:double(0), time:double(0), date:double(0), $
    Bxraw:double(0), Byraw:double(0), Bzraw:double(0), Bxmean:double(0), Bymean:double(0), Bzmean:double(0), $
    Bx:double(0), By:double(0), Bz:double(0), B:double(0), $
    Bscx:double(0), Bscy:double(0), Bscz:double(0), Bgeix:double(0), Bgeiy:double(0), Bgeiz:double(0), $
    alt:double(0), lat:double(0), lon:double(0), thetaE:double(0), phiE:double(0), $
    IGRFgeix:double(0), IGRFgeiy:double(0), igrfgeiz:double(0), igrfmag:double(0), $
    spikeX:float(0), spikeY:float(0), spikeZ:float(0)}
    
  rocketTime = double(0)
  time = double(0)
  Bxraw = double(0)
  Byraw = double(0)
  Bzraw = double(0)
  
  sampleArray = replicate({magSample}, worstCaseMagLength)
  
  sampleCounter = ulong(0)
  
  line = ''
  
  openr, inputFileUnit, InputMagFileName, /GET_LUN
  readf, inputFileUnit, line ; Skip the header line
  
  if keyword_set(notimestamp) then begin
  
    while not eof(inputFileUnit) do begin
      readf, inputFileUnit, line
      
      reads, strmid(line,0, 400), rocketTime, Bxraw, Byraw, Bzraw
    
      sampleArray[sampleCounter].index = sampleCounter
      sampleArray[sampleCounter].Bxraw = Bxraw
      sampleArray[sampleCounter].Byraw = Byraw
      sampleArray[sampleCounter].Bzraw = Bzraw
      
      sampleCounter++
      
    endwhile
  endif else begin
    while not eof(inputFileUnit) do begin
      readf, inputFileUnit, line
      
      
      magString = STRMID(line,STREGEX(line, STRING(9b)))
      timeString = STRMID(line, 0, STREGEX(line, STRING(9b)))    

      reads, strmid(magString,0, 400), Bxraw, Byraw, Bzraw

      sampleArray[sampleCounter].index = sampleCounter
      sampleArray[sampleCounter].Bxraw = Bxraw
      sampleArray[sampleCounter].Byraw = Byraw
      sampleArray[sampleCounter].Bzraw = Bzraw
      
      if STRLEN(timeString) lt 10 then begin
        sampleArray[sampleCounter].date = !VALUES.F_NAN
      endif else begin
        sampleArray[sampleCounter].date = JULDAY(loMonth,loDay,loYear,STRMID(timeString, 0, 2),STRMID(timeString, 3, 2),STRMID(timeString, 6, 9))
      endelse
      
      sampleCounter++

    endwhile
  endelse    
        
  close, inputFileUnit
  free_lun, inputFileUnit

  ;; Cull to only the read samples
  sampleArray=sampleArray[0:sampleCounter-1]
  
  ;; Set transformed data to NaN for now
  sampleArray.Bscx = !VALUES.F_NAN
  sampleArray.Bscy = !VALUES.F_NAN
  sampleArray.Bscz = !VALUES.F_NAN
  sampleArray.Bgeix = !VALUES.F_NAN
  sampleArray.Bgeiy = !VALUES.F_NAN
  sampleArray.Bgeiz = !VALUES.F_NAN
  sampleArray.IGRFgeix = !VALUES.F_NAN
  sampleArray.IGRFgeiz = !VALUES.F_NAN
  sampleArray.IGRFgeiz = !VALUES.F_NAN
    
  ;; Deal with timestamp creation
  print, "Parsed " + string(sampleCounter) + " magnetic samples or " + string(sampleCounter/magSamplingCadence) + " seconds of data"
  
  magStartDate = MIN(sampleArray.date, /NAN)
  magStopDate = MAX(sampleArray.date, /NAN)
  print, "Parsed magnetic timestamps from " + STRING(magStartDate, FORMAT=dateFormat  ) + " to " + STRING( magStopDate, FORMAT=dateFormat ) + "."

  ;; Linear interpolation of timestamps to all samples using julian date. Note: High sample rate means this has slight numerical precision issues.
  dateIndex = where(finite(sampleArray.date))
  dateFit = ladfit(sampleArray[dateIndex].index, sampleArray[dateIndex].date, /DOUBLE)
  sampleArray.date = dateFit[0] + dateFit[1]*sampleArray.index
  
  ;; Date is often referenced to Seconds After Liftoff. Create this time index as well.
  
  loDate = JULDAY(loMonth,loDay,loYear,loHour,loMinute,loSecond)
  print, "Liftoff declared to be " + STRING(loDate, FORMAT=dateFormat) + "."
  
  loIndex = Value_Locate(sampleArray.date, loDate)
  sampleArray.time = 24D*60D*60D*( (dateFit(0) - sampleArray(loIndex).date) + dateFit(1)*sampleArray.index )
  
  magStartTime = MIN(sampleArray.time)
  magStopTime = MAX(sampleArray.time)
  print, "Data spans range of " + STRING(magStartTime) + " to " + string(magStopTime) + " seconds from liftoff."

  print, "Culling data to liftoff minus " + STRING(tMinus, FORMAT='(I0)') + "s and later."
  sampleArray = sampleArray[ where(sampleArray.time ge -1*tMinus) ]
  
  ;; Data was logged as 32 bit binary but is actualy 32 bit two's complement
  SampleArray.Bxraw = long(ulong(sampleArray.Bxraw))
  SampleArray.Byraw = long(ulong(sampleArray.Byraw))
  SampleArray.Bzraw = long(ulong(sampleArray.Bzraw))
  
  
  FilePrefix = "_" + STRING(magStartTime,FORMAT='(I0)') + "s_to_" + STRING(magStopTime,FORMAT='(I0)') + "s_" + STRING(cadence,FORMAT='(I0)') + "sps_" + VersionString + "_"

  ;; Look for and remove spikes (probably telemetery errors?)
  if ~keyword_set(disabledespike) then begin

    print, "Despiking magnetic data with a window of " + string(magDespikeInterval) + " and a threshold of " + string(magdespikeThreshold) + " counts"

    temp = sampleArray.Bxraw
    spikeIndexX = where(abs(temp-median(temp,magDespikeInterval)) gt magdespikeThreshold)
    if spikeIndexX eq [-1] then spikeCountX = 0 else spikeCountX = N_ELEMENTS(spikeIndexX)
    SampleArray[spikeIndexX].Bxraw = !Values.D_NaN
    sampleArray.spikeX = 5
    sampleArray[spikeIndexX].spikeX = 6

    temp = sampleArray.Byraw
    spikeIndexY = where(abs(temp-median(temp,magDespikeInterval)) gt magdespikeThreshold)
    if spikeIndexY eq [-1] then spikeCountY = 0 else spikeCountY= N_ELEMENTS(spikeIndexY)
    SampleArray[spikeIndexY].Byraw = !Values.D_NaN
    sampleArray.spikeY = 3
    sampleArray[spikeIndexy].spikeY = 4

    temp = sampleArray.Bzraw
    spikeIndexZ = where(abs(temp-median(temp,magDespikeInterval)) gt magdespikeThreshold)
    if spikeIndexZ eq [-1] then spikeCountZ = 0 else spikeCountZ = N_ELEMENTS(spikeIndexZ)
    SampleArray[spikeIndexZ].Bzraw = !Values.D_NaN
    sampleArray.spikeZ = 1
    sampleArray[spikeIndexZ].spikeZ = 2

    print, "Despiking found " + string(spikeCountX) + " in Channel X. Replaced with NAN for future interpolation."
    print, "Despiking found " + string(spikeCountY) + " in Channel Y. Replaced with NAN for future interpolation."
    print, "Despiking found " + string(spikeCountZ) + " in Channel Z. Replaced with NAN for future interpolation."
  endif else begin
    print, "Skipping despike of magnetic data."
  endelse

  if ~keyword_set(disableinterpol) then begin
    sampleArray[where( finite(sampleArray.Bxraw, /NAN) )].Bxraw = interpol( sampleArray[where( finite(sampleArray.Bxraw) )].Bxraw, sampleArray[where( finite(sampleArray.Bxraw) )].index, sampleArray[where( finite(sampleArray.Bxraw, /NAN) )].index)
    sampleArray[where( finite(sampleArray.Byraw, /NAN) )].Byraw = interpol( sampleArray[where( finite(sampleArray.Byraw) )].Byraw, sampleArray[where( finite(sampleArray.Byraw) )].index, sampleArray[where( finite(sampleArray.Byraw, /NAN) )].index)
    sampleArray[where( finite(sampleArray.Bzraw, /NAN) )].Bzraw = interpol( sampleArray[where( finite(sampleArray.Bzraw) )].Bzraw, sampleArray[where( finite(sampleArray.Bzraw) )].index, sampleArray[where( finite(sampleArray.Bzraw, /NAN) )].index)

  endif else begin
    print, "Skipping interpolation of magnetic data"
  endelse
  
;  spinFreq = findPeakFreq(bx=sampleArray.bxraw, by=sampleArray.byraw, bz=sampleArray.bzraw, samplingCadence=MagSamplingCadence, lowFreq=1, highFreq=5)
;  smoothWidth = round(magSamplingCadence / spinFreq)
;
;  sampleArray.Bxmean = smooth(sampleArray.Bxraw, smoothwidth, /edge_truncate, /nan)
;  sampleArray.Bymean = smooth(sampleArray.Byraw, smoothwidth, /edge_truncate, /nan)
;  sampleArray.Bzmean = smooth(sampleArray.Bzraw, smoothwidth, /edge_truncate, /nan)
;
;  sampleArray = sampleArray[smoothWidth/2:-smoothWidth/2]
    
        
;  spinCals = calSpinPlane(sampleArray = sampleArray[where(sampleArray.time ge fitstart and sampleArray.time lt fitstop)], magSamplingCadence=magSamplingCadence)

  
  ;; Transform to Physical Units
  ; (nT) from (counts) = Ax^3 + Bx^2 + Cx + d

  sampleArray.Bx = XCoeffA*(sampleArray.Bxmean)^3 + XCoeffB*(sampleArray.Bxmean)^2 + XCoeffC*(sampleArray.Bxmean) + XCoeffD
  sampleArray.By = YCoeffA*(sampleArray.Bymean)^3 + YCoeffB*(sampleArray.Bymean)^2 + YCoeffC*(sampleArray.Bymean) + YCoeffD
  sampleArray.Bz = ZCoeffA*(sampleArray.Bzmean)^3 + ZCoeffB*(sampleArray.Bzmean)^2 + ZCoeffC*(sampleArray.Bzmean) + ZCoeffD

  ;; Downsampling

  if cadence ne magSamplingCadence then begin

    print, "Downsampling from " + STRING(magSamplingCadence) + " to " + STRING(cadence) + " sps."

    ;; Align downsampled data
    outputStartTime = round(min(sampleArray.time))
    outputStopTime = round(max(sampleArray.time))

    ;  Low pass antialiasing filter
    aaFilterFreq = aaCutoff * 0.5 * cadence / magSamplingCadence

    print, "Antialias filter running at " + STRING(aaFilterFreq*magSamplingCadence) + " Hz."
    Coeff = DIGITAL_FILTER(0., (1.0*aaCutoff*cadence)/magSamplingCadence, 100., 1000, /DOUBLE)

    if n_elements(Coeff) ge (sampleCounter) then begin
      print, "ERROR: Number of sample (" + STRING(sampleCounter) + ")  is less than the length of the filter kernel (" +  STRING(N_ELEMENTS(Coeff), FORMAT='(I0)') + ") aborting."
      stop
    endif

    filtered_Bx = CONVOL(sampleArray.Bx, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_By = CONVOL(sampleArray.By, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_Bz = CONVOL(sampleArray.Bz, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_BxRaw = CONVOL(sampleArray.BxRaw, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_ByRaw = CONVOL(sampleArray.ByRaw, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_BzRaw = CONVOL(sampleArray.BzRaw, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_Bxmean = CONVOL(sampleArray.Bxmean, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_Bymean = CONVOL(sampleArray.Bymean, Coeff, /NAN, /EDGE_TRUNCATE)
    filtered_Bzmean = CONVOL(sampleArray.Bzmean, Coeff, /NAN, /EDGE_TRUNCATE)

    lengthDownsampled = (outputStopTime - outputStartTime) * cadence

    if lengthDownsampled le 0 then begin
      print, "Length of downsampled data is less than zero. Aborting."
      stop
    endif

    resampledArray = replicate({magSample},lengthDownsampled)
    resampledArray.time = dindgen(lengthDownsampled) / (lengthDownsampled) * (outputStopTime - outputStartTime) + outputStartTime
    resampledArray.date = INTERPOL( sampleArray.date, sampleArray.time, resampledArray.time)
    resampledArray.Bx = INTERPOL( filtered_Bx, sampleArray.time, resampledArray.time)
    resampledArray.By = INTERPOL( filtered_By, sampleArray.time, resampledArray.time)
    resampledArray.Bz = INTERPOL( filtered_Bz, sampleArray.time, resampledArray.time)
    resampledArray.BxRaw = INTERPOL( filtered_BxRaw, sampleArray.time, resampledArray.time)
    resampledArray.ByRaw = INTERPOL( filtered_ByRaw, sampleArray.time, resampledArray.time)
    resampledArray.BzRaw = INTERPOL( filtered_BzRaw, sampleArray.time, resampledArray.time)
    resampledArray.Bxmean = INTERPOL( filtered_Bxmean, sampleArray.time, resampledArray.time)
    resampledArray.Bymean = INTERPOL( filtered_Bymean, sampleArray.time, resampledArray.time)
    resampledArray.Bzmean = INTERPOL( filtered_Bzmean, sampleArray.time, resampledArray.time)

    print, "WARNING: Despike data invalidated by downsampling!"
    sampleArray = resampledArray
    sampleArray.Bscx = !VALUES.F_NAN
    sampleArray.Bscy = !VALUES.F_NAN
    sampleArray.Bscz = !VALUES.F_NAN
    sampleArray.Bgeix = !VALUES.F_NAN
    sampleArray.Bgeiy = !VALUES.F_NAN
    sampleArray.Bgeiz = !VALUES.F_NAN
    sampleArray.IGRFgeix = !VALUES.F_NAN
    sampleArray.IGRFgeiz = !VALUES.F_NAN
    sampleArray.IGRFgeiz = !VALUES.F_NAN
  endif
  
  ;; Calculate |B|
  sampleArray.B = sqrt( sampleArray.Bx^2 + sampleArray.By^2 + sampleArray.Bz^2 )
  
  ;; Create data rotated into the rocket (SC) frame
  magArray = rotateMagData(magx=sampleArray.Bx, magy=sampleArray.By, magz=sampleArray.Bz, anglex=xAngleBoom, angley=yAngleBoom, anglez=zAngleBoom)
  sampleArray.Bscx = reform(magArray(0,*))
  sampleArray.Bscy = reform(magArray(1,*))
  sampleArray.Bscz = reform(magArray(2,*))
  
  
  ;; Read in housekeeping magnetometers
  
  hkb = {hkmagSample, index:double(0), time:double(0), bx:double(0), by:double(0), bz:double(0) }
  
  time = double(0)
  bx = double(0)
  by = double(0)
  bz = double(0)
  
  hkMagFilename = "mag_data.txt"
  worstCaseHKMagLength = file_lines(hkMagFilename)
  print, hkMagFilename + " contains approximately " + STRING(worstCaseHKMagLength) + " records."
  hkMagArray = replicate({hkmagSample}, worstCaseHKMagLength)

  hkCounter = ulong(0)
  
  
  ;; Read in flight solution
  
  b = {trajectorySample, index:double(0), time:double(0), alt:double(0), lat:double(0), lon:double(0), thetaE:double(0), phiE:double(0) }
  time= double(0)
  thetaEbody = double(0)
  phiEbody = double(0)
  thetaEmag = double(0)
  phiEmag = double(0)
  lat = double(0)
  lon = double(0)
  alt = double(0)
  magx = double(0)
  magy = double(0)
  magz = double(0)
  eddy = double(0)
  
  
  
  if (cadence eq 100) then begin
    ;; Parse the raw attitude telemetry file
    trajectoryFileName = "ICI420151022.dat"
    worstCaseTrajectoryLength = file_lines(trajectoryFileName)
    print, trajectoryFileName + " contains approximately " + STRING(worstCaseTrajectoryLength) + " records or " + STRING(worstCaseTrajectoryLength/1000) + " seconds of data."
      
    trajectoryArray = replicate({trajectorySample}, worstCaseTrajectoryLength)
  
    trajectoryCounter = ulong(0)
  
    line = ''
  
    openr, inputFileUnit2, trajectoryFileName, /GET_LUN
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
    readf, inputFileUnit2, line ; Skip the header line
  
    while not eof(inputFileUnit2) do begin
      readf, inputFileUnit, line
  
      reads, strmid(line,0, 400), time, thetaEbody, phiEbody, thetaEmag, pieEmag, lat, lon, alt, magx, magy, magz, eddy 
  
      trajectoryArray[trajectoryCounter].index = trajectoryCounter
      trajectoryArray[trajectoryCounter].time = time
      trajectoryArray[trajectoryCounter].alt = alt
      trajectoryArray[trajectoryCounter].lat = lat
      trajectoryArray[trajectoryCounter].lon = lon
      trajectoryArray[trajectoryCounter].thetaE = thetaEbody
      trajectoryArray[trajectoryCounter].phiE = phiEbody  
  
      trajectoryCounter++
  
    endwhile
    
    ;; Cull to only the read samples
    trajectoryArray = trajectoryArray[0:trajectoryCounter-1]
    
    ; Align sampleArray to start where we have trajectory data
    sampleArray = sampleArray(where(sampleArray.time eq trajectoryArray(0).time):*)
    
    tempEnd = where(trajectoryArray.time eq sampleArray(-1).time)
    
    sampleArray.alt = trajectoryArray(0:tempEnd).alt
    sampleArray.lat = trajectoryArray(0:tempEnd).lat
    sampleArray.lon = trajectoryArray(0:tempEnd).lon
    sampleArray.thetaE = trajectoryArray(0:tempEnd).thetaE
    sampleArray.phiE = trajectoryArray(0:tempEnd).phiE
    
  ;  dummyRoeE = MAKE_ARRAY(N_ELEMENTS(sampleArray.thetaE), /DOUBLE, VALUE=0)
  ;   
  ;  rotatedB = rotateMagDataTimeVarying(magx=sampleArray.Bscx, magy=sampleArray.Bscy, magz=sampleArray.Bscz, theta=sampleArray.thetaE, roe=dummyRoeE, phi=sampleArray.phiE)
endif else begin
  print, "ERROR:  needs to be at same cadence as trajectory data (currently 100 Hz)
  sampleArray.alt = !VALUES.F_NAN
  sampleArray.lat = !VALUES.F_NAN
  sampleArray.lon = !VALUES.F_NAN
  sampleArray.thetaE = !VALUES.F_NAN
  sampleArray.phiE = !VALUES.F_NAN
endelse
  

  
  ;; Use geopack to calculate IGRF at each sample
  print, "Calculating IGRF at each sample"
  for j = 0UL, N_ELEMENTS(sampleArray) - 1 do begin 
    caldat, sampleArray[j].date, Month, Day, Year, Hour, Minute, Second
    geopack_recalc, Year, Month, Day, Hour, Minute, Second, /date
    geopack_igrf_geo, (sampleArray[j].alt+Re)/Re, 90-sampleArray[j].lat, sampleArray[j].lon, br, btheta, bphi, /DEGREE
    sampleArray[j].IGRFgeix = br
    sampleArray[j].IGRFgeiy = btheta
    sampleArray[j].IGRFgeiz = bphi
    sampleArray[j].igrfmag = (br^2 + btheta^2 + bphi^2)^0.5
;    geopack_conv_coord, br, btheta, bphi, bgeix, bgeiy, bgeiz, /FROM_GEO, /TO_GEI
;    sampleArray[j].IGRFgeix = bgeix
;    sampleArray[j].IGRFgeiy = bgeiy
;    sampleArray[j].igrfgeiz = bgeiz
;    sampleArray[j].igrfmag = (bgeix^2 + bgeiy^2 + bgeiz^2)^0.5
  endfor

  ;; Fit magnetometer calibrations to flight data
 
  if keyword_set(fitcalibrations) then begin
    print, 'Fitting calibration coefficients.'
    sampleArray2 = sampleArray[where(sampleArray.time ge fitstart and sampleArray.time le fitstop)]
    spinCals2 = ici4FitCalibrationWithSpinCals(bx=sampleArray2.bxraw, by=sampleArray2.byraw, bz=sampleArray2.bzraw, igrfbx=sampleArray2.IGRFgeix, igrfby=sampleArray2.IGRFgeiy, igrfbz=sampleArray2.IGRFgeiz, spinCals=spinCals)
    
    print, 'Looking for spin freqency'
    spinFreq = findPeakFreq(bx=sampleArray.bxraw, by=sampleArray.byraw, bz=sampleArray.bzraw, samplingCadence=cadence, lowFreq=1, highFreq=5)

    print, 'Looking for coning frequency'
    coneFreq = findPeakFreq(bx=sampleArray.bxraw, by=sampleArray.byraw, bz=sampleArray.bzraw, samplingCadence=cadence, lowFreq=1.0/20, highFreq=1.0/5)
    
    BxCal = spinCals2(0)*sampleArray.bxraw+spinCals2(3)
    ByCal = spinCals2(1)*sampleArray.byraw+spinCals2(4)
    BzCal = spinCals2(2)*sampleArray.bzraw+spinCals2(5)
    
    Levels1 = checkSpinLevels(bx=bxcal,by=bycal,bz=bzcal,samplerate=cadence,spinfreq=spinfreq,conefreq=conefreq)
    
    BxCal = spinCals2(0)*sampleArray.bxraw-spinCals2(3)
    ByCal = spinCals2(1)*sampleArray.byraw-spinCals2(4)
    BzCal = spinCals2(2)*sampleArray.bzraw-spinCals2(5)

    Levels2 = checkSpinLevels(bx=bxcal,by=bycal,bz=bzcal,samplerate=cadence,spinfreq=spinfreq,conefreq=conefreq)
    
    
    
    ici4FitCalibration, sampleArray=sampleArray2
  endif
  
  ;; Write out processed data
  outputFileName = "ICI4 FGM Data " + STRING(cadence, FORMAT='(I0)') + "sps " + VersionString + ".asc"
  print, "Writing processed data to: " + outputFileName
  openw, outputFile, outputFileName, /get_lun
  printf, outputFile, "ICI-4 FGM Data at "  + STRING(cadence, FORMAT='(I0)') + "sps. Processed by " + ScriptName + " " + VersionString + " at " + systime(/utc) + " UTC"
  printf, outputFile, "Liftoff (0s) set to " + STRING(loDate, FORMAT=dateFormat) + ". Boom considered fully deployed at " + STRING(boomDeployedTime, FORMAT='(F14.6)') + " s."
  printf, outputFile, "Data provided in the frames of the sensor (SNSR), rocket (SC), and Geocentric Equatorial Inertial (GEI)
  printf, outputFile, "Data before liftoff minus " + STRING(tMinus, FORMAT='(I0)') + "s discarded"
  printf, outputFile, "Data follows this line. "
  printf, outputFile, "     Date (Julian)      Time (s)   B_XSNSR (nT)   B_YSNSR (nT)   B_ZSNSR (nT)     B_XSC (nT)     B_YSC (nT)     B_ZSC (nT)    B_XGEI (nT)    B_YGEI (nT)    B_ZGEI (nT) IGRF_XGEI (nT) IGRF_YGEI (nT) IGRF_ZGEI (nT)"
  for j = 0UL, size(sampleArray, /N_ELEMENTS) - 1 do begin
    printf, outputFile, sampleArray[j].date, sampleArray[j].time, sampleArray[j].Bx, sampleArray[j].By, sampleArray[j].Bz, $
      sampleArray[j].Bscx, sampleArray[j].Bscy, sampleArray[j].Bscz, sampleArray[j].Bgeix, sampleArray[j].Bgeiy, sampleArray[j].Bgeiz, $
      sampleArray[j].IGRFgeix, sampleArray[j].IGRFgeiy, sampleArray[j].IGRFgeiz, $
      FORMAT='(F18.10, F14.6, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3, F15.3)'
  endfor
  close, outputFile
  free_lun, outputFile
    
  ;; Truncate data to specified start and stop times

  if keyword_set( startSecond ) then begin
    if (startSecond gt magStartTime) and (startSecond lt magStopTime) then begin
      magStartTime = startSecond
    endif else begin
      print, "StartSecond not possible in magnetic times." + STRING(startSecond)
    endelse
  endif

  if keyword_set( stopSecond ) then begin
    if (stopSecond gt magStartTime) and (stopSecond lt magStopTime) then begin
      magStopTime = stopSecond
    endif else begin
      print, "StopSecond not possible in magnetic times. " + STRING(stopSecond)
    endelse
  endif

  print, "Setting magnetic timespan set to ", + STRING(magStartTime) + " to " + STRING(magStopTime)

  sampleArray = sampleArray[where(sampleArray.time ge magStartTime and sampleArray.time le magStopTime)]

  ;; Timeseries plotting
  print, "Creating timeseries summary plot"

  ici4timeseries, time=sampleArray.time, bx1=sampleArray.bx, by1=sampleArray.by, bz1=sampleArray.bz, bx2=sampleArray.IGRFgeix, by2=sampleArray.IGRFgeiy, bz2=sampleArray.IGRFgeiz, spikex=sampleArray.spikex, spikey=sampleArray.spikeY, spikeZ=sampleArray.SpikeZ, coordinates='Sensor', ScriptName=ScriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer
 
 
  ;; Extract deviation of |B| from |B| IGRF
  
  
  Bresidual = sampleArray.B - sampleArray.igrfmag
  
  Bresidual2 = ici4_fft_bandstop(data=Bresidual, lowStopFreq=3, highStopFreq=4, cadence=100)
  
 
  
;  ici4timeseries, time=sampleArray.time, bx1=sampleArray.bscx, by1=sampleArray.bscy, bz1=sampleArray.bscz, bx2=sampleArray.IGRFgeix, by2=sampleArray.IGRFgeiy, bz2=sampleArray.IGRFgeiz, spikex=sampleArray.spikex, spikey=sampleArray.spikeY, spikeZ=sampleArray.SpikeZ, coordinates='Spacecraft', ScriptName=ScriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer
;
;
;  dummyRoeE = MAKE_ARRAY(N_ELEMENTS(sampleArray.thetaE), /DOUBLE, VALUE=0)
;
;  rotatedB = rotateMagDataTimeVarying(magx=sampleArray.Bscx, magy=sampleArray.Bscy, magz=sampleArray.Bscz, theta=sampleArray.thetaE, roe=dummyRoeE, phi=sampleArray.phiE)
;  
;  ici4timeseries, time=sampleArray.time, bx1=rotatedB(0,*), by1=rotatedB(1,*), bz1=rotatedB(2,*), bx2=sampleArray.IGRFgeix, by2=sampleArray.IGRFgeiy, bz2=sampleArray.IGRFgeiz, spikex=sampleArray.spikex, spikey=sampleArray.spikeY, spikeZ=sampleArray.SpikeZ, coordinates='Rotated', ScriptName=ScriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer


;  ;; Create bandpass filtered data
;bpLowFreq = 0.125
;bpHighFreq = 1.75
;
;
;print, "Bandstop filter between " + STRING(bpLowFreq) + " to " + STRING(bpHighFreq) + " Hz."
;Coeff = DIGITAL_FILTER(bpLowFreq*2.0/cadence, bpHighFreq*2/cadence, 100., 1000, /DOUBLE)
;
;if n_elements(Coeff) ge (sampleCounter) then begin
;  print, "ERROR: Number of sample (" + STRING(sampleCounter) + ")  is less than the length of the filter kernel (" +  STRING(N_ELEMENTS(Coeff), FORMAT='(I0)') + ") aborting."
;  stop
;endif
;
;BpBx = CONVOL(sampleArray.Bx, Coeff, /NAN, /EDGE_TRUNCATE)
;BpBy = CONVOL(sampleArray.By, Coeff, /NAN, /EDGE_TRUNCATE)
;BpBz = CONVOL(sampleArray.Bz, Coeff, /NAN, /EDGE_TRUNCATE)
;BB = CONVOL(sampleArray.B, Coeff, /NAN, /EDGE_TRUNCATE)
;  
;ici4timeseries, time=sampleArray.time, bx1=BpBz, by1=BpBy, bz1=BpBZ, spikey=sampleArray.spikeY, spikeZ=sampleArray.SpikeZ, coordinates='Sensor BP Filtered', ScriptName=ScriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer 
; 
  
  
  if ~keyword_set(disableslowplots) then begin
     
;    ;; Bandstop the spin frequency
;    
;    print, "Bandstop filter to suppress spin running at " + STRING(3.0) + " to " + STRING(3.5) + " Hz."
;    Coeff = DIGITAL_FILTER((3.5*2.0/cadence), (3.0*cadence), 100., 1000, /DOUBLE)
;
;    if n_elements(Coeff) ge (sampleCounter) then begin
;      print, "ERROR: Number of sample (" + STRING(sampleCounter) + ")  is less than the length of the filter kernel (" +  STRING(N_ELEMENTS(Coeff), FORMAT='(I0)') + ") aborting."
;      stop
;    endif
;
;    sampleArray.Bx = CONVOL(sampleArray.Bx, Coeff, /NAN, /EDGE_TRUNCATE)
;    sampleArray.By = CONVOL(sampleArray.By, Coeff, /NAN, /EDGE_TRUNCATE)
;    sampleArray.Bz = CONVOL(sampleArray.Bz, Coeff, /NAN, /EDGE_TRUNCATE)
;    sampleArray.B = CONVOL(sampleArray.B, Coeff, /NAN, /EDGE_TRUNCATE) 
;     
     
    ;; Welch's periodogram plotting
    
    wx = welch_periodogram(xin=sampleArray.Bx, tin=0, nfft=welchNFFT, cadence=cadence, overlap=welchOverlap, detrend=1)
    wy = welch_periodogram(xin=sampleArray.By, tin=0, nfft=welchNFFT, cadence=cadence, overlap=welchOverlap, detrend=1)
    wz = welch_periodogram(xin=sampleArray.Bz, tin=0, nfft=welchNFFT, cadence=cadence, overlap=welchOverlap, detrend=1)
    wm = welch_periodogram(xin=sampleArray.B, tin=0, nfft=welchNFFT, cadence=cadence, overlap=welchOverlap, detrend=1)
  
    freq = findgen(n_elements(wx))/n_elements(wx)*round(cadence/2)
  
    xmin = 10.0^floor(alog10(min( freq(1:-1) ) ) )
    xmax = 10.0^ceil(alog10(max( freq(1:-1) ) ) )
    ymin = 10.0^floor(alog10(min( [ wx, wy, wz, wm ] ) ) )
    ymax = 10.0^ceil(alog10(max( [ wx, wy, wz, wm ] ) ) )
    
    xrange = [ xmin, xmax ]
    yrange = [ ymin, ymax ]
    
    p2 = WINDOW(WINDOW_TITLE="ICI-4 Fluxgate Welch's Periodogram", DIMENSIONS=[11*140, 8.5*140], LAYOUT=[7,1,1], buffer=buffer)
  
    title = TEXT(0.5, 0.96, "ICI-4 Fluxgate Welch's Periodogram", FONT_SIZE=16, ALIGNMENT=0.5)
    
    position = [0.02,0.95,0.2,0.98]
  
    logo=read_image('UA-COLOUR.png')
    im1 = image(logo, POSITION=position, /OVERPLOT, buffer=buffer)
  
    position = [0.12,0.73,0.89,0.93]
    plot_wx = plot(freq, wx, XSHOWTEXT=0, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Sensor Coordinates!cBx (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,2], NAME='Bx', buffer=buffer)
  
    position = [0.12,0.52,0.89,0.72]
    plot_wx = plot(freq, wy, XSHOWTEXT=0, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Sensor Coordinates!cBy (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,3], NAME='By', buffer=buffer)
  
    position = [0.12,0.31,0.89,0.51]
    plot_wx = plot(freq, wz, XSHOWTEXT=0, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Sensor Coordinates!cBz (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,4], NAME='Bz', buffer=buffer)
  
  
    position = [0.12,0.10,0.89,0.30]
  
    ; Dummy plot to create the custom x axes
    plot_wx = plot(freq, wm, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Sensor Coordinates!c|B} (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,5], NAME='Bm', buffer=buffer)
  
    ; Footer string showing version and file information
    t7 = TEXT(0.03, 0.005, "Plotted using " + ScriptName + " " + VersionString + ' using ' + InputMagFileName, FONT_SIZE=PlotFontSize)
  
    ;; Get output filename
    outputFileName = InputMagFileName + FilePrefix + "WelchPeriodogram_Summary" + ".eps"
  
    print, "Saving ", outputFileName
    p2.save, outputFileName, RESOLUTION=PlotResolution
    if buffer then p2.close
  
  
    ;; Dynamic Spectra plotting
    
    dx = sonogram(xin=sampleArray.Bx, tin=0, nfft=DynSpectraNFFT, cadence=cadence, overlap=DynSpectraOverlap)
    dy = sonogram(xin=sampleArray.By, tin=0, nfft=DynSpectraNFFT, cadence=cadence, overlap=DynSpectraOverlap)
    dz = sonogram(xin=sampleArray.Bz, tin=0, nfft=DynSpectraNFFT, cadence=cadence, overlap=DynSpectraOverlap)
    dmag = sqrt(dx^2 + dy^2 + dz^2)
    
    freq = findgen(n_elements(dmag[0,*]))/n_elements(dmag[0,*])*round(cadence/2)
    
    numlevels = 20
    levels = 20*dindgen(numlevels+1)/numlevels-15
    levelnames=['-15','','','','','-10','','','','','-5','','','','','0','','','','','5']  
    
    p3 = WINDOW(WINDOW_TITLE="ICI-4 Fluxgate Dynamic Spectra", DIMENSIONS=[11*140, 8.5*140], LAYOUT=[6,1,1], buffer=buffer)
    
    title = TEXT(0.5, 0.96, "ICI-4 Fluxgate Dynamic Spectra", FONT_SIZE=16, ALIGNMENT=0.5)
    
    position = [0.02,0.95,0.2,0.98]
    
    logo=read_image('UA-COLOUR.png')
    im1 = image(logo, POSITION=position, /OVERPLOT, buffer=buffer)
    
    xrange = [min(SampleArray.time), max(sampleArray.time)]
    
    position = [0.12,0.78,0.89,0.94]
    c1 = CONTOUR(alog(dx), findgen(n_elements(dx[*,0])), freq, /FILL, RGB_TABLE=39, C_VALUE=levels, POSITION=position, /CURRENT, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Bx Sensor Coordinates!cFrequency (Hz)', YTICKFORMAT='(I0)', XMAJOR=6, XSTYLE=1, YSTYLE=1, LAYOUT=[6,1,2], buffer=buffer)
      
    position = [0.12,0.61,0.89,0.77]
    c1 = CONTOUR(alog(dy), findgen(n_elements(dy[*,0])), freq, /FILL, RGB_TABLE=39, C_VALUE=levels, POSITION=position, /CURRENT, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'By Sensor Coordinates!cFrequency (Hz)', YTICKFORMAT='(I0)', XMAJOR=6, XSTYLE=1, YSTYLE=1, LAYOUT=[6,1,3], buffer=buffer)
    
    position = [0.12,0.44,0.89,0.60]
    c1 = CONTOUR(alog(dz), findgen(n_elements(dz[*,0])), freq, /FILL, RGB_TABLE=39, C_VALUE=levels, POSITION=position, /CURRENT, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Bz Sensor Coordinates!cFrequency (Hz)', YTICKFORMAT='(I0)', XMAJOR=6, XSTYLE=1, YSTYLE=1, LAYOUT=[6,1,4], buffer=buffer)
    
    
    position = [0.12,0.27,0.89,0.43]
    
    ; Dummy plot to create the custom x axes
    yrange = [min(sampleArray.B), max(sampleArray.B)]
    c1 = CONTOUR(alog(dmag), findgen(n_elements(dmag[*,0])), freq, /FILL, RGB_TABLE=39, C_VALUE=levels, POSITION=position, /CURRENT, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = '|B|!cFrequency (Hz)', YTICKFORMAT='(I0)', XMAJOR=6, XSTYLE=1, YSTYLE=1, LAYOUT=[6,1,6], buffer=buffer)
      
    position = [0.12,0.10,0.89,0.26]
    plot_dummy = plot([0], [0], XRANGE=xrange, YRANGE=[0, 1], XSHOWTEXT=0, YSHOWTEXT=0, XMAJOR=0, XMINOR=0, YMAJOR=0, YMINOR=0, /NODATA, $
      YTICKFONT_SIZE=PlotFontSize, YTICKFORMAT='(I0)', POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[6,1,5], NAME='Dummy', buffer=buffer)
      
    yrange = [ 0 , 7 ]
    plot_spikex = plot(sampleArray.time, sampleArray.spikeX, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Removed Spikes', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='Spike X', buffer=buffer)
  
    plot_spikey = plot(sampleArray.time, sampleArray.spikeY, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Removed Spikes', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,7], NAME='Spike Y', buffer=buffer)
  
    plot_spikez = plot(sampleArray.time, sampleArray.spikeZ, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
      YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Removed Spikes', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,7], NAME='Spike Z', buffer=buffer)
      
    yReference = -0.02
  
    t1 = TEXT(0.03, 0.08, 'Time (s)', FONT_SIZE=PlotFontSize)
    xaxis_time = AXIS( 'X',  AXIS_RANGE=xrange, MAJOR=6, TICKFONT_SIZE=PlotFontSize, TARGET=plot_dummy, LOCATION=yreference, TICKLAYOUT=1)
    
    ; Footer string showing version and file information
    t7 = TEXT(0.03, 0.005, "Plotted using " + ScriptName + " " + VersionString + ' using ' + InputMagFileName, FONT_SIZE=PlotFontSize)
  
    ;; Get output filename
    outputFileName = InputMagFileName + FilePrefix + "DynamicSpectra_Summary" + ".eps"
  
    print, "Saving ", outputFileName
    p3.save, outputFileName, RESOLUTION=PlotResolution
    if buffer then p3.close
  endif
  
endif else begin
  print, "ERROR: No magnetic data specified - skipping"
endelse


if keyword_set(InputHKFileName) then begin
  print, "Parsing file: " + InputHKFileName + " for housekeeping data."
  
  ;; Parse the raw magnetic telemtry file
  worstCaseHKLength = file_lines(InputHKFileName)

  print, InputHKFileName + " contains approximately " + STRING(worstCaseHKLength) + " records or " + STRING(worstCaseHKLength/hkSamplingCadence) + " seconds of data."

  a = {hkSample, index:double(0), date:double(0), time:double(0), hk0:float(0), hk1:float(0), hk2:float(0), hk3:float(0), hk4:float(0), hk5:float(0), hk6:float(0), hk7:float(0), hk8:float(0), hk9:float(0), $
    spike0:float(0), spike1:float(0), spike2:float(0), spike3:float(0), spike4:float(0), spike5:float(0), spike6:float(0), spike7:float(0), spike8:float(0), spike9:float(0) }

  date = double(0)
  hk0 = float(0)
  hk1 = float(0)
  hk2 = float(0)
  hk3 = float(0)
  hk4 = float(0)
  hk5 = float(0)
  hk6 = float(0)
  hk7 = float(0)
  hk8 = float(0)
  hk9 = float(0)

  hkArray = replicate({hkSample}, worstCaseHKLength)

  hkCounter = ulong(0)

  line = ''

  openr, inputFileUnit, InputHKFileName, /GET_LUN
  readf, inputFileUnit, line ; Skip the header line


  while not eof(inputFileUnit) do begin
    readf, inputFileUnit, line
    
    
    hkString = STRMID(line,STREGEX(line, STRING(9b)))
    timeString = STRMID(line, 0, STREGEX(line, STRING(9b)))

    if STRLEN(timeString) lt 10 then begin
      hkArray[hkCounter].date = !VALUES.F_NAN
    endif else begin
      hkArray[hkCounter].date = JULDAY(loMonth,loDay,loYear,STRMID(timeString, 0, 2),STRMID(timeString, 3, 2),STRMID(timeString, 6, 9))
    endelse

    reads, strmid(hkString,0, 400), hk0, hk1, hk2, hk3, hk4, hk5, hk6, hk7, hk8, hk9

    hkArray[hkCounter].index = hkCounter
    hkArray[hkCounter].hk0 = hk0
    hkArray[hkCounter].hk1 = hk1
    hkArray[hkCounter].hk2 = hk2
    hkArray[hkCounter].hk3 = hk3
    hkArray[hkCounter].hk4 = hk4
    hkArray[hkCounter].hk5 = hk5
    hkArray[hkCounter].hk6 = hk6
    hkArray[hkCounter].hk7 = hk7
    hkArray[hkCounter].hk8 = hk8
    hkArray[hkCounter].hk9 = hk9


    hkCounter++

  endwhile

  close, inputFileUnit
  free_lun, inputFileUnit

  hkArray=hkArray[0:hkCounter-1]

  print, "Parsed " + string(hkCounter) + " housekeeping samples or " + string(hkCounter/hkSamplingCadence) + " seconds of data"
  
  ;; Deal with timestamp creation
  hkStartDate = MIN(hkArray.date, /NAN)
  hkStopDate = MAX(hkArray.date, /NAN)
  print, "Parsed housekeeping timestamps from " + STRING(hkStartDate, FORMAT=dateFormat  ) + " to " + STRING( hkStopDate, FORMAT=dateFormat ) + "."

  ;; Linear interpolation of timestamps to all samples using julian date.
  dateIndex = where(finite(hkArray.date))
  dateFit = ladfit(hkArray[dateIndex].index, hkArray[dateIndex].date, /DOUBLE)
  hkArray.date = dateFit[0] + dateFit[1]*hkArray.index

  ;; Date is often referenced to Seconds After Liftoff. Create this time index as well.

  loDate = JULDAY(loMonth,loDay,loYear,loHour,loMinute,loSecond)
  print, "Liftoff declared to be " + STRING(loDate, FORMAT=dateFormat) + "."

  loIndex = Value_Locate(hkArray.date, loDate)
  hkArray.time = 24D*60D*60D*( (dateFit(0) - hkArray(loIndex).date) + dateFit(1)*hkArray.index )

  hkStartTime = MIN(hkArray.time)
  hkStopTime = MAX(hkArray.time)
  print, "Data spans range of " + STRING(hkStartTime) + " to " + string(hkStopTime) + " seconds from liftoff."
  
  
  ;; Plot the sequence numbers to check for an instrument reset
  p99 = WINDOW(WINDOW_TITLE="ICI-4 Fluxgate Sequence", DIMENSIONS=[11*140, 8.5*140], buffer=buffer)
  plot_sequence = plot(hkArray.time, hkArray.hk9, YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Sequence Number', YTICKFORMAT='(I0)', XMAJOR=6, /CURRENT, COLOR='red', NAME='Sequence Number')
  outputFileName = InputHKFileName + FilePrefix  + "Sequence" + ".eps"
  print, "Saving ", outputFileName
  p99.save, outputFileName, RESOLUTION=PlotResolution
  if buffer then p99.close
  
  

  print, "Parsed housekeeping timestamps from " + STRING( hkStartTime ) + " to " + STRING( hkStopTime ) + "."

  if keyword_set( startSecond ) then begin
    if (startSecond gt hkStartTime) and (startSecond lt hkStopTime) then begin
      hkStartTime = startSecond
    endif else begin
      print, "StartSecond not possible in housekeeping data." + STRING(startSecond)
    endelse
  endif

  if keyword_set( stopSecond ) then begin
    if (stopSecond gt hkStartTime) and (stopSecond lt hkStopTime) then begin
      hkStopTime = stopSecond
    endif else begin
      print, "StopSecond not possible in housekeeping data." + STRING(stopSecond)
    endelse
  endif

  print, "Setting housekeeping timespan set to ", + STRING(hkStartTime) + " to " + STRING(hkStopTime)

  hkArray = hkArray[where(hkArray.time ge hkStartTime and hkArray.time le hkStopTime)]
  
  FilePrefix = "_" + STRING(hkStartTime,FORMAT='(I0)') + "s_to_" + STRING(hkStopTime,FORMAT='(I0)') + "s_"  + VersionString + "_"
  
  
  ; Look for and remove spikes (probably telemetery errors?)
  if ~keyword_set(disabledespike) then begin

    print, "Despiking housekeeping data with a window of " + string(hkDespikeInterval) + " and a threshold of " + string(hkDespikeThreshold) + " counts"

    temp = hkArray.hk0
    spikeIndexHK0 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK0 eq [-1] then spikeCountHK0 = 0 else spikeCountHK0 = N_ELEMENTS(spikeIndexHK0)
    hkArray[spikeIndexHK0].hk0 = !Values.D_NaN
    hkArray.spike0 = 1
    hkArray[spikeIndexHK0].spike0 = 2
    
    temp = hkArray.hk1
    spikeIndexHK1 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK1 eq [-1] then spikeCountHK1 = 0 else spikeCountHK1 = N_ELEMENTS(spikeIndexHK1)
    hkArray[spikeIndexHK1].hk1 = !Values.D_NaN
    hkArray.spike1 = 3
    hkArray[spikeIndexHK1].spike1 = 4
        
    temp = hkArray.hk2
    spikeIndexHK2 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK2 eq [-1] then spikeCountHK2 = 0 else spikeCountHK2 = N_ELEMENTS(spikeIndexHK2)
    hkArray[spikeIndexHK2].hk2 = !Values.D_NaN
    hkArray.spike2 = 5
    hkArray[spikeIndexHK2].spike2 = 6
    
    temp = hkArray.hk3
    spikeIndexHK3 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK3 eq [-1] then spikeCountHK3 = 0 else spikeCountHK3 = N_ELEMENTS(spikeIndexHK3)
    hkArray[spikeIndexHK3].hk3 = !Values.D_NaN
    hkArray.spike3 = 7
    hkArray[spikeIndexHK3].spike3 = 8
    
    temp = hkArray.hk4
    spikeIndexHK4 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK4 eq [-1] then spikeCountHK4 = 0 else spikeCountHK4 = N_ELEMENTS(spikeIndexHK4)
    hkArray[spikeIndexHK4].hk4 = !Values.D_NaN
    hkArray.spike4 = 9
    hkArray[spikeIndexHK4].spike4 = 10
    
    temp = hkArray.hk5
    spikeIndexHK5 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK5 eq [-1] then spikeCountHK5 = 0 else spikeCountHK5 = N_ELEMENTS(spikeIndexHK5)
    hkArray[spikeIndexHK5].hk5 = !Values.D_NaN
    hkArray.spike5 = 11
    hkArray[spikeIndexHK5].spike5 = 12
    
    temp = hkArray.hk6
    spikeIndexHK6 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK6 eq [-1] then spikeCountHK6 = 0 else spikeCountHK6 = N_ELEMENTS(spikeIndexHK6)
    hkArray[spikeIndexHK6].hk6 = !Values.D_NaN
    hkArray.spike6 = 13
    hkArray[spikeIndexHK6].spike6 = 14
    
    temp = hkArray.hk7
    spikeIndexHK7 = where(abs(temp-median(temp,hkDespikeInterval)) gt hkDespikeThreshold)
    if spikeIndexHK7 eq [-1] then spikeCountHK7 = 0 else spikeCountHK7 = N_ELEMENTS(spikeIndexHK7)
    hkArray[spikeIndexHK7].hk7 = !Values.D_NaN
    hkArray.spike7 = 15
    hkArray[spikeIndexHK7].spike7 = 16

    print, "Despiking found " + string(spikeCountHK0) + " in HK0. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK1) + " in HK1. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK2) + " in HK2. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK3) + " in HK3. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK4) + " in HK4. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK5) + " in HK5. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK6) + " in HK6. Replaced with NAN."
    print, "Despiking found " + string(spikeCountHK7) + " in HK7. Replaced with NAN."
  endif else begin
    print, "Skipping despike of Housekeeping data."
  endelse
  
  if ~keyword_set(disableinterpol) then begin
    hkArray[where( finite(hkArray.hk0, /NAN) )].hk0 = interpol( hkArray[where( finite(hkArray.hk0) )].hk0, hkArray[where( finite(hkArray.hk0) )].index, hkArray[where( finite(hkArray.hk0, /NAN) )].index)
    hkArray[where( finite(hkArray.hk1, /NAN) )].hk1 = interpol( hkArray[where( finite(hkArray.hk1) )].hk1, hkArray[where( finite(hkArray.hk1) )].index, hkArray[where( finite(hkArray.hk1, /NAN) )].index)
    hkArray[where( finite(hkArray.hk2, /NAN) )].hk2 = interpol( hkArray[where( finite(hkArray.hk2) )].hk2, hkArray[where( finite(hkArray.hk2) )].index, hkArray[where( finite(hkArray.hk2, /NAN) )].index)
    hkArray[where( finite(hkArray.hk3, /NAN) )].hk3 = interpol( hkArray[where( finite(hkArray.hk3) )].hk3, hkArray[where( finite(hkArray.hk3) )].index, hkArray[where( finite(hkArray.hk3, /NAN) )].index)
    hkArray[where( finite(hkArray.hk4, /NAN) )].hk4 = interpol( hkArray[where( finite(hkArray.hk4) )].hk4, hkArray[where( finite(hkArray.hk4) )].index, hkArray[where( finite(hkArray.hk4, /NAN) )].index)
    hkArray[where( finite(hkArray.hk5, /NAN) )].hk5 = interpol( hkArray[where( finite(hkArray.hk5) )].hk5, hkArray[where( finite(hkArray.hk5) )].index, hkArray[where( finite(hkArray.hk5, /NAN) )].index)
    hkArray[where( finite(hkArray.hk6, /NAN) )].hk6 = interpol( hkArray[where( finite(hkArray.hk6) )].hk6, hkArray[where( finite(hkArray.hk6) )].index, hkArray[where( finite(hkArray.hk6, /NAN) )].index)
    hkArray[where( finite(hkArray.hk7, /NAN) )].hk7 = interpol( hkArray[where( finite(hkArray.hk7) )].hk7, hkArray[where( finite(hkArray.hk7) )].index, hkArray[where( finite(hkArray.hk7, /NAN) )].index)
  endif else begin
    print, "Skipping interpolation of housekeeping data"
  endelse
  
  hkArray.hk0 = HK0CoeffC * (hkArray.hk0) + HK0CoeffD
  hkArray.hk1 = HK1CoeffC * (hkArray.hk1) + HK1CoeffD
  hkArray.hk2 = HK2CoeffC * (hkArray.hk2) + HK2CoeffD
  hkArray.hk3 = HK3CoeffC * (hkArray.hk3) + HK3CoeffD
  hkArray.hk4 = HK4CoeffC * (hkArray.hk4) + HK4CoeffD
  hkArray.hk5 = HK5CoeffC * (hkArray.hk5) + HK5CoeffD
  hkArray.hk6 = HK6CoeffC * (hkArray.hk6) + HK6CoeffD
  hkArray.hk7 = HK7CoeffC * (hkArray.hk7) + HK7CoeffD
  
  hk8Errors = where(hkArray.hk8 ne 0)
  hkArray.spike8 = 17
  hkArray[hk8Errors].spike8 = 18

  hkArray[0].hk9 = 0
  hkArray[1:*].hk9 = hkArray[1:-1].hk9 - hkArray[0:-2].hk9
  hk9NonSequence = where(hkArray.hk9 ne 1 and hkArray.hk9 ne -65536)
  hkArray.spike9 = 19
  hkArray[hk9NonSequence].spike9 = 20
  
  
  
  
  ;; Plot timeseries of housekeeping data
  p4 = WINDOW(WINDOW_TITLE="ICI-4 Fluxgate Housekeeping", DIMENSIONS=[11*140, 8.5*140], LAYOUT=[7,1,1], buffer=buffer)

  title = TEXT(0.5, 0.96, "ICI-4 Fluxgate Housekeeping", FONT_SIZE=16, ALIGNMENT=0.5)

  position = [0.02,0.95,0.2,0.98]

  logo=read_image('UA-COLOUR.png')
  im1 = image(logo, POSITION=position, /OVERPLOT, buffer=buffer)

  xrange = [min(hkArray.time), max(hkArray.time)]

  position = [0.12,0.78,0.89,0.94]
  yrange = [min(0), max(5.5)]
  
  plot_hk0 = plot(hkArray.time, hkArray.hk0, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Internal Voltages (V)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,2], NAME='+1V5 Supply')
  
  plot_hk3 = plot(hkArray.time, hkArray.hk3, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, /OVERPLOT, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Internal Voltages (V)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,2], NAME='+1V78 Ref')    
  
  plot_hk4 = plot(hkArray.time, hkArray.hk4, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, /OVERPLOT, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Internal Voltages (V)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,2], NAME='+5V0 Supply')
    
  position = [1.125,0.86]
  leg1 = LEGEND(TARGET=[plot_hk0,plot_hk3,plot_hk4], FONT_SIZE=PlotFontSize, VERTICAL_ALIGNMENT=0.5, POSITION=position, SHADOW=0, LINESTYLE='non', SAMPLE_WIDTH=0.05)
  leg1.rotate, 90    
    

  position = [0.12,0.61,0.89,0.77]
  yrange = [min([ -10, hkArray.hk1, hkArray.hk2, hkArray.hk6 ]), max([ hkArray.hk1, hkArray.hk2, hkArray.hk6, 40] )]
  
  plot_hk1 = plot(hkArray.time, hkArray.hk1, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Temperatures (!Eo!NC)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,3], NAME='Sensor Temp')

  plot_hk2 = plot(hkArray.time, hkArray.hk2, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, /OVERPLOT, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Temperatures (!Eo!NC)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,3], NAME='Reference Temp')

  plot_hk6 = plot(hkArray.time, hkArray.hk6, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, /OVERPLOT, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Temperatures (!Eo!NC)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,3], NAME='PCB Temp')
    
  position = [1.16,0.69]
  leg1 = LEGEND(TARGET=[plot_hk1,plot_hk2,plot_hk6], FONT_SIZE=PlotFontSize, VERTICAL_ALIGNMENT=0.5, POSITION=position, SHADOW=0, LINESTYLE='non', SAMPLE_WIDTH=0.05)
  leg1.rotate, 90

  position = [0.12,0.44,0.89,0.60]
  ;;yrange = [min([ 28-6, hkArray.hk5, 28+6] ), max( [ 28-6, hkArray.hk5, 28+6] )]
  yrange = [25, 30]
  plot_hk5 = plot(hkArray.time, hkArray.hk5, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Input Voltage (V)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,4], NAME='Input Voltage')
    
  position = [1.125,0.52]
  leg1 = LEGEND(TARGET=[plot_hk5], FONT_SIZE=PlotFontSize, VERTICAL_ALIGNMENT=0.5, POSITION=position, SHADOW=0, LINESTYLE='non', SAMPLE_WIDTH=0.05)
  leg1.rotate, 90
    
  
  position = [0.12,0.27,0.89,0.43]
  yrange = [min( [ 50, hkArray.hk7, 60 ] ), max([ 50, hkArray.hk7, 60 ] )]
  plot_hk7 = plot(hkArray.time, hkArray.hk7, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Input Current (mA)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,5], NAME='Input Current')

  position = [1.125,0.35]
  leg1 = LEGEND(TARGET=[plot_hk5], FONT_SIZE=PlotFontSize, VERTICAL_ALIGNMENT=0.5, POSITION=position, SHADOW=0, LINESTYLE='non', SAMPLE_WIDTH=0.05)
  leg1.rotate, 90


  position = [0.12,0.10,0.89,0.26]
  ; Dummy plot to create the custom x axes
  plot_dummy = plot([0], [0], XRANGE=xrange, YRANGE=[0, 1], XSHOWTEXT=0, YSHOWTEXT=0, XMAJOR=0, XMINOR=0, YMAJOR=0, YMINOR=0, /NODATA, $
    YTICKFONT_SIZE=PlotFontSize, YTICKFORMAT='(I0)', POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[7,1,6], NAME='Dummy', buffer=buffer)


  yrange = [ 0, 21 ]
  plot_spike0 = plot(hkArray.time, hkArray.spike0, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK0 Spikes')
  plot_spike1 = plot(hkArray.time, hkArray.spike1, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK1 Spikes')
  plot_spike2 = plot(hkArray.time, hkArray.spike2, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK2 Spikes')
  plot_spike3 = plot(hkArray.time, hkArray.spike3, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK3 Spikes')
  plot_spike4 = plot(hkArray.time, hkArray.spike4, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK4 Spikes')
  plot_spike5 = plot(hkArray.time, hkArray.spike5, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK5 Spikes')
  plot_spike6 = plot(hkArray.time, hkArray.spike6, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK6 Spikes')
  plot_spike7 = plot(hkArray.time, hkArray.spike7, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK7 Spikes')
  plot_spike8 = plot(hkArray.time, hkArray.spike8, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK8 Spikes')
  plot_spike9 = plot(hkArray.time, hkArray.spike9, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Check Values (counts)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='HK9 Spikes')
    
  yReference = -0.02

  t1 = TEXT(0.03, 0.08, 'Time (s)', FONT_SIZE=PlotFontSize)
  xaxis_time = AXIS( 'X',  AXIS_RANGE=xrange, MAJOR=6, TICKFONT_SIZE=PlotFontSize, TARGET=plot_dummy, LOCATION=yreference, TICKLAYOUT=1)


  ; Footer string showing version and file information
  t7 = TEXT(0.03, 0.005, "Plotted using " + ScriptName + " " + VersionString + ' using ' + InputHKFileName, FONT_SIZE=PlotFontSize)

  ;; Get output filename
  outputFileName = InputHKFileName + FilePrefix  + "Housekeeping" + ".eps"

  print, "Saving ", outputFileName
  p4.save, outputFileName, RESOLUTION=PlotResolution
  if buffer then p4.close
  
  
endif else begin
  print, "ERROR: No housekeeping data specified - skipping"
endelse
  
  


toc

end