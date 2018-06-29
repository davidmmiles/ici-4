pro ici4timeseries, time=time, bx1=bx1, by1=by1, bz1=bz1, bx2=bx2, by2=by2, bz2=bz2, spikex=spikex, spikey=spikey, spikez=spikez, coordinates=coordinates, scriptName=scriptName, VersionString=VersionString, InputMagFileName=InputMagFileName, FilePrefix=FilePrefix, buffer=buffer

  numSamples = N_ELEMENTS(time)

  if ~keyword_set(bx2) then bx2 = MAKE_ARRAY(numSamples, /DOUBLE, VALUE=!Values.D_NaN)
  if ~keyword_set(by2) then by2 = MAKE_ARRAY(numSamples, /DOUBLE, VALUE=!Values.D_NaN)
  if ~keyword_set(bz2) then bz2 = MAKE_ARRAY(numSamples, /DOUBLE, VALUE=!Values.D_NaN)
  
  bx1 = reform(bx1)
  by1 = reform(by1)
  bz1 = reform(bz1)
  bx2 = reform(bx2)
  by2 = reform(by2)
  bz2 = reform(bz2)

  B1 = sqrt(Bx1^2 + By1^2 + Bz1^2)
  B2 = sqrt(Bx2^2 + By2^2 + Bz2^2)

  p1 = WINDOW(WINDOW_TITLE="ICI-4 Fluxgate Summary Plot" + coordinates, DIMENSIONS=[11*140, 8.5*140], LAYOUT=[7,1,1], buffer=buffer)

  title = TEXT(0.5, 0.96, "ICI-4 Fluxgate Summary Plot " + coordinates, FONT_SIZE=16, ALIGNMENT=0.5)

  position = [0.02,0.95,0.2,0.98]

  logo=read_image('UA-COLOUR.png')
  im1 = image(logo, POSITION=position, /OVERPLOT, buffer=buffer)

  xrange = [min(time), max(time)]

  maxComponent = 65000
  maxMagnitude = 65000

  position = [0.12,0.78,0.89,0.94]
yrange = [-maxComponent, maxComponent];  yrange = [ max( [ min([Bx1, Bx2]), -maxComponent ] ), min( [ max([Bx1, Bx2]), maxComponent ] ) ]
  plot_bx1 = plot(time, Bx1, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!cBx (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,2], NAME='Bx', buffer=buffer)
  plot_bx2 = plot(time, Bx2, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!c|B| (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[7,1,6], NAME='IGRFBx', buffer=buffer)

  position = [0.12,0.61,0.89,0.77]
yrange = [-maxComponent, maxComponent];  yrange = [ max( [ min([By1, By2]), -maxComponent ] ), min( [ max([By1, By2]), maxComponent ] ) ]
  plot_by1= plot(time, By1, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!cBy (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,3], NAME='By', buffer=buffer)
  plot_by2 = plot(time, By2, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!c|B| (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[7,1,6], NAME='IGRFBy', buffer=buffer)

  position = [0.12,0.44,0.89,0.60]
yrange = [-maxComponent, maxComponent];  yrange = [ max( [ min([Bz1,Bz2]), -maxComponent ] ), min( [ max([Bz1, Bz2]), maxComponent ] ) ]
  plot_bz1 = plot(time, Bz1, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!cBz (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,4], NAME='Bz', buffer=buffer)
  plot_bz2 = plot(time, Bz2, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!c|B| (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[7,1,6], NAME='IGRFBz', buffer=buffer)

  position = [0.12,0.27,0.89,0.43]
  yrange = [ max( [ min([B1, B2]), -maxMagnitude ] ), min( [ max([B1, B2]), maxMagnitude] ) ]
  plot_b1 = plot(time, B1, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!c|B| (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='black', LAYOUT=[7,1,6], NAME='B', buffer=buffer)
  plot_b2 = plot(time, B2, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Sensor Coordinates!c|B| (nT)', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[7,1,6], NAME='IGRFB', buffer=buffer)

  position = [0.12,0.10,0.89,0.26]
  ; Dummy plot to create the custom x axes
  plot_dummy = plot([0], [0], XRANGE=xrange, YRANGE=[0, 1], XSHOWTEXT=0, YSHOWTEXT=0, XMAJOR=0, XMINOR=0, YMAJOR=0, YMINOR=0, /NODATA, $
    YTICKFONT_SIZE=8, YTICKFORMAT='(I0)', POSITION=position, /CURRENT, COLOR='pink', LAYOUT=[7,1,5], NAME='Dummy', buffer=buffer)


  yrange = [ 0 , 7 ]
  plot_spikex = plot(time, spikeX, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Removed Spikes', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,7], NAME='Spike X', buffer=buffer)

  plot_spikey = plot(time, spikeY, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Removed Spikes', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,7], NAME='Spike Y', buffer=buffer)

  plot_spikez = plot(time, spikeZ, XRANGE=xrange, YRANGE=yrange, XSHOWTEXT=0, $
    YTICKFONT_SIZE=8, YTITLE = 'Removed Spikes', YTICKFORMAT='(I0)', XMAJOR=6, POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,7], NAME='Spike Z', buffer=buffer)

  yReference = -0.02

  t1 = TEXT(0.03, 0.08, 'Time (s)', FONT_SIZE=8)
  xaxis_time = AXIS( 'X',  AXIS_RANGE=xrange, MAJOR=6, TICKFONT_SIZE=8, TARGET=plot_dummy, LOCATION=yreference, TICKLAYOUT=1)


  ; Footer string showing version and file information
  t7 = TEXT(0.03, 0.005, "Plotted using " + ScriptName + " " + VersionString + ' using ' + InputMagFileName, FONT_SIZE=8)

  ;; Get output filename
  outputFileName = InputMagFileName + FilePrefix + coordinates + "Timeseries_Summary" + ".eps"

  print, "Saving ", outputFileName
  p1.save, outputFileName, RESOLUTION=PlotResolution
  if buffer then p1.close
  
  end