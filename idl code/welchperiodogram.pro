pro welchPeriodogram, bx=bx, by=by, bz=bz, bmag=bmag, nfft=nfft, cadence=cadence, overlap=overlap, detrend=detrend, scriptname=scriptname, versionstring=vesionstring, InputMagFileName



;; Welch's periodogram plotting

wx = welch_periodogram(xin=bx, tin=0, nfft=nfft, cadence=cadence, overlap=overlap, detrend=1)
wy = welch_periodogram(xin=by, tin=0, nfft=nfft, cadence=cadence, overlap=overlap, detrend=1)
wz = welch_periodogram(xin=bz, tin=0, nfft=nfft, cadence=cadence, overlap=overlap, detrend=1)
wm = welch_periodogram(xin=bmag, tin=0, nfft=nfft, cadence=cadence, overlap=overlap, detrend=1)

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
  YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Bx (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,2], NAME='Bx', buffer=buffer)

position = [0.12,0.52,0.89,0.72]
plot_wx = plot(freq, wy, XSHOWTEXT=0, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
  YTICKFONT_SIZE=PlotFontSize, YTITLE = 'By (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='green', LAYOUT=[7,1,3], NAME='By', buffer=buffer)

position = [0.12,0.31,0.89,0.51]
plot_wx = plot(freq, wz, XSHOWTEXT=0, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
  YTICKFONT_SIZE=PlotFontSize, YTITLE = 'Bz (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='blue', LAYOUT=[7,1,4], NAME='Bz', buffer=buffer)


position = [0.12,0.10,0.89,0.30]

; Dummy plot to create the custom x axes
plot_wx = plot(freq, wm, XLOG=1, YLOG=1, XRANGE=xrange, YRANGE=yrange, $
  YTICKFONT_SIZE=PlotFontSize, YTITLE = '|B|} (PSD nT)', YTICKUNITS='scientific', POSITION=position, /CURRENT, COLOR='red', LAYOUT=[7,1,5], NAME='Bm', buffer=buffer)

; Footer string showing version and file information
t7 = TEXT(0.03, 0.005, "Plotted using " + ScriptName + " " + VersionString + ' using ' + InputMagFileName, FONT_SIZE=PlotFontSize)

;; Get output filename
outputFileName = InputMagFileName + FilePrefix + "WelchPeriodogram_Summary" + ".png"

print, "Saving ", outputFileName
p2.save, outputFileName, RESOLUTION=PlotResolution
if buffer then p2.close