
RESTORE, FILENAME = 'bresidual.sav'
Bresidual2 = ici4_fft_bandstop(data=Bresidual, lowStopFreq=3, highStopFreq=7, cadence=100)
Bresidual3 = ici4_fft_bandstop(data=Bresidual2, lowStopFreq=1.0/20.0, highStopFreq=1.0/16.0, cadence=100)
plot(sampleArray[where(sampleArray.time gt 225 and sampleArray.time lt 325)].time,Bresidual3[where(sampleArray.time gt 225 and sampleArray.time lt 325)], 'blue')