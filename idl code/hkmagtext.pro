;; Read in housekeeping magnetometers

pro hkmagtext

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

line = ''

openr, inputFileUnit3, hkMagFileName, /GET_LUN
readf, inputFileUnit3, line ; Skip the header line

stop

while not eof(inputFileUnit3) do begin
  readf, inputFileUnit3, line

  reads, strmid(line,0, 400), time, bx, by, bz
  
  hkMagArray.index = hkCounter
  hkMagArray.time = time
  hkMagArray.bx = bx
  hkMagArray.bx = by
  hkMagArray.bx = bz

  hkMagCounter++

endwhile

print, hkMagCounter

;; Cull to only the read samples
hkMagArray = hkMagArray[0:hkMagCounter-1]

end