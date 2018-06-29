function spinHodogram, b1=b1, b2=b2, title=title

len=max( [abs(b1), abs(b2)] )
len = 10*round(len/10)

p = plot(b1,b2, title=title)

ratio = (max(b1) - min(b1))/(max(b2) - min(b2))
zeros = [mean(b1), mean(b2)]

print, Title + ' H/W ratio: ' + string(ratio) + ' zeroes: ' + string(zeros(0)) + ', ' + string(zeros(1))

return, [ratio, zeros]

end