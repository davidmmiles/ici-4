function spinHodogramCompare, b1a=b1a, b1b=b1b, b2a=b2a, b2b=b2b

len=max( [abs(b1a), abs(b1b), abs(b2a), abs(b2b)] )
len = 10*round(len/10)

p = plot(b1a,b1b, 'red')
p = plot(b2a,b2b, xaxis=[-len, len], yrange=[-len, len], 'blue', /overplot)

b1h2w = (max(b1b) - min(b1b))/(max(b1a) - min(b1a))
b2h2w = (max(b2b) - min(b2b))/(max(b2a) - min(b2a))
b1zero = [mean(b1a), mean(b1b)]
b2zero = [mean(b2a), mean(b2b)]


print, 'b1 has height/width ratio:' + string(b1h2w)
print, 'b2 has height/width ratio:' + string(b2h2w)
print, 'b1 has zeros of: [' + string(b1zero[0]) + ', ' + string(b1zero[1]) + ']'
print, 'b2 has zeros of: [' + string(b2zero[0]) + ', ' + string(b2zero[1]) + ']' 

return, [b1h2w, b2h2w, b1zero, b2zero]

end