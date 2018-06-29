
	load_colortable, /them

	ici4_read_magn, juls, secs, x, y, z

	ji = where( secs gt 100 and secs lt 500, jc )

	xfft = fft( x[ji]*hanning(jc), -1 ) & nf = n_elements(xfft)
	si = 700 & fi = 29999 & xfft[si:fi] = 1e-4 & xfft[nf-fi-1:nf-si-1] = 1e-4
	si = 0 & fi = 50 & xfft[si:fi] = 1e-4 & xfft[nf-fi-1:nf-si-1] = 1e-4

	yfft = fft( y[ji]*hanning(jc), -1 ) & nf = n_elements(yfft)
	si = 700 & fi = 29999 & yfft[si:fi] = 1e-4 & yfft[nf-fi-1:nf-si-1] = 1e-4
	si = 0 & fi = 50 & yfft[si:fi] = 1e-4 & yfft[nf-fi-1:nf-si-1] = 1e-4

	zfft = fft( z[ji]*hanning(jc), -1 ) & nf = n_elements(zfft)
	si = 700 & fi = 29999 & zfft[si:fi] = 1e-4 & zfft[nf-fi-1:nf-si-1] = 1e-4
	si = 0 & fi = 50 & zfft[si:fi] = 1e-4 & zfft[nf-fi-1:nf-si-1] = 1e-4

	hfft = fft( sqrt( x[ji]^2 + z[ji]^2 )*hanning(jc), -1 ) & nf = n_elements(hfft)
	si = 600 & fi = 29999 & hfft[si:fi] = 1e-4 & hfft[nf-fi-1:nf-si-1] = 1e-4
	si = 0 & fi = 50 & hfft[si:fi] = 1e-4 & hfft[nf-fi-1:nf-si-1] = 1e-4
	hfft[where(hfft gt 0.1)] = 1e-4

	xf = fft( xfft, 1 ) & yf = fft( yfft, 1 ) & zf = fft( zfft, 1 ) & hf = fft( hfft, 1 )

	ps_open, 'ici4_plot_magn.ps'

	trange = [150,450]
	set_format, /port, /sard
	cp

	plot, secs[ji], xf, chars=1, /xstyle, thick=2, pos=df(1,5,0,0,/no_tit), $
		yrange=[-3,3], /ystyle, xrange=trange, xtickname=replicate(' ', 30), ytitle='X [nT]'

	plot, secs[ji], yf, chars=1, /xstyle, thick=2, pos=df(/next,/no_tit), $
		yrange=[-3,3], /ystyle, xrange=trange, xtickname=replicate(' ', 30), ytitle='Y [nT]'

	plot, secs[ji], hf, chars=1, /xstyle, thick=2, pos=df(/next,/no_tit), $
		yrange=[-3,3], /ystyle, xrange=trange, xtickname=replicate(' ', 30), ytitle='H [nT]'

	plot, secs[ji], zf, chars=1, /xstyle, thick=2, pos=df(/next,/no_tit), $
		yrange=[-3,3], /ystyle, xrange=trange, xtickname=replicate(' ', 30), ytitle='Z [nT]'

	plot, secs[ji], sqrt( xf^2 + yf^2 + zf^2 ), chars=1, /xstyle, thick=2, pos=df(/next,/no_tit), $
		yrange=[0,3], /ystyle, xrange=trange, ytitle='|B| [nT]'

	ps_close, /no_f

	dt = 0.01
	fft_dur = 10.
	skip_dur = fft_dur/4
	fft_len = fft_dur/dt
	skip_len = skip_dur/dt

	wave_characteristics, z, x, y, juls, fft_len=fft_len, skip=skip_len, dt=dt, max_freq=4e3, sm_length=3, $
			freqs=freqs, pjuls=pjuls, $
			ellip=ellip, polarization=pol, $
			xpower=xpower, ypower=ypower, zpower=zpower, $
			polarized_power=polarized_power, e1=e1, e2=e2, e3=e3
	;if n_elements(idx) gt 1 then begin
	;	_x = median(_x,2) & y = mean(_y,dim=2) & _z = mean(_y,dim=2)
	;endif

; polar angle
phi = atan( e3[*,*,1], e3[*,*,0] )
theta = atan( sqrt(e3[*,*,0]^2+e3[*,*,1]^2), e3[*,*,2] )

ps_open, 'ici4_plot_magn_spec.ps'

	cp
	draw_image, alog10(ypower), (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[0,8], $
		ytitle='f [Hz]', xtickname=replicate(' ', 40), pos=df(1,3,0,0,/bar)
	plot_colorbar, panel_pos=df(/same,/bar), scale=[0,8], nlev=4, leg='X Power', colorstep=245
	draw_image, alog10(zpower), (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[0,8], $
		ytitle='f [Hz]', xtickname=replicate(' ', 40), pos=df(1,3,0,1,/bar)
	plot_colorbar, panel_pos=df(/same,/bar), scale=[0,8], nlev=4, leg='Y Power', colorstep=245
	draw_image, alog10(xpower), (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[0,8], $
		ytitle='f [Hz]', xtitle='Time after launch [s]', pos=df(1,3,0,2,/bar)
	plot_colorbar, panel_pos=df(/same,/bar), scale=[0,8], nlev=4, leg='Z Power', colorstep=245

	cp
	draw_image, phi[*,*]*!radeg, (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[-180,180], $
		ytitle='f [Hz]', xtickname=replicate(' ', 40), pos=df(1,3,0,0,/bar)
	plot_colorbar, panel_pos=df(/same,/bar), scale=[-180,180], nlev=4, leg=textoidl('\phi [\circ]'), $
		colorstep=245
	draw_image, theta[*,*]*!radeg, (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[0,180], $
		ytitle='f [Hz]', xtitle='Time after launch [s]', pos=df(1,3,0,1,/bar)
	plot_colorbar, panel_pos=df(/same,/bar), scale=[0,180], nlev=4, leg=textoidl('\theta [\circ]'), $
		colorstep=245

	;cp
	;draw_image, pol[*,*], (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[0,1], $
	;	ytitle='f [Hz]', xtickname=replicate(' ', 40), pos=df(1,3,0,0,/bar)
	;plot_colorbar, panel_pos=df(/same,/bar), scale=[0,1], nlev=4, leg=textoidl('pol'), $
	;	colorstep=245
	;draw_image, ellip[*,*], (pjuls-ici4_launch())*86400.d, freqs/1e3, chars=1, scale=[-1,1], $
	;	ytitle='f [Hz]', xtitle='Time after launch [s]', pos=df(1,3,0,1,/bar)
	;plot_colorbar, panel_pos=df(/same,/bar), scale=[-1,1], nlev=4, leg=textoidl('ellip'), $
	;	colorstep=245

	idx = [31,32,33]
	_x = fltarr( n_elements( phi[*,idx[0],0] ), n_elements(idx) )
	_y = fltarr( n_elements( phi[*,idx[0],0] ), n_elements(idx) )
	_z = fltarr( n_elements( phi[*,idx[0],0] ), n_elements(idx) )
	for f=0,n_elements(idx)-1 do begin
		_x[*,f] = cos( phi[*,idx[f],0] )*sin( theta[*,idx[f],0] )
		_y[*,f] = sin( phi[*,idx[f],0] )*sin( theta[*,idx[f],0] )
		_z[*,f] = cos( theta[*,idx[f],0] )
	endfor

	set_format, /gup
	cp
	plot, median(_x,dim=2), median(_y,dim=2), psym=1, pos=df(1,2,0,0,asp=1), xrange=[-1,1], yrange=[-1,1], $
		title='Spin at 3.2 Hz', xtitle='Z', ytitle='X'
	plots, median(_x), median(_y), psym=2, color=get_red()
	plot, median(_x,dim=2), median(_z,dim=2), psym=1, pos=df(1,2,0,1,asp=1), xrange=[-1,1], yrange=[-1,1], $
		xtitle='Z', ytitle='Y'
	plots, median(_x), median(_z), psym=2, color=get_red()

	idx = [3,4]
	_x = fltarr( n_elements( phi[*,idx[0],0] ), n_elements(idx) )
	_y = fltarr( n_elements( phi[*,idx[0],0] ), n_elements(idx) )
	_z = fltarr( n_elements( phi[*,idx[0],0] ), n_elements(idx) )
	for f=0,n_elements(idx)-1 do begin
		_x[*,f] = cos( phi[*,idx[f],0] )*sin( theta[*,idx[f],0] )
		_y[*,f] = sin( phi[*,idx[f],0] )*sin( theta[*,idx[f],0] )
		_z[*,f] = cos( theta[*,idx[f],0] )
	endfor

	cp
	plot, median(_x,dim=2), median(_y,dim=2), psym=1, pos=df(1,2,0,0,asp=1), xrange=[-1,1], yrange=[-1,1], $
		title='Coning at 0.3 Hz', xtitle='Z', ytitle='X'
	plots, median(_x), median(_y), psym=2, color=get_red()
	plot, median(_x,dim=2), median(_z,dim=2), psym=1, pos=df(1,2,0,1,asp=1), xrange=[-1,1], yrange=[-1,1], $
		xtitle='Z', ytitle='Y'
	plots, median(_x), median(_z), psym=2, color=get_red()

ps_close, /no_f
	end
