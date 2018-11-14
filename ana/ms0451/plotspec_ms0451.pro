pro plotspec_panel,science1,science2,epsname,lambdalim=lambdalim,ylim=ylim,nofit=nofit,smooth=smooth,annote=annote
nspec = 2
if ~keyword_set(lambdalim) then lambdalim = [3700,5500]
if ~keyword_Set(ylim) then ylim = [0,2.]
if ~keyword_set(nofit) then nofit = intarr(nspec)
if ~keyword_set(smooth) then smooth = intarr(nspec)

ysize = 3*nspec+2
xleft = 0.1
xright = 0.95
panelsize = 0.8/nspec
errmax=0.3
Delta = '!9'+string("104B)+'!x'
flambda = 'f!D!9'+string("154B)+'!x!n'

if keyword_set(idealmask) then readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indstart,indend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,format='I,D,D,D,D,D,D,I,A',comment='#',/silent

readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'

set_plot,'ps'
device, filename = epsname,xsize = 25,ysize = ysize, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.multi=[0,1,1]
!p.font = 0
sunsym=sunsymbol()
for i=0,nspec-1 do begin
   if i eq 0 then science = science1
   if i eq 1 then science = science2
   ;;;masking stuff;;;;;;;;;;;;;;;;; 
   ;read in ideal mask
   if keyword_set(idealmask) then begin
      mask = science.fitmask
      mask = bytarr(n_elements(mask))+1
      use_indices = ['Mg_b','Mg_2','Mg_1']
      for j=0,n_elements(use_indices)-1 do begin
           indnow = where(indname eq use_indices(j),cindnow)
           if cindnow eq 0 then stop
           w = where(science.lambda/(1.+science.zspec) gt indstart(indnow[0]) and science.lambda/(1.+science.zspec) lt indend(indnow[0]),cw)
           if cw gt 0 then mask[w]=0
      endfor
      w = where(science.lambda/(1.+science.zspec) lt 3650. or science.lambda/(1.+science.zspec) gt 5440.,c)
      if c gt 0 then mask[w]=0
      ;mark telluric regions
      n = n_elements(tellstart)
      for j=0,n-1 do begin
          w= where(science.lambda gt (tellstart)[j] and science.lambda lt (tellend)[j],c)
          if c gt 0 then mask[w]=0
      endfor

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      n = n_elements(science.lambda)
      nhalf = round(double(n) / 2.)
      mask[0:4] = 0
      mask[nhalf-4:nhalf+4] = 0
      mask[n-5:n-1] = 0
      science.fitmask = mask
   endif

   t = round(-1*ts_diff(science.fitmask, 1))
   wstart = where(t eq 1, cstart)+1
   wend = where(t eq -1, cend)
   if science.fitmask[0] eq 1 then begin
     if cstart eq 0 then begin
        wstart = 0
     endif else begin
        wstart = [0, wstart]
     endelse
     cstart += 1
  endif
   if science.fitmask[n_elements(t)-1] eq 1 then begin
      if cend eq 0 then begin
         wend = n_elements(t)-1
      endif else begin
        wend = [wend, n_elements(t)-1]
     endelse
      cend += 1
   endif
   if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if science.zfit ne 0 and finite(science.zfit) then z=science.zfit else z=science.zspec
   goodrange = where(science.lambda/(1d + z) gt lambdalim[0] and science.lambda/(1d + z) lt lambdalim[1])
   lambda = science.lambda(goodrange)/(1d + z)
   spec = science.contdiv(goodrange)/science.spscont(goodrange)
   ivar = science.contdivivar(goodrange)*(science.spscont(goodrange))^2
   model = science.spsspecfull;*science.spscont

   dspec = abs(deriv(lambda,spec))
   dylim = cgpercentiles(dspec,percentiles=0.97)
   dspecdiff = abs(spec/sqrt(ivar))
   errmax = cgpercentiles(dspecdiff(where(finite(dspecdiff))),percentiles=0.985)

;   badpix = where(dspecdiff gt errmax or spec lt 0 or dspec gt dylim and lambda lt 5400.,cbadpix)
;   if cbadpix gt 0 then spec(badpix) = 1./0.

   modeldiv = (spec-model(goodrange))/model(goodrange)
   modeldiverr = (1./sqrt(ivar))/model(goodrange)

   yhigh  = 0.95-panelsize*i
   ylow   = 0.95-panelsize*(i+1)
   ymid   = ylow+(yhigh-ylow)*0.35
   ymid = ylow
   ;;;;;;;;;;;;;;;;;;plot the spec (top panel);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if i eq 0 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=5, ystyle=5, color=fsc_color('black'),position=[xleft,ymid,xright,yhigh],/nodata
   if i gt 0 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=5, ystyle=5, color=fsc_color('black'),position=[xleft,ymid,xright,yhigh],/nodata,/noerase

   if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
      for j=0,cstart-1 do begin
         x = ([science.lambda[wstart[j]], science.lambda[wstart[j]], science.lambda[wend[j]], science.lambda[wend[j]]]/(1d + science.zspec) > lambdalim[0]) < lambdalim[1]
         y = [ylim[0], ylim[1], ylim[1], ylim[0]]
         polyfill, x, y, color=fsc_color('light cyan'),/data
      endfor
   endelse
   oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
   ;oplot, lambda,spec, color=fsc_color('black') 

   readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
   tellthick = [5, 2, 5, 2, 2,2]
   linewaves = [2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85]
   linenames = ['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', ' ', ' ', 'H8', 'CaK', ' ', 'CaH', ' ', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', ' ', 'Mgb', ' ', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]']
   linecolors = ['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue']
   ;mark telluric regions
   n = n_elements(tellstart)
   for j=0,n-1 do begin
      oplot, [(tellstart)[j], (tellend)[j]] / (1d + science.zspec), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(tellthick)[j]
   endfor
   ;mark and label lines
   n = n_elements(linewaves)
   for j=0,n-1 do begin
      if (linewaves)[j] le !X.CRANGE[0] or (linewaves)[j] ge !X.CRANGE[1] then continue
      oplot, [(linewaves)[j], (linewaves)[j]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((linecolors)[j])
      xyouts, (linewaves)[j]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (linenames)[j], orientation=90, alignment=1, color=fsc_color((linecolors)[j]),charsize=0.8
   endfor
   ;plot spec 
   plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), /noerase,position=[xleft,ymid,xright,yhigh],xtickformat='(A1)',yticks=5,ytickv=[0,0.5,1.,1.5],ytitle=flambda,charsize=1.3
  ; if i eq nspec-1 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!3rest wavelength (!sA!r!u !9o!n)!3', /noerase,position=[xleft,ymid,xright,yhigh],yticks=5,ytickv=[0,0.5,1.,1.5],ytitle='flux'
   ;plot model
  if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) and nofit[i] eq 0 then oplot, science.lambda/(1d + science.zspec), model, color=fsc_color('red')

  if science.logmstar lt 10. then sigfigmass = 2 else sigfigmass =3
  annote = 'log(M!D*!N/M'+sunsym+')='+sigfig(science.logmstar,sigfigmass)+'  S/N='+sigfig(science.snfit/sqrt(0.65),3)
  ;xyouts,4800,0.3,annote,charsize=1.3
  ;xyouts,5300,1.6,string(objnumarr[i])
   ;;;;;;;;;;;;;;;;;;;;;;;;;;plot model diviation(bottom panel);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;i

;  if i lt nspec-1 then plot,lambda,modeldiv,xrange=lambdalim,yrange=[-0.5,0.5],xstyle=1,ystyle=1,/noerase,position=[xleft,ylow,xright,ymid],ytitle=Delta+flambda+'/'+flambda,/nodata,xtickformat='(A1)',yticks=2,ytickv=[-0.3,0,0.3],charsize=1.3
;;  polyfill,[lambda,reverse(lambda)],[modeldiv+modeldiverr,reverse(modeldiv-modeldiverr)],color=fsc_Color('lightgray')
;;  if i eq nspec-1 then plot,lambda,modeldiv,xrange=lambdalim,yrange=[-0.5,0.5],xstyle=1,ystyle=1,/noerase,position=[xleft,ylow,xright,ymid],ytitle=Delta+flambda+'/'+flambda,/nodata,xtitle='!3Rest-frame wavelength (!sA!r!u !9o!n)!3',yticks=2,ytickv=[-0.3,0,0.3],charsize=1.3
;  if i eq nspec-1 then plot,lambda,modeldiv,xrange=lambdalim,yrange=[-0.5,0.5],xstyle=1,ystyle=1,/noerase,position=[xleft,ylow,xright,ymid],ytitle=Delta+flambda+'/'+flambda,/nodata,xtitle='Rest-frame wavelength (A)',yticks=2,ytickv=[-0.3,0,0.3],charsize=1.3
;  oplot,lambda,modeldiv
endfor
  axis,xaxis=0,xtitle='Rest-frame wavelength (A)',xrange=lambdalim,xstyle=1
  xyouts,0.625,0.06,string('!3!S!R!U!9 o!3'),/normal,charsize=1.2
device,/close
end

pro plotspec,science,epsname,lambdalim=lambdalim,ylim=ylim,nofit=nofit,smooth=smooth,annote=annote
  if ~keyword_set(lambdalim) then lambdalim = [3700,5500]
  if ~keyword_Set(ylim) then ylim = [0,2.]
  if ~keyword_set(nofit) then nofit = 0
  if ~keyword_set(smooth) then smooth = 0

  t = round(-1*ts_diff(science.fitmask, 1))
  wstart = where(t eq 1, cstart)+1
  wend = where(t eq -1, cend)
  if science.fitmask[0] eq 1 then begin
     if cstart eq 0 then begin
        wstart = 0
     endif else begin
        wstart = [0, wstart]
     endelse
     cstart += 1
  endif
  if science.fitmask[n_elements(t)-1] eq 1 then begin
     if cend eq 0 then begin
        wend = n_elements(t)-1
     endif else begin
        wend = [wend, n_elements(t)-1]
     endelse
     cend += 1
  endif
  if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
  set_plot,'ps'
  device, filename = epsname,xsize = 18,ysize = 6, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi=[0,1,1]
  !p.font = 0
  if science.zfit ne 0 and finite(science.zfit) then z=science.zfit else z=science.zspec
  goodrange = where(science.lambda/(1d + z) gt lambdalim[0] and science.lambda/(1d + z) lt lambdalim[1])
  lambda = science.lambda(goodrange)/(1d + z)
  spec = science.contdiv(goodrange)/science.spscont(goodrange)
  ivar = science.contdivivar(goodrange)
  if keyword_set(smooth) then begin
     mask = science.fitmask(goodrange)
     bad = where(spec gt 1.25 or spec lt 0.75 and lambda gt 5000. and lambda lt 5300.)
     ;bad = where(lambda gt 5000. and lambda lt 5500.)
     good = where(mask eq 1 and lambda gt 5000.)
     meandiv = stdev(spec(good))
     spec(Bad) = (randomn(seed,n_elements(bad))*meandiv)+1.
     ;spec(bad) = smooth(spec(bad),3,/nan)
  endif
  plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
  if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
     for i=0,cstart-1 do begin
        x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]]/(1d + science.zspec) > lambdalim[0]) < lambdalim[1]
        y = [ylim[0], ylim[1], ylim[1], ylim[0]]
        polyfill, x, y, color=fsc_color('light cyan')
     endfor
  endelse
  oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
  oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
  oplot, lambda,spec, color=fsc_color('black')

  readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
  tellthick = [5, 2, 5, 2, 2,2]
  linewaves = [2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85]
  linenames = ['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', ' ', ' ', 'H8', 'CaH', ' ', 'CaK', ' ', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', ' ', 'Mgb', ' ', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]']
  linecolors = ['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue']

  n = n_elements(tellstart)
  for i=0,n-1 do begin
;     oplot, [(tellstart)[i], (tellend)[i]] / (1d + science.zspec), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(tellthick)[i]
  endfor
  n = n_elements(linewaves)
  for i=0,n-1 do begin
     if (linewaves)[i] le !X.CRANGE[0] or (linewaves)[i] ge !X.CRANGE[1] then continue
     oplot, [(linewaves)[i], (linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((linecolors)[i])
     xyouts, (linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (linenames)[i], orientation=90, alignment=1, color=fsc_color((linecolors)[i])
  endfor

  plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!3rest wavelength (!sA!r!u !9o!n)!3', ytitle='!3flux (normalized)!3', /nodata, /noerase

  if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) and nofit eq 0 then $
       oplot, science.lambda/(1d + science.zspec), science.spsspecfull, color=fsc_color('red')

  if keyword_set(annote) then xyouts,0.15,0.26,annote,/normal

  device,/close
end

pro zoommgb,sci1,sci2,epsname

  science=sci2
  if science.zfit ne 0 and finite(science.zfit) then z=science.zfit else z=science.zspec
  goodrange = where(science.lambda/(1d + z) gt 5250. and science.lambda/(1d + z) lt 5310)
  lambda = science.lambda(goodrange)/(1d + z)
  spec = science.contdiv(goodrange);/science.spscont(goodrange)
  ivar = science.contdivivar(goodrange)

  set_plot,'ps'
  device, filename = epsname,xsize = 10,ysize = 10, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi=[0,1,1]
    plot,lambda,spec
    oplot,sci1.lambda/(1.+sci1.zspec),sci1.spsspec,color=fsc_color('red')
    oplot,sci2.lambda/(1.+sci2.zspec),sci2.spsspec,color=fsc_color('navy')
  device,/close

end

pro plotspec_ms0451
  sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit02.fits.gz',1)
  sci1 = sci[65]
  print, sci1.objname
;  plotspec,sci,'p40132_feh_age.eps'
  sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
  sci2 = sci[33]
  print, sci2.objname
  wnomask = where(sci2.lambda/(1.+sci2.zfit) gt 3800. and sci2.lambda/(1.+sci2.zfit) lt 4500.)
  sci2.fitmask(wnomask) = 1
  plotspec,sci2,'p40132_copyi3.eps'
;  plotspec_panel,sci1,sci2,'p40132_fitpanel.eps'
  ;zoom Mgb band
;  zoommgb,sci1,sci2,'p40132_zoommgb.eps'
end
