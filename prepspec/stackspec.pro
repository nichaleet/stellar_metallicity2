pro stackspec,strin,strout,saveplot=saveplot,plotname=plotname
   ;Input is an array of structure
   ;Each structure should contain lambda,contdiv,contdivivar,dlam,sn
   ;Note that the lambda should be already in rest frame 
   ;Each lambda in each structure can have a head or tail of zero's to make the number of elements
   ;in each structure the same. Cannot have zeros spotting around inside the array.

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;find the total wl range (lambdafull)
   ndup = n_Elements(strin)
   lambdaarr = strin.lambda ;[npix,ndup]
   wzero = where(lambdaarr eq 0.,czero)
   if czero gt 0 then lambdaarr(wzero) = 1./0.
   minwl = min(lambdaarr,wminwl,/nan)
   maxwl = max(lambdaarr,wmaxwl,/nan)
   wminwl = array_indices(lambdaarr,wminwl)
   wmaxwl = array_indices(lambdaarr,wmaxwl)
   maxsn = max(strin.sn,ref)

   ;if overlapped, use the lambda of the spec with max sn as ref
   ;attach the extra lambdas to the main lambdas
   lambdafull = lambdaarr[*,ref]
   wzero = where(~finite(lambdafull),czero,complement=wgood)
   if czero gt 0 then lambdafull = lambdafull(wgood)

   if wmaxwl[1] ne ref then begin
       wadd = where(lambdaarr[*,wmaxwl[1]] gt max(lambdafull) and $
                    finite(lambdaarr[*,wmaxwl[1]]),cwadd)
       if cwadd gt 0 then lambdafull = [lambdafull,lambdaarr[min(wadd):wmaxwl[0],wmaxwl[1]]]
   endif
   if wminwl[1] ne ref then begin
       wadd = where(lambdaarr[*,wminwl[1]] lt min(lambdafull) and $
                    finite(lambdaarr[*,wminwl[1]]),cwadd)
       if cwadd gt 0 then lambdafull = [lambdaarr[wminwl[0]:max(wadd),wminwl[1]],lambdafull]
   endif
   ;get data in and interpolate
   npix = n_elements(lambdafull)
   strout={lambda:lambdafull,dlam:fltarr(npix),contdiv:fltarr(npix),contdivivar:fltarr(npix),sn:0.,exptime:0.}
   contdivarr=fltarr(npix,ndup)+1./0.
   contdivivararr = fltarr(npix,ndup)+1./0.
   dlamarr = fltarr(npix,ndup)+1./0.
   framearr = intarr(npix,ndup)
   for k=0,ndup-1 do begin
      strnow = strin[k]
      thislambda = strnow.lambda
      good = where(thislambda ne 0 and finite(thislambda))
      thislambda = thislambda(good)
      loc = value_locate(lambdafull,thislambda)
      if loc[0] eq -1 then stop,'this should not have happened'
      contdivarr[loc,k]=interpol(strnow.contdiv(good),thislambda,lambdafull[loc])
      contdivivararr[loc,k]=interpol(strnow.contdivivar(good),thislambda,lambdafull[loc])
      dlamarr[loc,k]=interpol(strnow.dlam(good),thislambda,lambdafull[loc])
      framearr[loc,k] +=1
   endfor
   ;weighted average
   nframearr = total(framearr,2)
   wnodup = where(nframearr eq 1, cwnodup)
   if cwnodup gt 0 then begin
       strout.dlam(wnodup) = total(dlamarr[wnodup,*],2,/nan)
       strout.contdiv(wnodup) = total(contdivarr[wnodup,*],2,/nan)
       strout.contdivivar(wnodup) = total(contdivivararr[wnodup,*],2,/nan)
   endif
   wdup = where(nframearr gt 1, cwdup)
   for ii=0,cwdup-1 do begin
       wnan = where(~finite(contdivivararr[wdup[ii],*]),cnan)
       dcontdiv = 1./sqrt(abs(contdivivararr[wdup[ii],*]))
       if cnan gt 0 then dcontdiv(wnan) = 1./0.
       err = abs(dcontdiv/contdivarr[wdup[ii],*])
       strout.dlam[wdup[ii]] = wmean(dlamarr[wdup[ii],*],abs(err*dlamarr[wdup[ii],*]),/nan)
       strout.contdiv[wdup[ii]]=wmean(contdivarr[wdup[ii],*],dcontdiv,error=contdivivar,/nan)
       strout.contdivivar[wdup[ii]]=contdivivar
   endfor
   ;signal to noise
   ;create mask
   readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
   contmask = bytarr(npix)+1
   lambdarest = lambdafull
   for i=0,n_elements(linestart)-1 do begin
      w = where(lambdarest ge linestart[i] and lambdarest le lineend[i], c)
      if c gt 0 then contmask[w] = 0
   endfor
   w = where(~finite(strout.contdiv) or ~finite(strout.contdivivar), c)
   if c gt 0 then contmask[w]=0
   ;calculate the deviation
   wcont = where(contmask eq 1)
   dev = abs(strout.contdiv[wcont] - 1.)
   avgdev = median(dev)
   w = where(dev lt 3.0*avgdev, c)
   if c gt 0 then strout.sn = 1.0/mean(dev[w])
   print,'signal to noise is ', strout.sn   

   if keyword_set(saveplot) then begin
      ;ploting
      colorlist = ['WT3','TAN3','BLK3','GRN3','BLU3','ORG3','RED3','PUR3','PBG3','YGB3','RYB3','TG3','WT6','TAN6','BLK6','GRN6','BLU6','ORG6','RED6','PUR6','PBG6','YGB6','RYB6','TG6']
      set_plot,'ps'
      !p.font=0
      if keyword_set(plotname) then psname='plot/'+plotname else psname='plot/stackedspec.eps'
      device, filename = psname,xsize = 30,ysize = 10, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color   
      plot,strout.lambda,strout.contdiv,yrange=[-2,4],xrange=[3500,7000]
      for i=0,ndup-1 do begin
         colornow = i mod n_elements(colorlist)
         oplot,strin[i].lambda,strin[i].contdiv,color=fsc_color(colorlist(colornow))
      endfor
      oplot,strout.lambda,strout.contdiv,color=fsc_color('red')
      device,/close
      set_plot,'x'
   endif
end
