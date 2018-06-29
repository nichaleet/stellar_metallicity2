pro make_telluric_curve_may2018
   obsfile = '/scr2/nichal/WMKO/LRIS/2018may_cl1604/reduce_GD140_400grating/spredux1d/elsr180512_0046r_GD140.spec'
   datafile = '/scr2/nichal/workspace4/prepspec_cl1604/telluric/gd140.dat.txt'

   readcol,obsfile,wl,flux,skyflux,fluxerr,xpix,ypix,comment='#'
   goodwl = where(wl gt 6500. and wl lt 10400.)
   obs = {wl:wl(goodwl),flux:flux(goodwl),skyflux:skyflux(goodwl),fluxerr:fluxerr(goodwl),xpix:xpix(goodwl),ypix:ypix(goodwl)}
   
   readcol,datafile,wl,mag,dwl
   mag = interpol(mag,wl,obs.wl)
   dwl = interpol(dwl,wl,obs.wl)
   flux = 10.^(-1.*(mag+48.59)/2.5)*10.^25 ;multiply by 10^25 to avoid small numbers
   data = {wl:obs.wl,mag:mag,flux:flux,dwl:dwl}

   tell = obs.flux/data.flux
   tellerr = tell*abs(obs.fluxerr/obs.flux)
   ;fit continuum to the curve
   tellstart = [6540.,6860.,7590.,8150,8900.,9300]
   tellend   = [6580.,6900.,7680.,8250,9200.,9700]
   mask = bytarr(n_elements(obs.wl))+1
   for i=0,n_elements(tellstart)-1 do begin
      w = where(obs.wl ge tellstart(i) and obs.wl le tellend(i),c)
      if c gt 0 then mask(w) = 0
   endfor
   irange = [6500,8800,9800,11000]
   won1 = where(mask and obs.wl ge irange[0] and obs.wl lt irange[1],cwon1)
   won2 = where(mask and obs.wl ge irange[1] and obs.wl lt irange[2],cwon1)
   won3 = where(mask and obs.wl ge irange[2] and obs.wl lt irange[3],cwon1)
   bkpt1 = slatec_splinefit(obs.wl(won1), tell(won1), coeff1, invvar=1/tellerr(won1)^2, bkspace=120, upper=3, lower=3, /silent) ;DEIMOS
   bkpt2 = slatec_splinefit(obs.wl(won2), tell(won2), coeff2, invvar=1/tellerr(won2)^2, bkspace=300, upper=3, lower=3, /silent) ;DEIMOS
   bkpt3 = slatec_splinefit(obs.wl(won3), tell(won3), coeff3, invvar=1/tellerr(won3)^2, bkspace=120, upper=3, lower=3, /silent) ;DEIMOS
   if bkpt1[0] eq -1 or bkpt2[0] eq -1 or bkpt3[0] eq -1 then message, 'Could not fit a spline to spsspec.'

   won1 = where(obs.wl ge irange[0] and obs.wl lt irange[1],cwon1)
   won2 = where(obs.wl ge irange[1] and obs.wl lt irange[2],cwon1)
   won3 = where(obs.wl ge irange[2] and obs.wl lt irange[3],cwon1)
   confit1 = slatec_bvalu(obs.wl(won1), bkpt1, coeff1)
   confit2 = slatec_bvalu(obs.wl(won2), bkpt2, coeff2)
   confit3 = slatec_bvalu(obs.wl(won3), bkpt3, coeff3)
   confit = [confit1,confit2,confit3]
   if n_Elements(confit) ne n_Elements(obs.wl) then stop,'range is wrong: stopped'
   !p.multi=[0,1,2]
   plot,obs.wl,tell
   oplot,obs.wl,confit,color=fsc_color('red')

   plot,obs.wl,tell/confit
   oplot,obs.wl,tellerr/confit
   tellstr = {lambda:obs.wl,spec:tell/confit,ivar:1./(tellerr/confit)^2}
   mwrfits,tellstr,'tellcor_gd140_may2018.fits',/create,/silent
   
   stop
end
