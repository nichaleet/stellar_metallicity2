pro make_telluric_curve_may2018
   obsfile = '/scr2/nichal/WMKO/LRIS/2018may_cl1604/reduce_GD140_400grating/spredux1d/elsr180512_0046r_GD140.spec'
   datafile = '/scr2/nichal/workspace4/prepspec_cl1604/telluric/gd140.dat.txt'

   readcol,obsfile,wl,flux,skyflux,fluxerr,xpix,ypix,comment='#'
   obs = {wl:wl,flux:flux,skyflux:skyflux,fluxerr:fluxerr,xpix:xpix,ypix:ypix}
   
   readcol,datafile,wl,mag,dwl
   mag = interpol(mag,wl,obs.wl)
   dwl = interpol(dwl,wl,obs.wl)
   flux = 10.^(-1.*(mag+48.59)/2.5)
   data = {wl:obs.wl,mag:mag,flux:flux,dwl:dwl}

   tell = data.flux/obs.flux
   tellerr = tell*obs.fluxerr/obs.flux
   ;fit continuum to the curve
   concoeff = poly_fit(obs.wl,tell,7,measure_errors=telerr,yfit=confit)
   !p.multi=[0,1,2]
   plot,obs.wl,tell
   oplot,obs.wl,confit,color=fsc_color('red')

   plot,obs.wl,tell/confit

   tell = {wl:obs.wl,tell:tell/confit,tellerr:tellerr/confit}
   mwrfits,'tellcor_gd140_may2018.fits',tell
   
   stop
end
