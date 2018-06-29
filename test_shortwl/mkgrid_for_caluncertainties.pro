function make_spec,grid_feh,grid_age,full=full,short=short
common sps_spec, sps, spsz, spsage
common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
   ;setting
   vdisp = 300.

   ;get the obs setting    
   sample = mrdfits('/scr2/nichal/workspace4/test_shortwl/mockdata/mockdeimos_cl1604_zmetm0d05_1gyr_sn10_1.fits',1,/silent)
   if keyword_set(short) then begin
      reallambda = sample.lambda
      dlam = sample.dlam/2.35
   endif else if keyword_set(full) then begin
      reallambda = sample.lambdafull
      dlam = sample.dlamfull/2.35
   endif 
   znow = sample.z
   nwl = n_Elements(reallambda)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;create standard fit mask
   mask = bytarr(nwl)+1
   ;emission line
   readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, $
            format='D,D,A,X', /silent, comment='#'
   wem = where(linetype eq 'e', cemlines)
   for i=0,cemlines-1 do begin
       w = where(reallambda/(1d + sample.z) ge linestart[wem[i]] and $
           reallambda/(1d + sample.z) le lineend[wem[i]], c)
       if c gt 0 then mask[w] = 0
    endfor

   ;telluric region
   tellstart = [6864., 7591., 8938.]
   tellend = [6935., 7694., 9500.]
   for i=1,1 do begin
       w = where(reallambda ge tellstart[i] and reallambda le tellend[i], c)
       if c gt 0 then mask[w] = 0
    endfor
  ;Mg region
   readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indstart,$
           indend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,$
           format='I,D,D,D,D,D,D,I,A',comment='#',/silent
   use_indices = ['Mg_b','Mg_2','Mg_1']
   for i=0,n_elements(use_indices)-1 do begin
      indnow = where(indname eq use_indices(i),cindnow)
      if cindnow eq 0 then stop
      w = where(reallambda/(1.+sample.z) gt indstart(indnow[0]) and $
          reallambda/(1.+sample.z) lt indend(indnow[0]),cw)
      if cw gt 0 then mask[w]=0
   endfor
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;wavelength range
   w = where(reallambda/(1.+sample.z) lt 3650. or reallambda/(1.+sample.z) gt 6000.,c)
   if c gt 0 then mask[w]=0
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nhalf = round(double(nwl) / 2.)
   mask[0:4] = 0
   mask[nhalf-4:nhalf+4] = 0
   mask[nwl-5:nwl-1] = 0
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   won = where(mask eq 1, con)
   wonfit = won

   xmp = reallambda[won]
   datalam = reallambda
   normalize = 0
   rest = 0

   ;loop over grid_feh and grid_age
   nfeh = n_elements(grid_feh)
   nage = n_elements(grid_age)
   for iz=0,nfeh-1 do for ia=0,nage-1 do begin
           spsspec = get_sps_obs(xmp,[grid_feh(iz),grid_age(ia),vdisp/2.35,znow])
           str = {lambda:xmp,spec:spsspec,feh:grid_feh(iz),age:grid_age(ia),vdisp:vdisp,redshift:znow}
           if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage)
           strall[iz,ia] = str
           ;if iz mod 10 eq 0 and ia mod 10 eq 0 then print, iz, ia
   endfor
   return,strall
end


pro mkgrid_for_caluncertainties
   ;make spectra with these grids (Total 100 files at each fixed age)
   grid_feh = findgen(120)/99.4-1. ;range from [-1,0.197] dex ;for new_grid_spec.fits
   grid_age = findgen(100)/10.+0.1 ;range from [0.1,10] Gyr
   for ia = 0,99 do begin
      spec_arr = make_spec(grid_feh,grid_age(ia),/short)
      outfile = 'sspdegen/sspdegen_short/age'+strtrim(string(ia,format='(I02)'),2)+'_sspdegen_gridspec_shortwl'+'.fits'
      print, 'writing '+outfile
      mwrfits,spec_arr,outfile,/silent,/create
      spec_arr = make_spec(grid_feh,grid_age(ia),/full)
      outfile =  'sspdegen/sspdegen_full/age'+strtrim(string(ia,format='(I02)'),2)+'_sspdegen_gridspec_fullwl'+'.fits'
      print, 'writing '+outfile
      mwrfits,spec_arr,outfile,/silent,/create
   endfor
end
