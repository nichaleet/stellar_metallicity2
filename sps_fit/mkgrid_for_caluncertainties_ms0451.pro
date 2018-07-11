function make_spec,grid_feh,grid_age,grid_alpha=grid_alpha
common sps_spec, sps, spsz, spsage
common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
common get_sps_alpha, element

   ;setting
   vdisp = 300.

   ;setting of DEIMOS spectra
   redshift = 0.55
   dwl = 0.65 ;dispersion per pixel in angstrom
   iwave = 3800.
   nwl = 4096
   dlam = fltarr(nwl)+4.7 ;FWHM in Angstrom for Deimos600ZD grating 1" slit      
   obslambda = findgen(nwl)*dwl+iwave*(1.+redshift)
   restlambda = obslambda/(1.+redshift)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;create standard fit mask
   mask = bytarr(nwl)+1
   ;emission line
   readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, $
            format='D,D,A,X', /silent, comment='#'
   wem = where(linetype eq 'e', cemlines)
   for i=0,cemlines-1 do begin
       w = where(obslambda/(1d + redshift) ge linestart[wem[i]] and $
           obslambda/(1d + redshift) le lineend[wem[i]], c)
       if c gt 0 then mask[w] = 0
    endfor

   ;telluric region
   tellstart = [6864., 7591., 8938.]
   tellend = [6935., 7694., 9500.]
   for i=1,1 do begin
       w = where(obslambda ge tellstart[i] and obslambda le tellend[i], c)
       if c gt 0 then mask[w] = 0
    endfor
  ;Mg region
   if ~keyword_set(grid_alpha) then begin
      readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indstart,$
              indend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,$
              format='I,D,D,D,D,D,D,I,A',comment='#',/silent
      use_indices = ['Mg_b','Mg_2','Mg_1']
      for i=0,n_elements(use_indices)-1 do begin
         indnow = where(indname eq use_indices(i),cindnow)
         if cindnow eq 0 then stop
         w = where(obslambda/(1.+redshift) gt indstart(indnow[0]) and $
             obslambda/(1.+redshift) lt indend(indnow[0]),cw)
         if cw gt 0 then mask[w]=0
      endfor
   endif
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;wavelength range
   w = where(obslambda/(1.+redshift) lt 3650. or obslambda/(1.+redshift) gt 6000.,c)
   if c gt 0 then mask[w]=0
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nhalf = round(double(nwl) / 2.)
   mask[0:4] = 0
   mask[nhalf-4:nhalf+4] = 0
   mask[nwl-5:nwl-1] = 0
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   won = where(mask eq 1, con)
   wonfit = won

   xmp = obslambda[won]
   datalam = obslambda
   normalize = 0
   rest = 0

   ;loop over grid_feh and grid_age
   nfeh = n_elements(grid_feh)
   nage = n_elements(grid_age)
   if ~keyword_set(grid_alpha) then begin
      for iz=0,nfeh-1 do for ia=0,nage-1 do begin
         spsspec = get_sps_obs(xmp,[grid_feh(iz),grid_age(ia),vdisp/2.35,redshift])
         str = {lambda:xmp,spec:spsspec,feh:grid_feh(iz),age:grid_age(ia),alpha:0.0,vdisp:vdisp,redshift:redshift,element:element}
         if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage)
         strall[iz,ia] = str
         ;if iz mod 10 eq 0 and ia mod 10 eq 0 then print, iz, ia
      endfor
   endif else begin
      nalpha = n_elements(grid_alpha)
      element = ['Mg','O','Si','Ca','Ti']
      for  iz=0,nfeh-1 do for ia=0,nage-1 do for iap=0,nalpha-1 do begin
         spsspec = get_sps_alpha_obs(xmp,[grid_feh(iz),grid_age(ia),vdisp/2.35,redshift,grid_alpha(iap)])
         str = {lambda:xmp,spec:spsspec,feh:grid_feh(iz),age:grid_age(ia),alpha:grid_alpha(iap),$
                vdisp:vdisp,redshift:redshift,element:element}
         if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage,nalpha)
         strall[iz,ia,iap] = str
         ;if iz mod 10 eq 0 and ia mod 10 eq 0 then print, iz, ia
      endfor 
   endelse
   return,strall
end


pro mkgrid_for_caluncertainties_ms0451
   ;make spectra with these grids (Total 100 files at each fixed age)
   grid_feh = findgen(61)/50-1.003 ;range from [-1,0.197] dex 
   grid_age = findgen(100)/10.+0.1 ;range from [0.1,10] Gyr
   grid_alpha = findgen(41)/50.-0.4 ;range from [-0.4,0.4] dex
   for ia = 0,99 do begin
     for iap=0,40 do begin
      spec_arr = make_spec(grid_feh,grid_age(ia),grid_alpha=[grid_alpha(iap)])
      outfile = 'sspdegen/sspdegen_ms0451/age'+strtrim(string(ia,format='(I02)'),2)+'_sspdegen_gridspec_ms0451'+'.fits'
      if iap eq 0 then print, 'writing '+outfile
      if iap eq 0 then mwrfits,spec_arr,outfile,/silent,/create else mwrfits,spec_arr,outfile,/silent
      if spec_arr[0].alpha eq -999 then stop
     endfor
   endfor
end
