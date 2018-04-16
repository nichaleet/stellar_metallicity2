pro stackspec_ms0451,retrieve=retrieve
   ;;;;;;;;;;;;;;;;;;;;;;;
   ;getting info on individual spectra
   ;;;;;;;;;;;;;;;;;;;;;;
   dir = 'combined_uniq_ms0451_spline'
   headword = 'spline'
   strname = dir+'.info.fits'
   if keyword_set(retrieve) then begin
      path = '/scr2/nichal/workspace4/prepspec/'
      files = file_search(path+dir[idir]+'/*.fits', count=cspec)
      mass = fltarr(cspec)
      sn = fltarr(cspec) ;boxcar
      sngauss = fltarr(cspec)
      f814w = fltarr(cspec)
      redshift = fltarr(cspec)
      for i=0,cspec-1 do begin
         spec = mrdfits(files[i],1)
         fits_open,files[i],fitsinfo
         cat = mrdfits(files[i],fitsinfo.nextend)
         free_lun,fitsinfo.unit
         mass[i] = cat.mass_kcorrect
         sn[i] = spec.sn
         f814w[i] = cat.f814w_auto
         redshift[i] = cat.z
         if fitsinfo.nextend gt 4 then begin
            spec2 = mrdfits(files[i],3)
            sngauss[i] = spec2.sn
         endif
      endfor
      str = {mass:mass,sn:sn,sngauss:sngauss,f814w:f814w,redshift:redshift,$
             file:files}
      mwrfits,strsave,strname[idir],/silent,/create
   endif
   info = mrdfits(strname,1)

   ;;;;;;;;;;;;;;;;;;
   ;stacking spectra
   ;;;;;;;;;;;;;;;;;
   massbounds = [9,9.7,10.2,10.7,11.2]
   snmax = 6.
   snmin = 0.1
   nspec = n_elements(massbounds)-1
   for i=0,nspec-1 do begin
       sel = where(info.mass gt massbounds[i] and info.mass le massbounds[i+1]$
                   and info.sn gt snmin and info.sn le snmax,csel)
       if csel lt 3 then stop,'number of individual spectra to be stacked is less than 3'
       
   endfor
end
