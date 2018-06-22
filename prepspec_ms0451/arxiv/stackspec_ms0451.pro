pro stackspec_ms0451,retrieve=retrieve
   ;;;;;;;;;;;;;;;;;;;;;;
   ;getting info on individual spectra
   ;;;;;;;;;;;;;;;;;;;;;;
   dir = 'combined_uniq_ms0451_spline'
   headword = 'spline'
   strname = dir+'.info.fits'
   if keyword_set(retrieve) then begin
      path = '/scr2/nichal/workspace4/prepspec/'
      files = file_search(path+dir+'/*.fits', count=cspec)
      mass = fltarr(cspec)
      sn = fltarr(cspec) ;boxcar
      sngauss = fltarr(cspec)
      f814w = fltarr(cspec)
      redshift = fltarr(cspec)
      oiiew = fltarr(cspec)
      oiiewerr = fltarr(cspec)
      zsource = fltarr(cspec)
      for i=0,cspec-1 do begin
         spec = mrdfits(files[i],1)
         fits_open,files[i],fitsinfo
         cat = mrdfits(files[i],fitsinfo.nextend)
         free_lun,fitsinfo.unit
         print,string(i)+'  '+file_basename(files[i])
         print, 'SN:',spec.sn
         mass[i] = cat.mass_kcorrect
         sn[i] = spec.sn
         f814w[i] = cat.f814w_auto
         redshift[i] = cat.z
         zsource[i] = cat.zsource
         print,'zsource,zquality:',cat.zsource,cat.zquality
         oiiew[i] = oii_ew_new(spec.contdiv,spec.lambda,spec.contdivivar,cat.z,widtherr=widtherr,znew=znew)
         oiiewerr[i] = widtherr
         if znew ne redshift[i] then redshift[i] = znew
         if fitsinfo.nextend gt 4 then begin
            spec2 = mrdfits(files[i],3)
            sngauss[i] = spec2.sn
         endif
      endfor
      strsave = {mass:mass,sn:sn,sngauss:sngauss,f814w:f814w,redshift:redshift,$
             file:files,oiiew:oiiew,oiiewerr:oiiewerr}
      mwrfits,strsave,strname,/silent,/create
   endif
   info = mrdfits(strname,1)

   ;;;;;;;;;;;;;;;;;;
   ;stacking spectra
   ;;;;;;;;;;;;;;;;;
   massbounds = [9.5,10.2,10.7,11.2]
   snmax = 10.
   snmin = 3.
   nspec = n_elements(massbounds)-1
   npixtemp = 8195
   minoiiew = -5
   zmin = 0.52
   zmax = 0.56
   for i=0,nspec-1 do begin
       sel = where(info.mass gt massbounds[i] and info.mass le massbounds[i+1]$
                   and info.sn gt snmin and info.sn le snmax and info.oiiew gt minoiiew $
                   and info.redshift gt zmin and info.redshift lt zmax,csel)
       if csel lt 3 then stop,'number of individual spectra to be stacked is less than 3'
       template = {lambda:fltarr(npixtemp),dlam:fltarr(npixtemp),contdiv:fltarr(npixtemp),$
                   contdivivar:fltarr(npixtemp),sn:0.,exptime:0.}
       specall = replicate(template,csel)
       for j=0,csel-1 do begin
           spec = mrdfits(strtrim(info.file[sel[j]],2),1,/silent)
           ;change the lambda to restframe
           spec.lambda = spec.lambda/(1.+info.redshift(sel[j])) 
           strnow = template
           struct_assign,spec,strnow
           if total(strnow.lambda) eq 0. then stop, 'oh ooo'
           specall[j] = strnow
       endfor
       plot,specall[0].lambda,specall[0].contdiv,yrange=[0,4*csel],xrange=[3000,7000]
       for k=1,csel-1 do oplot,specall[k].lambda,specall[k].contdiv+4*k
       print,'doing mass', massbounds[i:i+1]
       print, 'stacking'+string(n_elements(specall))+'spectrum...'
       plotname= 'stackedspec_mass'+strjoin(strtrim(string(massbounds[i:i+1]),2),'_')+'.eps'
       stackspec,specall,specout,/saveplot,plotname=plotname
   endfor
end


