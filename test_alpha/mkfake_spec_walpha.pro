pro mkfake_spec_walpha,notput=notput
   ;The fake data would be similar to that of Deimos data of ms0451
   zmet = [-0.5,-0.2,0.1]  
   age = [3.]
   sn  = [5,10,15,20,30] ;per angstrom at observed wl
   alpha = [-0.25,-0.15,0.,0.15,0.25] 
   element=['Mg','O','Si','Ca','Ti']
   nelements = n_elements(element)
   nspec = 5
   redshift = 0.55
   vdisp = 300. ;km/s in FWHM
   nfinalspec = n_Elements(zmet)*n_Elements(age)*n_Elements(sn)*n_Elements(alpha)*nspec
   print,'Making '+strtrim(string(nfinalspec),2)+' spectra.....'

   ;setting of DEIMOS spectra
   dwl = 0.65 ;dispersion per pixel in angstrom
   iwave = 3800.
   nwl = 4096
   dlam = fltarr(nwl)+4.7 ;FWHM in Angstrom for Deimos600ZD grating 1" slit      
   lambdaobs = findgen(nwl)*dwl+iwave*(1.+redshift)
   lambdarest = lambdaobs/(1.+redshift)
   snperpix = sn*sqrt(dwl) 

   ;setup for output names
   dir = '/scr2/nichal/workspace4/test_alpha/mockdata_alpha/'
   prefix = 'mockms0451'
   zmetname = strarr(n_Elements(zmet))
   for i=0,n_Elements(zmet)-1 do begin
      tmpz = strsplit(strtrim(string(zmet(i)),2),'.',/extract)
      mop = (strmid(tmpz[0],0,1) eq '-') ? 'm' : 'p'
      fdg = (strmid(tmpz[0],0,1) eq '-') ? strmid(tmpz[0],1,1) : strmid(tmpz[0],0,1) 
      sdg = strmid(tmpz[1],0,2)
      zmetname[i] = '_zmet'+mop+fdg+'d'+sdg
   endfor
   agename = '_'+strtrim(string(age,format='(I)'),2)+'gyr' 
   snname =  '_sn'+strtrim(string(sn,format='(I)'),2)
   alphaname = strarr(n_Elements(alpha))
   for i=0,n_Elements(alpha)-1 do begin
      tmpp = strsplit(strtrim(string(alpha(i)),2),'.',/extract)
      mop = (strmid(tmpp[0],0,1) eq '-') ? 'm' : 'p'
      fdg = (strmid(tmpp[0],0,1) eq '-') ? strmid(tmpp[0],1,1) : strmid(tmpp[0],0,1)
      sdg = strmid(tmpp[1],0,2)
      alphaname[i] = '_alpha'+mop+fdg+'d'+sdg
   endfor
   idname =  '_'+strtrim(string(indgen(nspec)+1,format='(I)'),2)

   ;constant
   clight = 299792.458
   ;throughput
   readcol,'thr_g600_65_gg455.asc',wltput,tput,comment='#'
   tput = interpol(tput,wltput,lambdaobs)
   tput = tput/median(tput)

   ;setting continuum mask
   readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
   contmask = bytarr(nwl)+1
   for j=0,n_elements(linestart)-1 do begin
      w = where(lambdarest ge linestart[j] and lambdarest le lineend[j], c)
      if c gt 0 then contmask[w] = 0
   endfor
   tellstart = [6864., 7591., 8938.]
   tellend = [6935., 7694., 9500.]
   for j=0,n_elements(tellstart)-1 do begin
       w = where(lambdaobs ge tellstart[j] and lambdaobs le tellend[j], c)
       if c gt 0 then begin
           if j eq 0 then contmask[w] = 0
       endif
   endfor
   contmask[0:2] = 0
   contmask[nwl-3:nwl-1] = 0
   won = where(contmask eq 1, complement=woff, con)

   ;create spectrum!!!
   for iz=0,n_elements(zmet)-1 do for ia=0,n_elements(age)-1 do for iap=0,n_elements(alpha)-1 do begin
      spsstruct = sps_interp(zmet[iz], age[ia])
      spslambda = spsstruct.lambda
      spsspec = spsstruct.spec
      w = where(spslambda gt 3000 and spslambda lt 6000, c)
      if c lt 25 then message, 'Not enough pixels.'
      spslambda = spslambda[w]
      spsspec = spsspec[w]
      spsspec = spsspec*clight/spslambda^2      ;change fnu(Lsun/Hz) to flambda 
      ;add response function
      abund = fltarr(nelements)+alpha(iap)
      spsspec = add_response(spslambda,spsspec,[zmet[iz],age[ia]],element,abund)
      ;normalize to around 1
      spsmedian = median(spsspec)
      spsspec = spsspec/spsmedian     
      ;smooth to vdisp
      spsspec = smooth_gauss_wrapper(spslambda, spsspec, spslambda, vdisp/clight/2.35*spslambda)

      ;smooth to deimos resolution and get observed wavelength array
      spsspec = smooth_gauss_wrapper(spslambda*(1.+redshift), spsspec, lambdaobs, dlam/2.35)
      ;apply through put
      if ~keyword_set(notput) then spsspec = spsspec*tput

      ;add SN
      for isn=0,n_Elements(sn)-1 do for ispec=0,nspec-1 do begin
         sigmasq = spsspec^2/snperpix(isn)^2
         spec = spsspec+randomn(seed,nwl)*sqrt(sigmasq)
         ivar = 1./sigmasq
         ;continuum normalizei
         bkspace =  fix(200./median(abs(ts_diff(lambdaobs,1)))) ;approximately every 200 A
         bkpt = slatec_splinefit(lambdaobs[won], spec[won], coeff, invvar=ivar[won], $
                bkspace=bkspace, upper=3, lower=3, /silent)
         if bkpt[0] eq -1 then stop,'cannot do spline continnum fit'
         cont = slatec_bvalu(lambdaobs, bkpt, coeff)
         contdiv = spec/cont
         contdivivar = ivar*cont^2

         outstr = {lambda:lambdaobs,contdiv:contdiv,contdivivar:contdivivar,dlam:dlam,z:redshift,zmet:zmet(iz),age:age(ia),sn:sn(isn),element:element,element_abund:alpha(iap)}

         outname = dir+prefix+zmetname(iz)+agename(ia)+alphaname(iap)+snname(isn)+idname(ispec)+'.fits'
         mwrfits,outstr,outname,/create,/silent
      endfor
   endfor
end
