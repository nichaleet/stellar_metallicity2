pro mkfake_deimos_spec,notput=notput
   ;create fake deimos spectra of cl1604
   zmet = [-0.75,-0.55,-0.25,-0.05,0.15] 
   age = [1.,3.,5.]
   sn  = [5,10,15.,20.,30.] ;per angstrom 
   ;above were restricted to 2 decimal points, integer and integer
   ;accordingly (due to naming issue) but can be fixed easily if you want
   nspec = 5
   redshift = 0.9
   vdisp = 300. ;km/s in FWHM

   ;setting of DEIMOS spectra
   dwl = 0.16 ;dispersion per pixel in angstrom
   iwave = 3400.
   nwl = 12500
   nwlreal = 8500
   dlam = fltarr(nwl)+1.6 ;FWHM in Angstrom       
   lambdarest = findgen(nwl)*dwl+iwave
   lambdaobs = lambdarest*(1.+redshift)
   snperpix = sn*sqrt(dwl)
   ;setup for output names
   dir = '/scr2/nichal/workspace4/test_shortwl/mockdata/'
   prefix = 'mockdeimos_cl1604'
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
   idname =  '_'+strtrim(string(indgen(nspec)+1,format='(I)'),2)

   ;constant
   clight = 299792.458
   readcol,'sens_1200G_03OCT2017_0054.dat',wltput,tput,comment='#'
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
   for iz=0,n_elements(zmet)-1 do for ia=0,n_elements(age)-1 do begin
      spsstruct = sps_interp(zmet(iz), age(ia))
      spslambda = spsstruct.lambda
      spsspec = spsstruct.spec
      w = where(spslambda gt 3000 and spslambda lt 6000, c)
      if c lt 25 then message, 'Not enough pixels.'
      spslambda = spslambda[w]
      spsspec = spsspec[w]
      spsspec = spsspec*clight/spslambda^2      ;change fnu(Lsun/Hz) to flambda 
      spsmedian = median(spsspec)
      spsspec = spsspec/spsmedian      ;normalize to around 1
      ;smooth to vdisp and deimos resolution
      spsspec = smooth_gauss_wrapper(spslambda, spsspec, spslambda, vdisp/clight/2.35*spslambda)
      spsspec = smooth_gauss_wrapper(spslambda*(1.+redshift), spsspec, lambdaobs, dlam)
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

         outstr = {lambda:lambdaobs[0:nwlreal-1],contdiv:contdiv[0:nwlreal-1],contdivivar:contdivivar[0:nwlreal-1],dlam:dlam[0:nwlreal-1],lambdafull:lambdaobs,contdivfull:contdiv,contdivivarfull:contdivivar,dlamfull:dlam,z:redshift}
         outname = dir+prefix+zmetname(iz)+agename(ia)+snname(isn)+idname(ispec)+'.fits'
         mwrfits,outstr,outname,/create,/silent
      endfor
   endfor
end
