function get_sps_alpha_obs, xin, a
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
    common get_sps_alpha, element
    z = a[0]
    age = a[1]
    vdisp = a[2]
    redshift = a[3]
    nelements= n_elements(element)
    if n_Elements(a) ge 5 then begin
       if n_elements(a) eq 5 then begin
          alpha = a[4]
          nelements = n_elements(element)
          abund = rebin([alpha],nelements)
       endif else begin
          abund = a[4:n_elements(a)-1]
       endelse
    endif

    ;spsspec = dblarr(n_elements(xin)) - 99999.
    ;if z lt min(spsz) or z gt max(spsz) then return, spsspec
    ;if age lt min(spsage) or age gt max(spsage) then return, spsspec

    spsstruct = sps_interp(z, age)
    lambda = spsstruct.lambda
    spsspec = spsstruct.spec

    w = where(lambda gt 3000 and lambda lt 8000, c)
    if c lt 25 then message, 'Not enough pixels.'
    lambda = lambda[w]
    spsspec = spsspec[w]
    clight = 299792.458

    spsspec = spsspec*clight/lambda^2    ;change fnu(Lsun/Hz) to flambda
    spsspec = spsspec/median(spsspec)    ;normalize to around 1
    ;add response function
    spsspec = add_response(lambda,spsspec,[z,age],element,abund,/silent)
   ;smooth to data wavelengths
;    spsspec = smooth_gauss_wrapper(lambda, spsspec, lambda, vdisp/clight*lambda)
;    spsspec = smooth_gauss_wrapper(lambda*(1.+redshift), spsspec, datalam, dlam)
    spsspec = smooth_gauss_wrapper(lambda*(redshift+1.), spsspec, datalam, sqrt(dlam^2+(vdisp/clight*datalam)^2))
;function smooth_gauss_wrapper, lambda1, spec1, lambda2, dlam_in, ivar1=ivar1, ivar2=ivar2
    lambda  = datalam ;datalam is science.lambda
  
    if normalize eq 1 or rest eq 1 then begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;fit continuum to synthetic spectra
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       spsspec_original = spsspec
       n = n_elements(lambda)
       nhalf = round(double(n)/2.)
       if n gt 8190 and n lt 8195 then begin
          wwhole = lindgen(nhalf)
          ccd2 = 2
       endif else begin
          wwhole = lindgen(n)
          ccd2 = 1
       endelse
       for ccd=1,ccd2 do begin
          won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
          woff += wwhole[0]   
          if con lt 25 then message, 'Not enough pixels.'
          invvar = dataivar[won]/(median(spsspec[won]))^2
          ;;fit continuum to rest wavelength
          bkpt = slatec_splinefit(lambda[won]/(1.+redshift), spsspec[won], coeff, invvar=invvar, bkspace=150, upper=3, lower=3, /silent) ;DEIMOS
          if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
          cont = slatec_bvalu((lambda[wwhole])/(1.+redshift), bkpt, coeff)
          
          if ccd eq 1 then contb = cont
          if ccd eq 2 then contr = cont
          wwhole += nhalf
       endfor
       if ccd2 eq 2 then cont = [contb, contr] else cont = contb
       
       spsspec /= cont
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if rest eq 0 then begin
       ;;only return the unmasked range (same elements as xin)
       spsspec = spsspec[wonfit]
       if total(xin-lambda(wonfit)) ne 0 then stop
       return, spsspec
    endif

    if rest eq 1 then begin ; This is to replace the get_sps_rest.pro
       if total(xin-lambda) ne 0 then stop
       return,[[spsspec],[spsspec_original],[cont]]
    endif
end
