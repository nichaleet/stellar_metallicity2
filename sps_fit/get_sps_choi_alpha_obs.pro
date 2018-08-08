function get_sps_choi_alpha_obs, xin, a
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    common get_sps_alpha, element
    z = a[0]
    age = a[1]
    vdisp = a[2]
    redshift = a[3]
    if n_Elements(a) ge 5 then begin
       yesalpha=1
       if n_elements(a) eq 5 then begin
          alpha = a[4]
          nelements = n_elements(element)
          abund = rebin([alpha],nelements)
       endif
       if n_elements(a) eq 6 then begin
          abund = a[4:5]
       endif
    endif else yesalpha=0
 
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

    ;add response function to [a/Fe]
    if yesalpha then begin
       spsspec = add_response(lambda,spsspec,[z,age],element,abund,/silent)
    endif
    ;smooth to data wavelengths
    spsspec = smooth_gauss_wrapper(lambda*(redshift+1.), spsspec, datalam, sqrt(dlam^2+(vdisp/clight*datalam)^2))
    lambda  = datalam           ;datalam is science.lambda

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;fit continuum to synthetic spectra
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if normalize eq 1 or rest eq 1 then begin
       spline = 1
       poly = 0
       simpoly = 0

       spsspec_original = spsspec
       won = where(contmask eq 1, con)
       if con lt 25 then message, 'Not enough pixels.'
       invvar=dataivar[won]

       case 1 of
          spline: begin
             bkpt = slatec_splinefit(lambda[won]/(1.+redshift), spsspec[won], coeff, bkspace=165,invvar=invvar,upper=3, lower=3, /silent,/everyn,mask=mask) ;SDSS
             
             if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
             cont = slatec_bvalu(lambda/(1.+redshift), bkpt, coeff)
          end
          poly: begin
             degree = 7
             norm = median(spsspec[won])
             a = [norm, replicate(0.0, degree-1)]
             p = lmfit(lambda[won]/(1.+redshift),spsspec[won], a, measure_errors=(invvar)^(-0.5), /double, function_name='legendre_poly')
             cont = legendre_poly(lambda/(1.+redshift), a, /noderiv)
          end
          simpoly:begin
             npoly = 6
             degree = npoly
             p=poly_fit(lambda[won]/(1.+redshift),spsspec[won],degree)
             cont = poly(lambda/(1.+redshift),p)
          end
       endcase
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
