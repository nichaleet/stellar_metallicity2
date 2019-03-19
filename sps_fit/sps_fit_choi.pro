pro poly6,x,a,f,pder
   f=a[0]+a[1]*x+a[2]*x^2+a[3]*x^3+a[4]*x^4+a[5]*x^5
end
; =================
pro sps_iterproc, funcname, p, iter, fnorm, functargs=functargs, parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if iter gt 1 then begin
       print, contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       if contiter mod 10 eq 0 then begin
          printf,long(copynum),contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       endif
    endif
end

pro sps_iterproc_alpha, funcname, p, iter, fnorm, functargs=functargs, $
                         parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if iter gt 1 then begin
       print, contiter, p[0], p[1], p[2],p[3],p[4], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D5.2,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       if contiter mod 10 eq 0 then begin
       printf,long(copynum),contiter, p[0], p[1], p[2],p[3],p[4], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D5.2,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       endif
    endif
end

pro sps_iterproc_omg, funcname, p, iter, fnorm, functargs=functargs, $
                         parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if iter gt 1 then begin
       print, contiter, p[0], p[1], p[2],p[3],p[4],p[5], fnorm/dof, dof,format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D5.2,2x,D6.3,1x,D6.3,1X,D10.5,1X,I4)'
       if contiter mod 10 eq 0 then begin
       printf,long(copynum),contiter, p[0], p[1], p[2],p[3],p[4],p[5], fnorm/dof, dof,format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D5.2,2x,D6.3,1x,D6.3,1X,D10.5,1X,I4)'
       endif
    endif
end

pro sps_fit::fit, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, npoly,contmask, normalize, rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if science.npoly eq 0 then npoly = 4 else npoly = science.npoly
    
    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting ...'

    znow = science.zspec
    reallambda = science.lambda
    nlambda = n_elements(reallambda)
    
    specresfile = self.directory+'specres_poly.sav'
    if file_test(specresfile) then begin
       restore, specresfile
       dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
    endif else begin
       if science.mask eq 'NGC6528' then dlam_all = reallambda*0.+3.1/2.35 else $
          dlam_all = reallambda/1000.
    endelse 
    
    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
    pi[0].limits = [-0.5,0.02]
    pi[0].limits = [-0.1,0.]
    pi[1].limits = [min(spsage),alog10(galage(znow,1000))<max(spsage)]
;    pi[1].limits = [min(spsage),4.0]
    ;pi[1].limits = minmax(spsage)
    pi[2].limits = [0., 500.0]
    pi[3].limits = [-0.1,0.1]+znow

   ;;make the initial guesses unfix but within limits except redshift
    if science.feh eq -999 then pi.value = randomu(seed,4)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*] else pi.value=[science.feh,science.age,science.vdisp,znow]
    pi.value=[science.fehchoi,science.agechoi,200,znow]
    pi[3].value = znow
    firstguess = pi.value    
    pi[0].limits = minmax(spsz)
    pi.step = double([0.1, 0.5, 25.0,0.0002])
    pi.parname = ['    Z', '  age', 'vdisp','redshift']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
    print, 'prior range:',pi.limits

    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d]
        perror = [-999d, -999d, -999d]
        science.spsspec = -999d
        goto, done
    endif    

    xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0

    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 150
    maxnloop = 100
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
    print, '--- ------- ----- ------- ---------  -------- ----'

    openw,copynum,'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_choi'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF   ZDIFF  AGEDIFF'
    printf,long(copynum), '--- ------- ----- ------- ---------  -------- ---- -----  ------'

    ;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.]
    besterror = [99.,99.,99.,99.]
    chisqarr = fltarr(maxnloop)
    valuearr = fltarr(4,maxnloop)
    errorarr = fltarr(4,maxnloop)

    ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop lt 30 do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam  = science.lambda
        wonfit   = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_choi_alpha_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        if status eq 0 then begin
           perror=pi.value*0.-99.
           stop
        endif
        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        vdispdiff = (pi[2].value-pars[2])/pi[2].value
        redshfdiff = (pi[3].value-pars[3])/pi[3].value
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)
        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_choi_alpha_obs(reallambda, pars) 
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

        ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fit 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont  
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.spec/median(science.spec)
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')

        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=165, upper=3, lower=3, /silent)
        
        if bkpt[0] eq -1 then begin
           pi.value = [-999d, -999d, -999d,-999d]
           perror = [-999d, -999d, -999d, -999d]
           science.spsspec = -999d
           science.spscont = -999d 
           break
        endif

        cont = slatec_bvalu(restlambda, bkpt, coeff)
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
        
        plot,xmp/(1.+pi[3].value),ympold
        oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
        oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        curchisq = bestnorm/dof

        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        chisqarr[nloop] = curchisq
        valuearr[*,nloop] = pars
        errorarr[*,nloop] = perror
        nloop +=1
     endwhile
    close,copynum
    ;check if the last chisq is the best chisq
    if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
       print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
       print,'The values used are:'
       print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3)'
       science.goodfit = 1.
       spsbestfitarr = bestspsbestfitarr
       spsbestfit=spsbestfitarr[*,1]
       pi.value = bestvalue
       perror   = besterror
       science.spscont = bestcont
    endif

    science.spsspec = spsbestfitarr[*,0] 
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '
    
    done:
    science.nloop = nloop
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    science.alphafe = 0.
    science.alphafeerr = -999.
    science.alphafeupper = -999.
    science.alphafelower = -999.

    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
 end

pro sps_fit::fitalphaband, science,noredraw=noredraw
    common sps_spec, sps, spsz, spsage
    common mask_in, mask_in, copynum
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting Mg...'
    znow = science.zfit
    reallambda = science.lambda
    nlambda = n_elements(reallambda)
    restlambda = reallambda / (1d + science.zfit)

    specresfile = self.directory+'specres_poly.sav'
    if file_test(specresfile) then begin
       restore, specresfile
       dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
    endif else dlam_all = reallambda/1000.
    dlam = dlam_all

    maskold = science.fitmask
    self->mask,science,/includemg 
    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and restlambda gt 5100. and restlambda lt 5200., con)
    if con lt 10 then begin
       science.alphafe = -999d
       science.alphafeerr = -999d
       goto, done
    endif
    ;get the best alpha=0 values
    spsspec = science.spsspecfull
    datalam = science.lambda
    z = science.feh
    age = science.age
    element = ['Mg']
    nelements = n_elements(element)
    ;get the response function    
    afearr = findgen(51)*0.016-0.4 ;-0.4 to 0.4
    nafearr = n_Elements(afearr)
    ;get chisq grid for each alpha
    afe_chisq = fltarr(nafearr)
    for i=0,nafearr-1 do begin
       alpha = afearr[i]
       spec = add_response(restlambda,spsspec,[z,age],element,rebin([alpha],nelements),/silent)
       afe_chisq[i] = total((spec[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won)-1)
    endfor
    ;return back to the previous mask
    science.fitmask = maskold
    ;get the probability function from chisq
    likelihood = -0.5*afe_chisq
    likelihood = likelihood-max(likelihood)
    probarr = exp(likelihood)
    area = int_tabulated(afearr,probarr)
    probarr = probarr/area
    bestafe = confidence_interval(afearr,probarr)

    science.alphafe = bestafe[0]
    science.alphafelower = bestafe[4]
    science.alphafeupper = bestafe[5]
    science.alphafeerr = (bestafe[5]-bestafe[4])/2.

    ;Printing
    print,'[Mg/Fe] results: Best, first qt, median, 3rd qt,lower limit, upper limit'
    print, strjoin(string(bestafe,format='(F6.3)'),' ')

    ;get the best model
    spec = add_response(restlambda,spsspec,[z,age],element,rebin([science.alphafe],nelements),/silent)
    science.spsspecfullmg = science.spsspecfull
    science.spsspecfullmg(won) = spec(won)

    ;calculate chisq
    science.chisqmg = total((science.spsspecfullmg[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    done:
    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
     endif
end

pro sps_fit::fitalpha, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting 5 parameters ...'
    widget_control, widget_info(self.base, find_by_uname='keepoldfit'), get_value=keepoldfit
    oldchisq = science.chisq
    if n_elements(element) eq 0 then element= ['Mg','O','Si','Ca','Ti']
    znow = science.zspec
    if znow le 0. then znow = science.z
    if znow le 0. then stop
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    ;check if the blue chip was failed
    neg = where(reallambda lt 0., cneg)
    if cneg gt 0 then reallambda(neg) = reallambda(neg)+10000.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    dlam_all = science.dlam*0.
    baddlam = where(~finite(dlam_all),cbaddlam,complement=gooddlam)
    if cbaddlam gt 0 then begin
       dlam_all(baddlam) = interpol(dlam_all(gooddlam),reallambda(gooddlam),reallambda(baddlam))
    endif

    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)

    pi[0].limits = [-0.6,0.19] ;this is just for initial parameters
    pi[1].limits = [min(spsage),6];(galage(znow,1000)/1.e9)<max(spsage)]
   ; pi[1].limits=[2,3.]
    pi[3].limits = [-0.05,0.05]+znow
    if pi[3].limits[0] lt 0. then pi[3].limits[0]=0
    pi[4].limits =[0,0.4]
    ;set the prior to velocity dispersion according to Faber-Jackson relation (Dutton2011)
;    if science.logmstar gt 5. then begin
;       logvdisp = 2.23+0.37*(science.logmstar-10.9)-0.19*alog10(0.5+0.5*(10.^science.logmstar/10.^10.9))
;       pi[2].limits = [10.^(logvdisp-0.4),10.^(logvdisp+0.4)]
;       print, 'velocity dispersion prior = ',pi[2].limits,' km/s'      
;    endif else pi[2].limits = [40.,400.]; [30.,200.]
    pi[2].limits = [200.,500.]
   ; pi[2].limits = [348.,352.]
    if mask_in eq 'stacked' then pi[2].limits = [50.,400.]
   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,5)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    firstguess = pi.value
    pi[0].limits = minmax(spsz) ;fix the limit of [Fe/H] back
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    pi[4].limits =[-0.4,0.4]
    pi.step = double([0.1, 0.5, 25.0,0.002,0.1])
    pi.parname = ['    Z', '  age', 'vdisp','redshift','alpha']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)']
    print, 'prior range:',pi.limits
    print, 'alpha elements are ',element
    if copynum eq '10' then bkspace = 300 else bkspace =165
    print, 'bkspace', bkspace
    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d,-999d,-999d]
        perror = [-999d, -999d, -999d, -999d,-999d]
        science.spsspec = -999d
        goto, done
    endif
   xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 100
    maxnloop = 100
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift  alpha   chi^2  DOF   ZDIFF  AGEDIFF'
    print, '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'

    openw,long(copynum),'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_choialpha'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), 'prior range:',pi.limits
    printf,long(copynum), 'alpha elements are ',element
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift  alpha   chi^2  DOF   ZDIFF  AGEDIFF'
    printf,long(copynum), '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'
;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.,99.]
    besterror = [99.,99.,99.,99.,99.]
   ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop le maxnloop do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_choi_alpha_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc_alpha')

        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_choi_alpha_obs(reallambda, pars)
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

       ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
;        wset, index
;        spectoplot = science.contdiv/science.spscont
;        !p.multi=[0,1,2]
;        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
;        oplot,restlambda,spectoplot
;        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')

;        npoly = 6
;        degree = npoly
;        p = poly_fit(restlambda[won],science.contdiv[won]/spsbestfit[won],degree)
;;        yfit= curvefit(restlambda[won],science.contdiv[won]/spsbestfit[won],science.contdivivar[won]*(spsbestfit[won])^2,p,function_name='poly6')
;        cont = poly(restlambda,p)

        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=bkspace, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d, -999d,-999d]
            perror = [-999d, -999d, -999d, -999d,-999d]
            science.spsspec = -999d
            science.spscont = -999d
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
       
;        plot,restlambda[won],science.contdiv[won]/spsbestfit[won]
;        oplot, restlambda,cont,color=fsc_color('purple')
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
 ;       plot,xmp/(1.+pi[3].value),ympold
 ;       oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
 ;       oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        nloop +=1
 ;       if nloop eq 8 and keepoldfit eq 0 and curchisq gt 1.5 then maxnloop = 10
 ;       if nloop eq 18 and keepoldfit eq 0 and curchisq-oldchisq gt 1. then maxnloop = 20
 ;       if nloop eq 48 and keepoldfit eq 0 and curchisq-oldchisq gt 0.2 then maxnloop = 50
 ;       if nloop eq 98 and keepoldfit eq 0 and curchisq-oldchisq gt 0.1 then maxnloop = 100
     endwhile
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'
    ;check if the last chisq is the best chisq
    if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
       print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
       print,'The values used are:'
       print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestvalue[4],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1x,D6.3,1X,D8.3)'
       science.goodfit = 1.
       spsbestfitarr = bestspsbestfitarr
       spsbestfit = spsbestfitarr[*,1]
       pi.value = bestvalue
       perror   = besterror
       science.spscont = bestcont
    endif

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0]
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '

    done:
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.zspec = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    science.alphafe = pi[4].value
    science.alphafeerr = perror[4]
    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    print, science.chisq
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    ;self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fitomg, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting 5 parameters ...'
    widget_control, widget_info(self.base, find_by_uname='keepoldfit'), get_value=keepoldfit
    oldchisq = science.chisq
    element= ['Mg','Fe']
    znow = science.zspec
    if znow le 0. then znow = science.z
    if znow le 0. then stop
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    ;check if the blue chip was failed
    neg = where(reallambda lt 0., cneg)
    if cneg gt 0 then reallambda(neg) = reallambda(neg)+10000.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    dlam_all = science.dlam*0.
    baddlam = where(~finite(dlam_all),cbaddlam,complement=gooddlam)
    if cbaddlam gt 0 then begin
       dlam_all(baddlam) = interpol(dlam_all(gooddlam),reallambda(gooddlam),reallambda(baddlam))
    endif

    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 6)

    pi[0].limits = [-0.6,0.19] ;this is just for initial parameters
    pi[1].limits = [min(spsage),6];(galage(znow,1000)/1.e9)<max(spsage)]
   ; pi[1].limits=[2,3.]
    pi[3].limits = [-0.05,0.05]+znow
    if pi[3].limits[0] lt 0. then pi[3].limits[0]=0
    pi[4].limits =[0,0.4]
    pi[5].limits = [0,0.8]
    ;set the prior to velocity dispersion according to Faber-Jackson relation (Dutton2011)
;    if science.logmstar gt 5. then begin
;       logvdisp = 2.23+0.37*(science.logmstar-10.9)-0.19*alog10(0.5+0.5*(10.^science.logmstar/10.^10.9))
;       pi[2].limits = [10.^(logvdisp-0.4),10.^(logvdisp+0.4)]
;       print, 'velocity dispersion prior = ',pi[2].limits,' km/s'      
;    endif else pi[2].limits = [40.,400.]; [30.,200.]
    pi[2].limits = [200.,500.]
   ; pi[2].limits = [348.,352.]
    if mask_in eq 'stacked' then pi[2].limits = [50.,400.]
   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,6)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    firstguess = pi.value
    pi[0].limits = minmax(spsz) ;fix the limit of [Fe/H] back
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    pi[4].limits =[-0.4,0.4]
    pi[5].limits =[-0.4,0.8]
    pi.step = double([0.1, 0.5, 25.0,0.002,0.1,0.1])
    pi.parname = ['    Z', '  age', 'vdisp','redshift','Mg','O']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)','(D6.3)']
    print, 'prior range:',pi.limits
    print, 'alpha elements are ',element
    bkspace =165
    print, 'bkspace', bkspace
    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d,-999d,-999d,-999d]
        perror = [-999d, -999d, -999d, -999d,-999d,-999d]
        science.spsspec = -999d
        goto, done
    endif
   xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 150
    maxnloop = 100
    print, 'maxnloop ', maxnloop
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift  Mg  N  chi^2  DOF'
    print, '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'

    openw,long(copynum),'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_choiomg'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), 'prior range:',pi.limits
    printf,long(copynum), 'alpha elements are ',element
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift  Mg  O   chi^2  DOF'
    printf,long(copynum), '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'
;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.,99.,99]
    besterror = [99.,99.,99.,99.,99.,99]
   ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop le maxnloop do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_choi_alpha_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc_omg')

        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_choi_alpha_obs(reallambda, pars)
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

       ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.contdiv/science.spscont
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')

;        npoly = 6
;        degree = npoly
;        p = poly_fit(restlambda[won],science.contdiv[won]/spsbestfit[won],degree)
;;        yfit= curvefit(restlambda[won],science.contdiv[won]/spsbestfit[won],science.contdivivar[won]*(spsbestfit[won])^2,p,function_name='poly6')
;        cont = poly(restlambda,p)

        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=bkspace, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d, -999d,-999d,-999d]
            perror = [-999d, -999d, -999d, -999d,-999d,-999d]
            science.spsspec = -999d
            science.spscont = -999d
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
       
;        plot,restlambda[won],science.contdiv[won]/spsbestfit[won]
;        oplot, restlambda,cont,color=fsc_color('purple')
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
 ;       plot,xmp/(1.+pi[3].value),ympold
 ;       oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
 ;       oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        nloop +=1
        if nloop eq 8 and keepoldfit eq 0 and curchisq gt 5 then maxnloop = 10
 ;       if nloop eq 18 and keepoldfit eq 0 and curchisq-oldchisq gt 1. then maxnloop = 20
 ;       if nloop eq 48 and keepoldfit eq 0 and curchisq-oldchisq gt 0.2 then maxnloop = 50
 ;       if nloop eq 98 and keepoldfit eq 0 and curchisq-oldchisq gt 0.1 then maxnloop = 100
     endwhile
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'
    ;check if the last chisq is the best chisq
    if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
       print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
       print,'The values used are:'
       print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestvalue[4],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1x,D6.3,1X,D8.3)'
       science.goodfit = 1.
       spsbestfitarr = bestspsbestfitarr
       spsbestfit = spsbestfitarr[*,1]
       pi.value = bestvalue
       perror   = besterror
       science.spscont = bestcont
    endif

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0]
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '

    done:
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.zspec = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    science.alphafe = pi[4].value
    science.alphafeerr = perror[4]
    science.pfit[0:5] = pi.value
    science.pfiterr[0:5] = perror
    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    print, science.chisq
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    ;self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end


pro sps_fit::fitmgfeo, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting 5 parameters ...'
    widget_control, widget_info(self.base, find_by_uname='keepoldfit'), get_value=keepoldfit
    oldchisq = science.chisq
    element= ['Mg','Fe','N']
    znow = science.zspec
    if znow le 0. then znow = science.z
    if znow le 0. then stop
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    ;check if the blue chip was failed
    neg = where(reallambda lt 0., cneg)
    if cneg gt 0 then reallambda(neg) = reallambda(neg)+10000.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    dlam_all = science.dlam*0.
    baddlam = where(~finite(dlam_all),cbaddlam,complement=gooddlam)
    if cbaddlam gt 0 then begin
       dlam_all(baddlam) = interpol(dlam_all(gooddlam),reallambda(gooddlam),reallambda(baddlam))
    endif

    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 7)

    pi[0].limits = [-0.6,0.19] ;this is just for initial parameters
    pi[1].limits = [min(spsage),6];(galage(znow,1000)/1.e9)<max(spsage)]
   ; pi[1].limits=[2,3.]
    pi[3].limits = [-0.05,0.05]+znow
    if pi[3].limits[0] lt 0. then pi[3].limits[0]=0
    pi[4].limits =[0,0.4]
    pi[5].limits = [-0.3,0.3]
    pi[6].limits = [0,0.8]
    pi[2].limits = [200.,500.]
   ; pi[2].limits = [348.,352.]
    if mask_in eq 'stacked' then pi[2].limits = [50.,400.]
   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,6)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    firstguess = pi.value
    pi[0].limits = minmax(spsz) ;fix the limit of [Fe/H] back
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    pi[4].limits =[-0.4,0.4]
    pi[6].limits =[-0.4,0.8]
    pi.step = double([0.1, 0.5, 25.0,0.002,0.1,0.1,0.1])
    pi.parname = ['    Z', '  age', 'vdisp','redshift',element]
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)','(D6.3)','(D6.3)']
    print, 'prior range:',pi.limits
    print, 'alpha elements are ',element
    bkspace =165
    print, 'bkspace', bkspace
    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d,-999d,-999d,-999d,-999d]
        perror = [-999d, -999d, -999d, -999d,-999d,-999d,-999d]
        science.spsspec = -999d
        goto, done
    endif
   xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 150
    maxnloop = 100
    print, 'maxnloop ', maxnloop
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift  Mg  Fe  chi^2  DOF'
    print, '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'

    openw,long(copynum),'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_choimgfeo'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), 'prior range:',pi.limits
    printf,long(copynum), 'alpha elements are ',element
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift  Mg  Fe   chi^2  DOF'
    printf,long(copynum), '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'
;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.,99.,99.,99]
    besterror = [99.,99.,99.,99.,99.,99.,99]
   ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop le maxnloop do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_choi_alpha_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc_omg')

        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_choi_alpha_obs(reallambda, pars)
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

       ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.contdiv/science.spscont
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')

;        npoly = 6
;        degree = npoly
;        p = poly_fit(restlambda[won],science.contdiv[won]/spsbestfit[won],degree)
;;        yfit= curvefit(restlambda[won],science.contdiv[won]/spsbestfit[won],science.contdivivar[won]*(spsbestfit[won])^2,p,function_name='poly6')
;        cont = poly(restlambda,p)

        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=bkspace, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d, -999d,-999d,-999d,-999d]
            perror = [-999d, -999d, -999d, -999d,-999d,-999d,-999d]
            science.spsspec = -999d
            science.spscont = -999d
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
       
;        plot,restlambda[won],science.contdiv[won]/spsbestfit[won]
;        oplot, restlambda,cont,color=fsc_color('purple')
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
 ;       plot,xmp/(1.+pi[3].value),ympold
 ;       oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
 ;       oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        nloop +=1
        if nloop eq 8 and keepoldfit eq 0 and curchisq gt 5 then maxnloop = 10
 ;       if nloop eq 18 and keepoldfit eq 0 and curchisq-oldchisq gt 1. then maxnloop = 20
 ;       if nloop eq 48 and keepoldfit eq 0 and curchisq-oldchisq gt 0.2 then maxnloop = 50
 ;       if nloop eq 98 and keepoldfit eq 0 and curchisq-oldchisq gt 0.1 then maxnloop = 100
     endwhile
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'
    ;check if the last chisq is the best chisq
    if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
       print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
       print,'The values used are:'
       print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestvalue[4],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1x,D6.3,1X,D8.3)'
       science.goodfit = 1.
       spsbestfitarr = bestspsbestfitarr
       spsbestfit = spsbestfitarr[*,1]
       pi.value = bestvalue
       perror   = besterror
       science.spscont = bestcont
    endif

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0]
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '

    done:
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.zspec = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    science.alphafe = pi[4].value
    science.alphafeerr = perror[4]
    science.pfit[0:6] = pi.value
    science.pfiterr[0:6] = perror
    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    print, science.chisq
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    ;self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fit_glorious, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting 13 parameters ...'
    widget_control, widget_info(self.base, find_by_uname='keepoldfit'), get_value=keepoldfit
    oldchisq = science.chisq
    element= ['Mg','N','Fe','O','C','N','Si','Ca','Ti'] ;9 elements
    znow = science.zspec
    if znow le 0. then znow = science.z
    if znow le 0. then stop
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    ;check if the blue chip was failed
    neg = where(reallambda lt 0., cneg)
    if cneg gt 0 then reallambda(neg) = reallambda(neg)+10000.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    dlam_all = science.dlam*0.
    baddlam = where(~finite(dlam_all),cbaddlam,complement=gooddlam)
    if cbaddlam gt 0 then begin
       dlam_all(baddlam) = interpol(dlam_all(gooddlam),reallambda(gooddlam),reallambda(baddlam))
    endif

    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 13)

    pi[0].limits = [-0.6,0.19] ;this is just for initial parameters
    pi[1].limits = [min(spsage),6];(galage(znow,1000)/1.e9)<max(spsage)]
   ; pi[1].limits=[2,3.]
    pi[3].limits = [-0.05,0.05]+znow
    if pi[3].limits[0] lt 0. then pi[3].limits[0]=0
    pi[4].limits =[0,0.4]
    pi[5].limits = [0,0.8]
    pi[6:12].limits = [-0.3,0.4]
    ;set the prior to velocity dispersion according to Faber-Jackson relation (Dutton2011)
;    if science.logmstar gt 5. then begin
;       logvdisp = 2.23+0.37*(science.logmstar-10.9)-0.19*alog10(0.5+0.5*(10.^science.logmstar/10.^10.9))
;       pi[2].limits = [10.^(logvdisp-0.4),10.^(logvdisp+0.4)]
;       print, 'velocity dispersion prior = ',pi[2].limits,' km/s'      
;    endif else pi[2].limits = [40.,400.]; [30.,200.]
    pi[2].limits = [200.,500.]
   ; pi[2].limits = [348.,352.]
    if mask_in eq 'stacked' then pi[2].limits = [50.,400.]
   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,13)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    firstguess = pi.value
    pi[0].limits = minmax(spsz) ;fix the limit of [Fe/H] back
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    pi[4].limits =[-0.4,0.5]
    pi.step = double([0.1, 0.5, 25.0,0.002,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1])
    pi.parname = ['    Z', '  age', 'vdisp','redshift','Mg','N','Fe','O','C','N','Si','Ca','Ti']
    pi.mpformat = ['(D6.3)','(D5.2)','(D6.1)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)']
    print, 'prior range:',pi.limits
    print, 'alpha elements are ',element
    bkspace =165
    print, 'bkspace', bkspace
    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
    if con lt 10 then begin
        pi.value = -999d
        perror = fltarr(13)-999. 
        science.spsspec = -999d
        goto, done
    endif
    xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 150
    maxnloop = 1000
    print, 'maxnloop ', maxnloop
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift  Mg  N  chi^2  DOF'
    print, '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'

    openw,long(copynum),'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_choiglorious'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), 'prior range:',pi.limits
    printf,long(copynum), 'alpha elements are ',element
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift  Mg  O   chi^2  DOF'
    printf,long(copynum), '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'
;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = fltarr(13)+99.
    besterror = fltarr(13)+99.
   ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop le maxnloop do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_choi_alpha_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc_omg')

        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_choi_alpha_obs(reallambda, pars)
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

       ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.contdiv/science.spscont
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')

;        npoly = 6
;        degree = npoly
;        p = poly_fit(restlambda[won],science.contdiv[won]/spsbestfit[won],degree)
;;        yfit= curvefit(restlambda[won],science.contdiv[won]/spsbestfit[won],science.contdivivar[won]*(spsbestfit[won])^2,p,function_name='poly6')
;        cont = poly(restlambda,p)

        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=bkspace, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = fltarr(13)-999.
            perror = fltarr(13)-999.
            science.spsspec = -999d
            science.spscont = -999d
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
       
;        plot,restlambda[won],science.contdiv[won]/spsbestfit[won]
;        oplot, restlambda,cont,color=fsc_color('purple')
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
 ;       plot,xmp/(1.+pi[3].value),ympold
 ;       oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
 ;       oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        nloop +=1
 ;       if nloop eq 8 and keepoldfit eq 0 and curchisq gt 5 then maxnloop = 10
 ;       if nloop eq 18 and keepoldfit eq 0 and curchisq-oldchisq gt 1. then maxnloop = 20
 ;       if nloop eq 48 and keepoldfit eq 0 and curchisq-oldchisq gt 0.2 then maxnloop = 50
 ;       if nloop eq 98 and keepoldfit eq 0 and curchisq-oldchisq gt 0.1 then maxnloop = 100
     endwhile
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'
    ;check if the last chisq is the best chisq
    if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
       print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
       print,'The values used are:'
       print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestvalue[4],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1x,D6.3,1X,D8.3)'
       science.goodfit = 1.
       spsbestfitarr = bestspsbestfitarr
       spsbestfit = spsbestfitarr[*,1]
       pi.value = bestvalue
       perror   = besterror
       science.spscont = bestcont
    endif

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0]
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '

    done:
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.zspec = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    science.alphafe = pi[4].value
    science.alphafeerr = perror[4]

    science.pfit = pi.value
    science.pfiterr = perror
    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    print, science.chisq
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    ;self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::cal_uncertainties, science
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Calculatign uncertainties ...'
   grid_file = *self.degen_file   ;nage
   grid_feh = *self.degen_fehgrid
   grid_age = *self.degen_agegrid
   grid_alpha = *self.degen_alphagrid ;nfeh x nage x nalpha
   loc_alpha = where(grid_alpha[0,0,*] eq 0.,cloc_alpha)
   loc_ext = loc_alpha+1
   if cloc_alpha ne 1 then stop,'Cannot find where alpha/Fe = 0'
   grid_Feh = grid_feh[*,*,loc_alpha]
   grid_age = grid_age[*,*,loc_alpha]

   ;find where is the spec closest to the input age and feh and alpha
   min_dist = min(sqrt((grid_feh-science.feh)^2+(grid_age-science.age)^2),iref)
   if min_dist gt 0.05 then stop,'halted because matching might be off grid in make_chisqarr_new.pro'
   loc = array_indices(grid_feh,iref)
   ;make noisy ref spectra
   specarr = mrdfits(grid_file(loc[1]),loc_ext[0],/silent)
   refspecstr = specarr(loc[0])
   if refspecstr.feh ne grid_feh(iref) or refspecstr.age ne grid_age(iref) then stop,'grid does not match with file'
   npix = n_elements(refspecstr.lambda)
   refspec_err = abs(refspecstr.spec/science.snfit)
   refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

   ;make chisqarr
   gridsize = size(grid_feh,/dimensions)
   chisqarr = fltarr(gridsize)

   for ia = 0,n_elements(grid_file)-1 do begin
      specarr = mrdfits(grid_file(ia),loc_ext[0],/silent)
      if gridsize[0] ne n_elements(specarr) then stop,'grid does not match size'
      for i=0,n_elements(specarr)-1 do begin
         loc = where(grid_feh eq specarr[i].feh and grid_age eq specarr[i].age,cloc)
         if cloc ne 1 then stop,'oops'
         ind = array_indices(grid_feh,loc)
         chisqarr[ind[0],ind[1]] = total((refspec-specarr[i].spec)^2/refspec_err)
      endfor
   endfor

   feharr = grid_feh[*,0]
   agearr = reform(grid_age[0,*])

   deltafeh = grid_feh[1,0]-grid_feh[0,0]
   deltaage = grid_age[0,1]-grid_age[0,0]
   ;read the chisq grid and calculate probability
   Lgrid = -0.5*(chisqarr-min(chisqarr))
   probgrid = exp(double(Lgrid))
   volume = 0
   arr_dimen = size(grid_feh,/dimension)
   for ii=0,arr_dimen(0)-1 do begin
      for jj=0,arr_dimen(1)-1 do begin
         volume = volume+deltafeh*deltaage*probgrid[ii,jj]
      endfor
   endfor
   probgrid = probgrid/volume

   probfeh = fltarr(arr_dimen(0))
   probage = fltarr(arr_dimen(1))
   for ii=0,arr_dimen(0)-1 do probfeh(ii) = int_tabulated(agearr,probgrid[ii,*])
   for jj=0,arr_dimen(1)-1 do probage(jj) = int_tabulated(feharr,probgrid[*,jj])
   probfeh = probfeh/int_tabulated(feharr,probfeh)
   probage = probage/int_tabulated(agearr,probage)

   cumprobfeh = fltarr(arr_dimen(0))
   for jj=1,arr_dimen(0)-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
   cumprobage = fltarr(arr_dimen(1))
   for jj=1,arr_dimen(1)-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])

   midagevalue = interpol(agearr,cumprobage,0.5)
   midfehvalue = interpol(feharr,cumprobfeh,0.5)
   science.ageupper = science.age+(interpol(agearr,cumprobage,0.84)-midagevalue)
   science.agelower = science.age+(interpol(agearr,cumprobage,0.16)-midagevalue)
   science.fehupper = science.feh+(interpol(feharr,cumprobfeh,0.84)-midfehvalue)
   science.fehlower = science.feh+(interpol(feharr,cumprobfeh,0.16)-midfehvalue)
   science.alphafeupper = -999.
   science.alphafelower = -999.
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready ...'
end
pro sps_fit::cal_uncertainties_alpha, science
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Calculatign uncertainties ...'
   grid_file = *self.degen_file   ;nage
   grid_feh = *self.degen_fehgrid ;nfeh x nage x nalpha
   grid_age = *self.degen_agegrid ;nfeh x nage x nalpha
   grid_alpha = *self.degen_alphagrid ;nfeh x nage x nalpha
   ;find where is the spec closest to the input age and feh and alpha
   min_dist = min(sqrt((grid_feh-science.feh)^2+(grid_age-science.age)^2+(grid_alpha-science.alphafe)^2),iref)
   if min_dist gt 0.1 then print,strtrim(string(min_dist))+' matching might be off grid'
   loc = array_indices(grid_feh,iref)

   ;make noisy ref spectra
   specarr = mrdfits(grid_file(loc[1]),loc[2]+1,/silent)
   refspecstr = specarr(loc[0])
   if refspecstr.feh ne grid_feh(iref) or refspecstr.age ne grid_age(iref) or refspecstr.alpha ne grid_alpha(iref) $
        then stop,'grid does not match with file'
   npix = n_elements(refspecstr.lambda)
   refspec_err = abs(refspecstr.spec/science.snfit)
   refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

   ;make chisqarr
   gridsize = size(grid_feh,/dimensions)
   chisqarr = fltarr(gridsize)

   for ia = 0,n_elements(grid_file)-1 do for iap = 1,gridsize[2] do begin
      specarr = mrdfits(grid_file(ia),iap,/silent)
      if gridsize[0] ne n_elements(specarr) then stop,'grid does not match size'
      for i=0,n_elements(specarr)-1 do begin
         loc = where(grid_feh eq specarr[i].feh and grid_age eq specarr[i].age and $
               grid_alpha eq specarr[i].alpha,cloc)
         if cloc ne 1 then stop,'oops'
         ind = array_indices(grid_feh,loc)
         chisqarr[ind[0],ind[1],ind[2]] = total((refspec-specarr[i].spec)^2/refspec_err)
      endfor
   endfor
   feharr = grid_feh[*,0,0]
   agearr = reform(grid_age[0,*,0])
   alphaarr = reform(grid_alpha[0,0,*])
   deltafeh = grid_feh[1,0,0]-grid_feh[0,0,0]
   deltaage = grid_age[0,1,0]-grid_age[0,0,0]
   deltaalpha = grid_alpha[0,0,1]-grid_alpha[0,0,0]
   ;read the chisq grid and calculate probability
   Lgrid = -0.5*(chisqarr-min(chisqarr))
   probgrid = exp(double(Lgrid))
   volume = 0
   arr_dimen = size(grid_feh,/dimension)
   for ii=0,arr_dimen(0)-1 do for jj=0,arr_dimen(1)-1 do for kk=0,arr_dimen(2)-1 do begin
         volume = volume+deltafeh*deltaage*deltaalpha*probgrid[ii,jj,kk]
   endfor
   probgrid = probgrid/volume

   probfeh = fltarr(arr_dimen(0))
   probage = fltarr(arr_dimen(1))
   probalpha = fltarr(arr_dimen(2))
   for ii=0,arr_dimen(0)-1 do probfeh(ii) = int_tabulated_2d(grid_age[ii,*,*],grid_alpha[ii,*,*],probgrid[ii,*,*])
   for jj=0,arr_dimen(1)-1 do probage(jj) = int_tabulated_2d(grid_feh[*,jj,*],grid_alpha[*,jj,*],probgrid[*,jj,*])
   for kk=0,arr_dimen(2)-1 do probalpha(kk) = int_tabulated_2d(grid_age[*,*,kk],grid_feh[*,*,kk],probgrid[*,*,kk])
   probfeh = probfeh/int_tabulated(feharr,probfeh)
   probage = probage/int_tabulated(agearr,probage)
   probalpha = probalpha/int_tabulated(alphaarr,probalpha)

   cumprobfeh = fltarr(arr_dimen(0))
   for jj=1,arr_dimen(0)-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
   cumprobage = fltarr(arr_dimen(1))
   for jj=1,arr_dimen(1)-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])
   cumprobalpha = fltarr(arr_dimen(2))
   for kk=1,arr_dimen(2)-1 do cumprobalpha(kk) = int_tabulated(alphaarr[0:kk],probalpha[0:kk])
   midagevalue = interpol(agearr,cumprobage,0.5)
   midfehvalue = interpol(feharr,cumprobfeh,0.5)
   midalphavalue = interpol(alphaarr,cumprobalpha,0.5)

   science.ageupper = science.age+(interpol(agearr,cumprobage,0.84)-midagevalue)
   science.agelower = science.age+(interpol(agearr,cumprobage,0.16)-midagevalue)
   science.fehupper = science.feh+(interpol(feharr,cumprobfeh,0.84)-midfehvalue)
   science.fehlower = science.feh+(interpol(feharr,cumprobfeh,0.16)-midfehvalue)
   science.alphafeupper = science.alphafe+(interpol(alphaarr,cumprobalpha,0.84)-midalphavalue)
   science.alphafelower = science.alphafe+(interpol(alphaarr,cumprobalpha,0.16)-midalphavalue)
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready ...'
end


pro sps_fit::fit_all,alpha=alpha,omg=omg,glorious=glorious,mgfeo=mgfeo
    common mask_in, mask_in, copynum
    widget_control, widget_info(self.base, find_by_uname='keepoldfit'), get_value=keepoldfit
    scienceall = *self.science
    curi = self.i
 for nwalkers=0,0 do begin
    nreplace = 0
    for i=0,self.nspec-1 do begin
        self.i = i
        self->default_range
        science = scienceall[self.i]
        if science.good eq 0 then continue
        if ~keyword_set(alpha) and ~keyword_set(omg) and ~keyword_Set(glorious) then self->fit, science
        if keyword_set(alpha) then begin
           ;self->mask, science ,/includemg
           self->fitalpha, science
        endif
        if keyword_set(omg) then begin
           self->fitomg, science
        endif
        if keyword_set(glorious) then begin
           self->fit_glorious, science
        endif
        if keyword_set(mgfeo) then begin
           self->fit_mgfeo, science
        endif
        if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
            print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
            scienceall[self.i] = science
            if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
            if keepoldfit eq 0 then printf,long(copynum),'chisq is smaller, replaced fit'
            nreplace += 1
        endif else begin
            if keepoldfit eq 0 then begin
            print, 'chisq is larger, use previous fit'
            printf,long(copynum),'chisq is larger, use previous fit'
            print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
            endif
        endelse
        close,copynum
        self->statusbox
        if i mod 30 eq 0 then begin
            ptr_free, self.science
            self.science = ptr_new(scienceall)
            self->writescience
            scienceall = *self.science
        endif
    endfor
    print,'nwalker',nwalkers,'total replace ', nreplace,' fits'
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
  endfor
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    self->redraw    
end

pro sps_fit::cal_uncertainties_all,alpha=alpha
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Calculatign uncertainties ...'
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        science = scienceall[self.i]
        if science.goodfit eq 0 then continue
        print,strtrim(string(i+1),2)+'/'+strtrim(string(self.nspec),2)
        if ~keyword_set(alpha) then self->cal_uncertainties, science else $
                                    self->cal_uncertainties_alpha,science
        scienceall[self.i] = science
    endfor
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready'
end


; =================
pro sps_fit_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)

    case (name) of
        'BUTT': obj->handle_button, ev
        'TEXT': obj->handle_text, ev
        'COMB': obj->handle_combobox, ev
        'DRAW': begin
            if ev.type eq 0 then obj->handle_draw_click, ev
            if ev.type eq 5 then obj->handle_draw_key, ev
        end
        'DONE': widget_control, ev.top, /destroy
        else: begin
            case (value) of
                'good': obj->toggle_good
                else: obj->redraw
            endcase
            widget_control, widget_info(ev.top, find_by_uname='spec'), /input_focus
        end
    endcase
end


; ================ BUTTONS ================
function line_event, ev
    self->redraw
end


pro sps_fit::handle_button, ev
    common mask_in, mask_in, copynum
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    
    case (uvalue) of
        'back': self->newspec, increment=-1
        'next': self->newspec, increment=1
        'default_range': begin
            self->default_range
            self->redraw
        end
        'reprepare': self->reprepare
        'fit': begin
            scienceall = *self.science
            science = scienceall[self.i]
            self->fit, science, /noredraw
            widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,$
                      scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,$
                      format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               scienceall[self.i] = science
               if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
           endif else begin
               if keepoldfit eq 0 then begin
               print, 'chisq is larger, use previous fit'
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               endif
           endelse
            ptr_free, self.science
            self.science = ptr_new(scienceall)
            self->redraw
            self->statusbox
         end
       'fitalpha':begin
           scienceall = *self.science
           ;for nw=0,5 do begin
           science = scienceall[self.i]
           ;self->mask, science ,/includemg
           self->fitalpha, science, /noredraw
           widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,$
                      scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,$
                      format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               scienceall[self.i] = science
               if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
               if keepoldfit eq 0 then printf,copynum,'chisq is smaller, replaced fit'
           endif else begin
               if keepoldfit eq 0 then begin
               print, 'chisq is larger, use previous fit'
               printf,copynum,'chisq is larger, use previous fit'
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               endif
           endelse
           close,copynum
           ;endfor
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
           self->statusbox
       end
       'fitomg':begin
           scienceall = *self.science
           ;for nw=0,5 do begin
           science = scienceall[self.i]
           ;self->mask, science ,/includemg
           self->fitomg, science, /noredraw
           widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,$
                      scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,$
                      format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               scienceall[self.i] = science
               if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
               if keepoldfit eq 0 then printf,copynum,'chisq is smaller, replaced fit'
           endif else begin
               if keepoldfit eq 0 then begin
               print, 'chisq is larger, use previous fit'
               printf,copynum,'chisq is larger, use previous fit'
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               endif
           endelse
           close,copynum
           ;endfor
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
           self->statusbox
       end
       'fitmgfeo':begin
           scienceall = *self.science
           ;for nw=0,5 do begin
           science = scienceall[self.i]
           ;self->mask, science ,/includemg
           self->fitomg, science, /noredraw
           widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,$
                      scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,$
                      format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               scienceall[self.i] = science
               if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
               if keepoldfit eq 0 then printf,copynum,'chisq is smaller, replaced fit'
           endif else begin
               if keepoldfit eq 0 then begin
               print, 'chisq is larger, use previous fit'
               printf,copynum,'chisq is larger, use previous fit'
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               endif
           endelse
           close,copynum
           ;endfor
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
           self->statusbox
       end

       'fitglorious':begin
           scienceall = *self.science
           ;for nw=0,5 do begin
           science = scienceall[self.i]
           ;self->mask, science ,/includemg
           self->fit_glorious, science, /noredraw
           widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,$
                      scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,$
                      format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               scienceall[self.i] = science
               if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
               if keepoldfit eq 0 then printf,copynum,'chisq is smaller, replaced fit'
           endif else begin
               if keepoldfit eq 0 then begin
               print, 'chisq is larger, use previous fit'
               printf,copynum,'chisq is larger, use previous fit'
               print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
               endif
           endelse
           close,copynum
           ;endfor
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
           self->statusbox
       end

        'fitalphaband':begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->fitalphaband, science, /noredraw
           scienceall[self.i] = science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
           self->statusbox
       end
        'cal_uncertainties': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->cal_uncertainties,science
           scienceall[self.i]=science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->statusbox, science=science
        end
        'cal_uncertainties_alpha': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->cal_uncertainties_alpha,science
           scienceall[self.i]=science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->statusbox, science=science
        end
        'reprepare_all': self->reprepare_all
        'fit_all': self->fit_all
        'fit_alpha_all':self->fit_all,/alpha
        'fit_glorious_all':self->fit_all,/glorious
        'fit_omg_all':self->fit_all,/omg
        'cal_uncertainties_all':self->cal_uncertainties_all
        'cal_uncertainties_alpha_all':self->cal_uncertainties_all,/alpha
        'default_cont': self->default_cont
        'default_mask': self->default_mask
        'default_maskall': self->default_maskall
        'default_indices_mask': self->default_indices_mask
        'default_goodspec':self->default_goodspec
        'backward': self->step, -1
        'forward': self->step, 1
        'good': self->toggle_good
        'blue': self->lambdarange, /blue
        'red': self->lambdarange, /red
        'save': self->writescience
        'exit': begin
            self->writescience
            widget_control, ev.top, /destroy
        end
        else:
    endcase
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


pro sps_fit::step, increment
    nmodes = 4
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 3 and increment gt 0: mode = 0
        mode eq 0 and increment lt 0: mode = 3
        else: mode += increment
    endcase
    curi = self.i
    if mode eq -1 then begin
        self->newspec, increment=-1, /noredraw
        if self.i eq curi then return
        mode = nmodes-1
        ;self->default_range
    endif
    if mode eq nmodes then begin
        self->newspec, increment=1, /noredraw
        if self.i eq curi then return
        mode = 1
        ;self->default_range
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), set_value=mode
    self.keystate = 0
    self->redraw
end


pro sps_fit::toggle_good
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
    scienceall = *self.science
    widget_control, widget_info(self.base, find_by_uname='good'), get_value=good
    scienceall[self.i].good = good[0]
    scienceall[self.i].goodfit = good[1]
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::output_objname
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Outputting object name to output.dat ...'
    comment = textbox(title='Provide a Comment (optional)', group_leader=self.base, label='Comment: ', cancel=cancelled, xsize=59, value='')
    if ~cancelled then begin
        science = (*self.science)[self.i]
        openw, lun, 'output.dat', /append, /get_lun
        printf, lun, science.objname, comment, format='(A-20,1X,A-59)'
        close, lun
        free_lun, lun
    endif
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::newspec, increment=increment, noredraw=noredraw
    newi = self.i
    newi += increment
    if newi lt 0 then begin
       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the first spectrum.'
       return
    endif
    if newi gt self.nspec-1 then begin
       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the last spectrum.'
       return
    endif

    widget_control, widget_info(self.base, find_by_uname='filelist'), set_combobox_select=newi
    self.i = newi
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec,/nan)
    ;self->default_range
    
    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw

end


pro sps_fit::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.science)[self.i].spec,/nan)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self.lambdalim = (minmax((*self.science)[self.i].lambda / (1d + (*self.science)[self.i].zspec),/nan) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d + (*self.science)[self.i].zspec)
     endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        mode eq 2: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.skyylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.skyylim[1], format='(D5.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
     endcase
    
    widget_control, widget_info(self.base, find_by_uname='npoly'), set_value=strcompress(string((*self.science)[self.i].npoly, format='(I2)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
end


pro sps_fit::lambdarange, red=red, blue=blue
    if keyword_set(red)+keyword_set(blue) ne 1 then message, 'You must specify red or blue.'
    lrange = self.lambdalim[1] - self.lambdalim[0]
    if keyword_set(blue) then begin
        if self.lambdalim[0] lt min((*self.science)[self.i].lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is bluest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] - 0.6*lrange
        lhighnew = self.lambdalim[1] - 0.6*lrange
    endif
    if keyword_set(red) then begin
        if self.lambdalim[1] gt max((*self.science)[self.i].lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is reddest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] + 0.6*lrange
        lhighnew = self.lambdalim[1] + 0.6*lrange
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 4 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    self->lambdalow, llownew, /noredraw
    self->lambdahigh, lhighnew
end


; ============== TEXT BOXES =============
pro sps_fit::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
       'npoly':self->npoly,val
       'ylow': self->ylow, val
       'yhigh': self->yhigh, val
       'lambdalow': self->lambdalow, val
       'lambdahigh': self->lambdahigh, val
       else: 
    end
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro sps_fit::lambdalow, lambdalow, noredraw=noredraw
    ;if lambdalow gt max((*self.science)[self.i].lambda) then begin
    ;    lambdalow = max((*self.science)[self.i].lambda) - (self.lambdalim[1] - lambdalow)
    ;    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    ;    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    ;    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    ;endif
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::lambdahigh, lambdahigh, noredraw=noredraw
    ;if lambdahigh lt min((*self.science)[self.i].lambda) then begin
    ;    lambdahigh = min((*self.science)[self.i].lambda) + (lambdahigh - self.lambdalim[0])
    ;    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    ;    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    ;    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    ;endif
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->redraw
end

pro sps_fit::npoly, npoly, noredraw=noredraw
    scienceall = *self.science
    scienceall[self.i].npoly = npoly
    self->statusbox, science=scienceall[self.i]
    ptr_free, self.science
    self.science = ptr_new(scienceall)

end


pro sps_fit::ylow, ylow, noredraw=noredraw
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[0] = ylow
        mode eq 2: self.skyylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[1] = yhigh
        mode eq 2: self.skyylim[1] = yhigh
        else: self.divylim[1] = yhigh
    endcase
    self->redraw
end


pro sps_fit::smooth, smooth    
    leftover = smooth mod 2
    smooth = fix(smooth)
    if leftover ge 1 then smooth += 1
    smooth >= 0
    self.smooth = smooth
    widget_control, widget_info(self.base, find_by_uname='smooth'), set_value=strcompress(self.smooth, /rem)
    self->redraw
end

; ============== COMBOBOX  =============
pro sps_fit::handle_combobox, ev
    ;widget_control, ev.id, get_uvalue=uvalue
    ;widget_control, ev.top, get_uvalue=obj
    self.i = ev.index
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec,/nan)
    self.keystate = 0
    self->redraw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; ============= DRAW CLICK =============
pro sps_fit::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 3 then click_coords /= 1d + (*self.science)[self.i].zspec
    case ev.modifiers of
        0: begin
            case ev.press of
                1: lrange = (self.lambdalim[1] - self.lambdalim[0]) / 2.
                4: lrange = (self.lambdalim[1] - self.lambdalim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.
            zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        1: begin
            if ev.press ne 1 then begin
                widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                return
            endif
            lrange = self.lambdalim[1] - self.lambdalim[0]
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.        
            zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 0: ylim = self.ylim
                mode eq 2: ylim = self.skyylim
                else: ylim = self.divylim
            endcase
            case ev.press of
                1: yrange = (ylim[1] - ylim[0]) / 2.
                4: yrange = (ylim[1] - ylim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            ylownew = click_coords[1] - yrange/2.
            yhighnew = click_coords[1] + yrange/2.
            case 1 of
                mode le 0: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(g8.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(g8.1)'), /rem)
                end
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(D5.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(D5.1)'), /rem)
                end
            endcase
            self->ylow, ylownew, /noredraw
            self->yhigh, yhighnew
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
            return            
        end
    endcase
end


; ============= DRAW KEYS ==============
pro sps_fit::handle_draw_key, ev
    if ev.press ne 1 then return
    key = string(ev.ch)
    coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    science = (*self.science)[self.i]
    case mode of
  
        0: begin        ;continuum fit
            case key of
                'g': begin
                    case self.keystate of
                        0: begin
                            widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                            self->toggle_good
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'b': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=-1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'n': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'f': begin
                    case self.keystate of
                        0: begin
                            self->default_range
                            self->redraw
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'q': begin
                    case 1 of
                        self.keystate eq 1 or self.keystate eq 2: begin
                            self.keystate = 0
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Continuum mask modification cancelled.'
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'e': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 1
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'e' again to exclude wavelength region from continuum mask."
                        end
                        1: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda gt lambda1 and scienceall[self.i].lambda lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].contmask[w] = 0
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'i': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 2
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'i' again to include wavelength region in continuum mask."
                        end
                        2: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda gt lambda1 and scienceall[self.i].lambda lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].contmask[w] = 1
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase

                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        1: begin        ;continuum division
            case key of
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        2: begin        ;sky line fit
            case key of
                's': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[~science.goodsky, science.good]
                    self->toggle_good
                end
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'a': self->add_skyline, coords[0]
                'd': self->delete_skyline, coords[0]
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        3: begin        ;pixel mask
            case key of
                'g': begin
                    case self.keystate of
                        0: begin
                            widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                            self->toggle_good
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'b': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=-1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'n': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'f': begin
                    case self.keystate of
                        0: begin
                            self->default_range
                            self->redraw
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'h': begin
                    self.lambdalim = [6513, 6613]
                    self->redraw
                end
                'q': begin
                    case 1 of
                        self.keystate eq 1 or self.keystate eq 2: begin
                            self.keystate = 0
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Pixel mask modification cancelled.'
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'e': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 1
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'e' again to exclude wavelength region from pixel mask."
                        end
                        1: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].fitmask[w] = 0
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'i': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 2
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'i' again to include wavelength region in pixel mask."
                        end
                        2: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].fitmask[w] = 1
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase

                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end
    endcase
end


; ============== DESTROY ===============
pro sps_fit_cleanup, ev    
    widget_control, ev, get_uvalue=obj
    obj_destroy, obj
end


pro sps_fit::cleanup
    ptr_free, self.science
end


pro sps_fit::reprepare_all
    update_phot = 1

    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        if ~update_phot and (*self.science)[self.i].good eq 0 then continue
        self->reprepare, /nostatusbar
    endfor
    scienceall = *self.science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = (*self.science)[self.i]
    self->statusbox, science=science
    self->redraw
end


pro sps_fit::default_cont
    scienceall = *self.science
    science = scienceall[self.i]
    science.contmask = 0
    self->continuum, science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end


pro sps_fit::default_mask
    scienceall = *self.science
    science = scienceall[self.i]
    self->mask, science, nmc=0
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
 end

pro sps_fit::default_maskall
  curi = self.i
  widget_control, widget_info(self.base, find_by_uname='status'), set_value='Defaulting masks'
  for i=0,self.nspec-1 do begin
     self.i=i
     scienceall = *self.science
     science = scienceall[self.i]
     self->mask, science;,/includemg
     scienceall[self.i] = science
     ptr_free, self.science
     self.science = ptr_new(scienceall)
  endfor
  self.i = curi
  widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
  self->redraw
end


pro sps_fit::default_goodspec
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
       self.i=i
       science = scienceall[self.i]
       sn = science.sn
       z  = science.zspec
       if sn gt 6. and z gt 0.0 then science.good = 1 else science.good = 0
       scienceall[self.i] = science
    endfor
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    self->redraw

end



pro sps_fit::reprepare, nostatusbar=nostatusbar
    common mask_in, mask_in, copynum
    update_else = 1

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Repreparing ...'
    scienceall = *self.science
    science = scienceall[self.i]
    contmask = science.contmask

    n = n_elements(science.lambda)
    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.spec[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])

    if update_else eq 1 then begin
        ;self->skytweak, science
        ;science.skylinemask = -1
        ;science.contmask = 0
        self->specres, science, /goodoverride
        self->continuum, science
        self->mask, science, /nomask
        science.spscont = 1.0
        self->indices, science, /noredraw
        ;self->fit, science, /noredraw, nostatusbar=nostatusbar
    endif

    self->statusbox, science=science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end


; ============== SPECRES ==============
pro fitskyspec_gauss, x, a, f, pder
    u = (x - a[0])/a[1]
    f = a[2]*exp(-0.5*u^2)
    pder = [[(a[0]-x)/a[1]*f], [(x-a[0])^2./a[1]^3.*f], [f/a[2]]]
    return
end


pro sps_fit::skytweak, science
    blue = mrdfits(strtrim(science.spec1dfile, 2), 1, headblue, /silent)
    red = mrdfits(strtrim(science.spec1dfile, 2), 2, headred, /silent)
    lambdablue = blue.lambda
    lambdared = red.lambda
    path = file_dirname(science.spec1dfile)
    lambdabluenew = applytweaks(lambdablue, headblue, path)
    lambdarednew = applytweaks(lambdared, headred, path)
    lambdanew = [lambdabluenew, lambdarednew]
    science.lambda = lambdanew
end


function sps_fit::fitskylines, science
    lambda = science.lambda
    skyspec = science.skyspec
    w = where(~finite(skyspec), c)
    if c gt 0 then skyspec[w] = 0
    skyspec = skyspec/max(skyspec)
    medskyspec = median(skyspec)

    deriv1skyspec = deriv(lambda, skyspec)
    deriv2skyspec = deriv(lambda, deriv1skyspec)

    thresh = 1.
    nlines = 1000
    while nlines gt 200 do begin
        w = where(abs(deriv1skyspec) lt 0.2 and deriv2skyspec lt -0.005 and skyspec gt thresh*medskyspec)
        w = [w, n_elements(lambda)+1]
        wstep = round(-1*ts_diff(w, 1))
        linestart = where(wstep gt 1, nlines)
        if nlines lt 5 then begin
            message, 'Fewer than 5 sky lines in this spectrum.', /info
            return, [[-1], [-1], [-1]]
        endif
        linepix = round(-1*ts_diff(linestart, 1))
        nlines -= 1
        thresh *= 2
     endwhile
    linepix = linepix[0:nlines-1]
    wlocmax = lindgen(nlines)
    sigma = dblarr(nlines)
    sigmaerr = dblarr(nlines)
    linelambda = dblarr(nlines)
    nskyspec = n_elements(skyspec)
    wgoodline = bytarr(nlines)+1
    for i=0,nlines-1 do begin        
        if linepix[i] eq 1 then wlocmax[i] = w[linestart[i]] else begin
            junk = min(abs(deriv1skyspec[w[linestart[i]]:w[linestart[i]]+linepix[i]-1]), wmin)
            wlocmax[i] = w[linestart[i]]+wmin
        endelse
        if wlocmax[i]-10 lt 0 or wlocmax[i]+10 ge nskyspec then begin
            wgoodline[i] = 0
            continue
        endif
        skyspecmax = max(skyspec[wlocmax[i]-10:wlocmax[i]+10], wlocmaxtemp)
        wlocmax[i] += wlocmaxtemp-10
        wfit = lindgen(20)+wlocmax[i]-10
        wfit = wfit[where(wfit ge 0 and wfit lt n_elements(lambda), nfit)]
        lambdafit = lambda[wfit]
        skyspecfit = skyspec[wfit]

        ;lambdarange = max(lambda[wfit])-min(lambda[wfit])
        ;lambdastep = 0.03
        ;lambdafit = dindgen(round(lambdarange/lambdastep))*lambdastep+min(lambda[wfit])
        ;skyspecfit = interpol(skyspec[wfit], lambda[wfit], lambdafit, /quadratic)

        skyspecmax = max(skyspecfit, wmax)

        ;a = [lambdafit[wmax], 1.3/2.35, skyspecmax]
        ;yfit = curvefit(lambdafit, skyspecfit, lindgen(n_elements(lambdafit))+1.0, a, aerr, function_name='fitskyspec_gauss', chisq=chisq, status=status, /double, tol=1d-12)
        ;a = [skyspecmax, lambdafit[wmax], 1.3/2.35]

        guess = [skyspecmax, lambdafit[wmax], 0, medskyspec]
        yfit = gaussfit(lambdafit, skyspecfit, a, estimates=guess, sigma=aerr, nterms=4, chisq=chisq)
        sigma[i] = abs(a[2])
        sigmaerr[i] = aerr[2]
        linelambda[i] = a[1]

        ;print, sigma[i], sigmaerr[i]/sigma[i], chisq
        ;plot, lambdafit, skyspecfit
        ;oplot, lambdafit, yfit, color=fsc_color('red')
        ;wait, 0.5

        if chisq gt 1d-3 or sigmaerr[i] gt 0.8 then wgoodline[i] = 0
    endfor
    wgood = where(sigma gt 0 and wgoodline eq 1, c)
    if c gt 0 then begin
        linelambda = linelambda[wgood]
        sigma = sigma[wgood]
        sigmaerr = sigmaerr[wgood]
        ;ploterr, linelambda, sigma, sigmaerr, psym=1
        
        ;plot, lambda, skyspec
        ;oplot, lambda[wlocmax], skyspec[wlocmax], psym=1, color=fsc_color('red')

        return, [[linelambda], [sigma], [sigmaerr]]
    endif else return, [[-1], [-1], [-1]]
end


pro sps_fit::specres, science, qf=qf, goodoverride=goodoverride
    common mask_in, mask_in, copynum
    lambda = science.lambda
    if mask_in eq 'NGC6528' then science.dlam = 3.1/2.35 else science.dlam = science.lambda/1000. ;R=1000

end


; ============= CONTINUUM =============
pro sps_fit::continuum, science
  common mask_in, mask_in, copynum
    fft = 0
    spline = 0
    poly = 0
    usesmooth = 0
    nofit = 1
    simpoly=0
   ; if mask_in eq 'NGC6528' then spline = 1 else nofit =1

    contmask = science.contmask
    w = where(contmask ne 0, c)
    lambda = science.lambda / (1d + science.zspec)
    n = n_elements(lambda)
    if c eq 0 then begin
        linestart = *self.linestart
        lineend = *self.lineend
        contmask = bytarr(n_elements(lambda))+1
        for i=0,n_elements(linestart)-1 do begin
            w = where(lambda ge linestart[i] and lambda le lineend[i], c)
            if c gt 0 then contmask[w] = 0
        endfor
        tellmask = bytarr(n_elements(lambda))
        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 100000.]
        for i=0,n_elements(tellstart)-1 do begin
            w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
            if c gt 0 then begin
                contmask[w] = 0
                tellmask[w] = 1
            endif
        endfor
        contmask[0:2] = 0
        contmask[n-3:n-1] = 0
    endif

    satbands = [6870, 7650]
    niter = 5

    wwhole = lindgen(n)
    spec = science.spec
    contmask(where(finite(spec) eq 0)) = 0
    contmask[wwhole[0:3]] = 0
    won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
    woff += wwhole[0]
    case 1 of
       spline: begin
          bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=science.ivar[won], bkspace=165, upper=5, lower=1.5,mask=mask, /silent,/everyn)
          if bkpt[0] eq -1 then return
          cont = slatec_bvalu(lambda[wwhole], bkpt, coeff)
       end
       poly: begin
          degree = 12
          norm = median(spec[won])
          a = [norm, replicate(0.0, degree-1)]
          p = lmfit(lambda[won], spec[won], a, measure_errors=(science.ivar[won])^(-0.5), /double, function_name='legendre_poly')
          cont = legendre_poly(lambda[wwhole], a, /noderiv)
       end
       fft: begin
          wfft = wwhole[10:nhalf-11]
          nfft = n_elements(wfft)
          nkeep = 10
          ft = fft(spec[wfft], -1, /double)
          
          ft[nkeep:nfft-1] = 0
          cont = real_part(fft(ft, 1, /double))
          plot, cont
          stop
       end
       usesmooth: begin
          ww = wwhole[3:nhalf-4]
          wcont = where(contmask[ww] eq 1, ccont) + ww[0]
          for i=1,5 do begin
             wcontold = wcont
             if ccont lt 50 then begin
                science.continuum = replicate(-999, n_elements(science.continuum))
                return
             endif
             cont = smooth_gauss_wrapper(lambda[wcont], spec[wcont], lambda[wwhole], 30.0, ivar1=science.ivar[wcont])
             wcont = where(abs(spec[ww]-cont[ww-wwhole[0]]) lt (science.ivar[ww])^(-0.5) and contmask[ww] eq 1, ccont) + ww[0]
             
             if array_equal(wcont, wcontold) then break
          endfor
       end
       nofit:begin
          cont =  lambda*0.+1.
          ;contmask = bytarr(n)+1
       end
       simpoly:begin
          degree=7
          p=poly_fit(lambda[won],spec[won],degree,measure_errors=(science.ivar[won])^(-0.5))
          cont = poly(lambda[wwhole],p)
       end
    endcase
    
    science.contmask = contmask
    science.continuum = cont
    science.contdiv = science.spec/cont
    science.contdivivar = science.ivar*cont^2.

    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.spec[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])
end


; =============== MASK ================
pro sps_fit::mask, science, nomask=nomask, zfind=zfind, nozfind=nozfind, nmc=nmc,includemg=includemg
    if ~keyword_set(nomask) then begin
      ;1) remove emission lines
       linestart = *self.linestart
       lineend = *self.lineend
       linetype = *self.linetype
       wem = where(linetype eq 'e', cemlines)
       mask = bytarr(n_elements(science.lambda))+1
       w = where(science.ivar le 0 or ~finite(science.ivar) or science.lambda/(1d + science.zspec) lt 3650, cw)
       if cw gt 0 then mask[w] = 0
       for i=0,cemlines-1 do begin
          w = where(science.lambda/(1d + science.zspec) ge linestart[wem[i]] and science.lambda/(1d + science.zspec) le lineend[wem[i]], c)
          if c gt 0 then mask[w] = 0
       endfor

       
       ;;3) where it peaks greater than 3 sigmas
       ;w = where((science.contdiv-science.continuum) gt 3.*science.contdivivar^(-0.5), cw)
       ;if cw gt 0 then mask[w] = 0
        
        ; 5) extreme noises
       w = where(science.contdiv gt 2. or science.contdiv lt -0.5, cw)
       nlambda = n_elements(science.lambda)
       if cw gt 0 then for i=0,cw-1 do mask[w[i]-2>0:w[i]+2<nlambda-1]=0
    
        ;6) where it's not finite
       w = where(finite(science.spec) eq 0 or finite(science.lambda) eq 0,c)
       if c gt 0 then mask[w]=0
       ;;

        ;7)where is is Mg region
        if ~keyword_set(includemg) then begin
           linestart = *self.linestart
           lineend  = *self.lineend
           linetype = *self.linetype
           indstart = *self.indstart 
           indend   = *self.indend   
           indname  = *self.indname 
   
           use_indices = ['Mg_b','Mg_2','Mg_1']
           for i=0,n_elements(use_indices)-1 do begin
              indnow = where(indname eq use_indices(i),cindnow)
              if cindnow eq 0 then stop
              w = where(science.lambda/(1.+science.zspec) gt indstart(indnow[0]) and science.lambda/(1.+science.zspec) lt indend(indnow[0]),cw)
              if cw gt 0 then mask[w]=0
           endfor
        endif
      ;8) where contmask is zero
       w = where(science.contmask eq 0, c)
       if c gt 0 then mask[w]=0

        ;8)wavelength range
        w = where(science.lambda/(1.+science.zspec) lt 3650. or science.lambda/(1.+science.zspec) gt 6000.,c)
        if c gt 0 then mask[w]=0


        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       
     endif else begin
        mask = science.fitmask
     endelse
     science.fitmask = mask

end

; ============== REDRAW ===============
pro sps_fit::redraw
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Redrawing ...'
    self->default_range, /update
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
    wset, index
    science = (*self.science)[self.i]

    case mode of
        -1: begin                ;smoothed
            specsmooth = smooth(science.spec, 100)
            plot, science.lambda, specsmooth, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6100-pixel boxcar smoothed flux (e!E-!N/hr)!3'
            yoff1 = 0.01*(self.ylim[1] - self.ylim[0])
            yoff2 = 0.05*(self.ylim[1] - self.ylim[0])
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + science.zspec) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + science.zspec) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + science.zspec), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + science.zspec))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end

        0: begin                ;continuum fit
            t = round(-1*ts_diff(science.contmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.contmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.contmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending continuum wavelengths.'
            if cstart eq 0 or cend eq 0 then message, 'There are no continuum regions.'

            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
            for i=0,cstart-1 do begin
                x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + science.zspec))) < (self.lambdalim[1] * (1d + science.zspec))
                y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
                polyfill, x, y, color=fsc_color('light cyan')
            endfor
            oplot, science.lambda, science.spec, color=fsc_color('black')
            oplot, science.lambda, science.continuum, color=fsc_color('green')
            if science.feh lt 3. then begin
               normfactor = median(science.spec)/median(science.spsspecfull)
               oplot, science.lambda,science.spsspecfull*normfactor,color=fsc_color('red')            
               oplot, science.lambda,science.spscontfull*normfactor,color=fsc_color('orange')            

            endif
            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + science.zspec) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + science.zspec) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + science.zspec), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + science.zspec))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end


        1: begin                ;continuum division
            plot, science.lambda, science.contdiv, xrange=self.lambdalim * (1d + science.zspec), yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata
            oplot, self.lambdalim, [1.0, 1.0], color=fsc_color('orange')
            oplot, science.lambda, science.contdiv, color=fsc_color('black')
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
        end

        2: begin        ;sky line fit
            plot, science.lambda, 2.35*science.dlam, xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6sky line FWHM (!sA!r!u!9 %!6!n)!3'
        end

        3: begin        ;rest frame
           ;;Highlight the fit mask regions
            t = round(-1*ts_diff(science.fitmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.fitmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.fitmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
             endif
            
            if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
            ;;make reference axis of contdiv with applied fake continuum
            plot, science.lambda/(1d + science.zfit), science.contdiv/science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
            if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
                for i=0,cstart-1 do begin
                    x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]]/(1d + science.zspec) > self.lambdalim[0]) < self.lambdalim[1]
                    y = [self.divylim[0], self.divylim[1], self.divylim[1], self.divylim[0]]
                    polyfill, x, y, color=fsc_color('light cyan')
                endfor
             endelse
           ;;finish highligting fit mask
            oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
            oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
            if science.zfit ne 0 and finite(science.zfit) then znow=science.zfit else znow = science.zspec
            ;;plot the contdiv with applied fake continuum
            oplot, science.lambda/(1d +znow), science.contdiv/science.spscont, color=fsc_color('black') 
            ;;highlighting the telluric bands and label line names
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]] / (1d + science.zspec), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
               if (*self.linewaves)[i] le !X.CRANGE[0] or (*self.linewaves)[i] ge !X.CRANGE[1] then continue
               oplot, [(*self.linewaves)[i], (*self.linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
               xyouts, (*self.linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
            ;;plot the models
            plot, science.lambda/(1d +znow), science.contdiv/science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6rest wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata, /noerase
            oplot,science.lambda/(1d + science.zfit), (1./sqrt(science.contdivivar))/science.spscont+0.2,color=fsc_color('gray')
            if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) then begin
               oplot, science.lambda/(1d + znow), science.spsspecfull, color=fsc_color('red')
            endif
            wmg = where(science.lambda/(1d + znow) gt 5065. and science.lambda/(1d + znow) lt 5250)
            oplot,science.lambda(wmg)/(1d + znow), science.spsspecfullmg(wmg), color=fsc_color('purple')
            ;;plot choi model
            if znow gt 0.3 and znow lt 0.4 then zreal = 0.35
            if znow gt 0.4 and znow lt 0.5 then zreal = 0.475
            if znow gt 0.6 and znow lt 0.7 then zreal = 0.625
           ; oplot,science.lambda/(1d + zreal), science.spsspecfullchoi, color=fsc_color('blue')
         end
    endcase        
    zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::statusbox, science=science
    if ~keyword_set(science) then science = (*self.science)[self.i]
    
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.good, science.goodfit]
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=science.zrange+science.massbin+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(science.zspec, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzfit'), set_value=strcompress(string(science.zfit, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curchisq'), set_value=strcompress(string(science.chisq, format='(D8.5)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curnloop'), set_value=strcompress(string(science.nloop, format='(I4)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curnpoly'), set_value=strcompress(string(science.npoly, format='(D4.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=science.sn gt 0 ? strcompress(string(science.sn, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curage'), set_value=science.age gt -100 ? strcompress(string(science.age, format='(D10.2)'), /rem)+(science.ageerr le 0 ? '' : ' +/- '+strcompress(string(science.ageerr, format='(D10.2)'), /rem))+' Gyr' : unknown
    widget_control, widget_info(self.base, find_by_uname='curageuncert'), set_value=science.agelower gt -100 ? strcompress(string(science.agelower, format='(D10.2)'), /rem)+(science.ageupper gt -100 ? ' : '+strcompress(string(science.ageupper, format='(D10.2)'), /rem):'') : unknown
    widget_control, widget_info(self.base, find_by_uname='curagechoi'), set_value=science.agechoi gt -100 ? strcompress(string(science.agechoi, format='(D10.2)'), /rem)+(science.agechoierr le 0 ? '' : ' +/- '+strcompress(string(science.agechoierr, format='(D10.2)'), /rem))+' Gyr' : unknown

    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=science.logmstar gt 0 ? strcompress(string(science.logmstar, format='(D10.2)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfeh'), set_value=science.feh gt -100 ? strcompress(string(science.feh, format='(D10.2)'), /rem)+(science.feherr le 0 ? '' : ' +/- '+strcompress(string(science.feherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfehuncert'), set_value=science.fehlower gt -100 ? strcompress(string(science.fehlower, format='(D10.2)'), /rem)+(science.fehupper ge -100 ? ' : '+strcompress(string(science.fehupper, format='(D10.2)'), /rem):'') : unknown

    widget_control, widget_info(self.base, find_by_uname='curmg'), set_value=science.alphafe gt -100 ? strcompress(string(science.alphafe, format='(D10.2)'), /rem)+(science.alphafeerr le 0 ? '' : ' +/- '+strcompress(string(science.alphafeerr, format='(D10.2)'), /rem)) : unknown
widget_control, widget_info(self.base, find_by_uname='curmguncert'), set_value=science.alphafelower gt -100 ? strcompress(string(science.alphafelower, format='(D10.2)'), /rem)+(science.alphafeupper ge -100 ? ' : '+strcompress(string(science.alphafeupper, format='(D10.2)'), /rem):'') : unknown

    widget_control, widget_info(self.base, find_by_uname='curfehchoi'), set_value=science.fehchoi gt -100 ? strcompress(string(science.fehchoi, format='(D10.2)'), /rem)+(science.fehchoierr le 0 ? '' : ' +/- '+strcompress(string(science.fehchoierr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curmgfechoi'), set_value=science.mgfechoi gt -100 ? strcompress(string(science.mgfechoi, format='(D10.2)'), /rem)+(science.mgfechoierr le 0 ? '' : ' +/- '+strcompress(string(science.mgfechoierr, format='(D10.2)'), /rem)) : unknown

    widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=strcompress(string(science.vdisp, format='(D10.1)'), /rem)+(science.vdisperr le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr, format='(D10.1)'), /rem))+' km/s'
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), set_value=science.nloop gt 0 ? strcompress(string(science.nloop, format='(I4)'), /rem) : unknown

end


pro sps_fit::getscience, files=files
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common mask_in, mask_in, copynum
    common npixcom, npix
    common sps_spec, sps, spsz, spsage
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest

;choi data ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
massbinval = [10.5,10.8,11.0,11.3,10.8,11.1,11.3,10.9,11.0,11.3]
fehchoi    = [-0.11,-0.05,-0.02,-0.03,-0.07,-0.04,-0.05,-0.15,-0.02,-0.05]
fehchoierr = [0.03,0.01,0.01,0.02,0.02,0.01,0.02,0.07,0.03,0.03]
agechoi    = [3.27,3.47,4.55,5.61,2.99,3.28,4.00,2.67,2.49,3.06]
agechoierr = [0.21,0.12,0.13,0.15,0.15,0.13,0.18,0.20,0.12,0.11]
CFe  = [0.19,0.16,0.18,0.19,0.18,0.18,0.14,0.16,0.24,0.26]
NFe  = [0.26,0.25,0.24,0.26,0.26,0.34,0.25,0.18,0.58,0.73]
NFeerr = [0.10,0.05,0.03,0.05,0.08,0.05,0.06,0.30,0.09,0.09]
MgFe = [0.18,0.22,0.21,0.29,0.23,0.23,0.30,0.05,0.09,0.19]
MgFeerr = [0.04,0.02,0.01,0.03,0.04,0.02,0.03,0.13,0.05,0.04]
CaFe = [0.04,0.03,0.04,0.04,0.05,-0.03,0.08,0.06,-0.03,0.00]
choi = {logmstar:massbinval,feh:fehchoi,feherr:fehchoierr,age:agechoi,ageerr:agechoierr,CFe:CFe,NFe:NFe,NFeerr:NFeerr,MgFe:MgFe,MgFeerr:MgFeerr}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    npix = 1973
    observatory, 'apo', obs
    sciencefits = self.directory+'sps_fit'+copynum+'.fits.gz'
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        c = n_elements(files)
        zranges = strarr(c)
        massbins = strarr(c)
        for i=0,c-1 do begin
            basefile = file_basename(files[i])
            extensions = strsplit(basefile, '_', /extract)
            zranges[i] = extensions[1]
            massbins[i] = extensions[2]
        endfor
        nspec = n_elements(zranges)
        self.nspec = nspec
        speclist = zranges+' '+massbins
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1
        for i=0,nspec-1 do begin
            science = {science}
            science.skyfit = -1
            science.skylinemask = -1
            science.mask = mask_in
            science.zrange = zranges[i]
            science.massbin = massbins[i]
            case science.zrange of 
               'z0.30.4': science.zspec = 0.35
               'z0.40.55': science.zspec = 0.475
               'z0.550.7': science.zspec = 0.625
               else: print,'no matched redshift'
            endcase
            science.logmstar = choi.logmstar[i]
            science.fehchoi = choi.feh[i]
            science.fehchoierr = choi.feherr[i]
            science.agechoi = choi.age[i]
            science.agechoierr = choi.ageerr[i]
            science.MgFechoi = choi.mgfe[i]
            science.MgFechoierr = choi.mgfeerr[i]
            readcol,files[i],lambda,flux,fluxerr,contmask
            data1={lambda:lambda,flux:flux,ivar:1./fluxerr^2,contmask:contmask}
            science.objname = speclist[i]


            nlamb = n_elements(lambda)
            if nlamb gt npix then data1=data1[0:npix-1]
            if nlamb lt npix then begin
               sdss_nan = data1[0]
               for ntag=0,n_elements(tag_names(data1))-1 do sdss_nan.(ntag) = 1./0.
               nanarr   = replicate(sdss_nan,npix-nlamb)
               nanarr.ivar = 0
               data1 = [data1,nanarr]
            endif
            lambda    = data1.lambda
            spec      = data1.flux
            ivar      = data1.ivar
            science.fitmask = data1.contmask
            science.contmask = data1.contmask
            
            science.spec1dfile = files[i]
            science.good = 1

            w = where(science.age le 0, c)
            if c gt 0 then begin
                science[w].age = -999d
                science[w].ageerr = -999d
            endif
            science.feh = -999d
            science.feherr = -999d
            science.vdisp = -999d
            science.vdisperr = -999d
            science.alphafe = 0.0

            self.i = i

            w = where(ivar gt 0 and finite(ivar), civar)
            if civar gt 10 then if w[civar-1] ne n_elements(ivar)-1 then ivar[w[civar-11:civar-1]] = 0
      
            n = n_elements(lambda)

            t = (-1*ts_diff(lambda, 1))[0:n-2]
            wt = where(t le 0, ct)
            if ct gt 0 then begin
                message, 'Wavelength array for '+strtrim(objnames[i], 2)+' is not monotonic.  Omitting.', /info
                wgood[i] = 0
                continue
            endif

            science.lambda = lambda*(1.+science.zspec)
            science.spec = spec
            science.ivar = ivar

            self->specres, science
            self->continuum, science
            if array_equal(science.continuum, replicate(-999, n_elements(science.continuum))) then begin
                message, 'Failed at finding continuum for '+strtrim(objnames[i], 2)+'.  Omitting.', /info
                wgood[i] = 0
                continue
            endif
            self->mask, science
            science.spscont = 1.0

            ;get choi's best fit spectrum
             dlam = science.dlam
             datalam = science.lambda
             wonfit = where(datalam)
             ;;get the model
             rest = 0
             normalize = 0
             spschoi = get_sps_choi_alpha_obs(datalam,[science.fehchoi,science.agechoi,250.,science.zspec,science.mgfechoi])
             science.spsspecfullchoi = spschoi
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            self->statusbox, science=science
            scienceall[i] = science
        endfor
        self.i = 0
        wgood = where(wgood eq 1, cgood)
        scienceall = scienceall[wgood]
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        self.nspec = cgood
        self->writescience
        objlist = scienceall.objname
        speclist = mask_in+' '+strtrim(string(objlist[wgood]), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
        self->writescience
     endif else begin
        scienceall = mrdfits(sciencefits, 1, /silent)

        self.nspec = n_elements(scienceall)
        speclist = mask_in+' '+strtrim(string(scienceall.objname), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
    endelse
  end

pro sps_fit::getscience_ngc, files=files
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common mask_in, mask_in, copynum
    common npixcom, npix

;Walcher09/Conroy2014 data ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nspec = n_Elements(files)
massbinval = 8.+fltarr(nspec)
fehchoi    = -0.1+fltarr(nspec)
fehchoierr = 0.05+fltarr(nspec)
agechoi    = 11.+fltarr(nspec)
agechoierr = 2.+fltarr(nspec)
CFe  = fltarr(nspec)
NFe  = fltarr(nspec)
NFeerr =  fltarr(nspec)
MgFe =  0.1+fltarr(nspec)
MgFeerr = 0.1+fltarr(nspec)
CaFe = fltarr(nspec)
choi = {logmstar:massbinval,feh:fehchoi,feherr:fehchoierr,age:agechoi,ageerr:agechoierr,CFe:CFe,NFe:NFe,NFeerr:NFeerr,MgFe:MgFe,MgFeerr:MgFeerr}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    npix = 3071
    observatory, 'apo', obs
    sciencefits = self.directory+'sps_fit'+copynum+'.fits.gz'
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        exposure = strarr(nspec)
        aperture = strarr(nspec)
        auxfile = strarr(nspec)
        for i=0,nspec-1 do begin
            basefile = file_basename(files[i])
            extensions = strsplit(basefile, '_', /extract)
            exposure[i] = extensions[1]
            aperture[i] = strmid(extensions[2],0,1)
            auxfile[i] = '/scr2/nichal/workspace2/Choi14_spec/NGC6528_aux/'+strmid(basefile,0,strpos(basefile,'.fits'))+'.aux.fits'
            if file_test(auxfile[i]) eq 0 then message,'no auxillary file found'
        endfor
        self.nspec = nspec
        speclist = exposure+aperture
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1
        for i=0,nspec-1 do begin
            science = {science}
            science.skyfit = -1
            science.skylinemask = -1
            science.mask = mask_in
            science.zrange = 'local GC'
            science.massbin = -99.
            science.zspec = 0.1               
            science.logmstar = choi.logmstar[i]
            science.fehchoi = choi.feh[i]
            science.fehchoierr = choi.feherr[i]
            science.agechoi = choi.age[i]
            science.agechoierr = choi.ageerr[i]
            science.MgFechoi = choi.mgfe[i]
            science.MgFechoierr = choi.mgfeerr[i]
            ;read data file
            flux = mrdfits(files[i],0)
            aux = mrdfits(auxfile[i],0)
            lambda = findgen(npix)+3360.
            ivar = (aux[*,3]/flux)^2  ;aux[0,*] is S/N
            contmask = bytarr(npix)+1
            contmask(where(lambda gt 4543 and lambda lt 4552)) = 0
            contmask(where(lambda gt 5045 and lambda lt 5051)) = 0
            contmask(where(lambda gt 5575 and lambda lt 5582)) = 0
            contmask(where(lambda gt 6055 and lambda lt 6075)) = 0

            data1={lambda:lambda,flux:flux,ivar:ivar,contmask:contmask}
            science.objname = speclist[i]


            nlamb = n_elements(lambda)
            if nlamb gt npix then data1=data1[0:npix-1]
            if nlamb lt npix then begin
               sdss_nan = data1[0]
               for ntag=0,n_elements(tag_names(data1))-1 do sdss_nan.(ntag) = 1./0.
               nanarr   = replicate(sdss_nan,npix-nlamb)
               nanarr.ivar = 0
               data1 = [data1,nanarr]
            endif
            lambda    = data1.lambda
            spec      = data1.flux
            ivar      = data1.ivar
            science.fitmask = data1.contmask
            science.contmask = data1.contmask
            
            science.spec1dfile = files[i]
            science.good = 1

            w = where(science.age le 0, c)
            if c gt 0 then begin
                science[w].age = -999d
                science[w].ageerr = -999d
            endif
            science.feh = -999d
            science.feherr = -999d
            science.vdisp = -999d
            science.vdisperr = -999d
            science.alphafe = 0.0

            self.i = i

            w = where(ivar gt 0 and finite(ivar), civar)
            if civar gt 10 then if w[civar-1] ne n_elements(ivar)-1 then ivar[w[civar-11:civar-1]] = 0
      
            n = n_elements(lambda)

            t = (-1*ts_diff(lambda, 1))[0:n-2]
            wt = where(t le 0, ct)
            if ct gt 0 then begin
                message, 'Wavelength array for '+strtrim(objnames[i], 2)+' is not monotonic.  Omitting.', /info
                wgood[i] = 0
                continue
            endif

            science.lambda = lambda*(1.+science.zspec)
            science.spec = spec
            science.ivar = ivar

            self->specres, science
            self->continuum, science
            if array_equal(science.continuum, replicate(-999, n_elements(science.continuum))) then begin
                message, 'Failed at finding continuum for '+strtrim(objnames[i], 2)+'.  Omitting.', /info
                wgood[i] = 0
                continue
            endif
            self->mask, science
            science.spscont = 1.0
            self->indices, science, /noredraw

            self->statusbox, science=science
            scienceall[i] = science
        endfor
        self.i = 0
        wgood = where(wgood eq 1, cgood)
        scienceall = scienceall[wgood]
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        self.nspec = cgood
        self->writescience
        objlist = scienceall.objname
        speclist = mask_in+' '+strtrim(string(objlist[wgood]), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
        self->writescience
     endif else begin
        scienceall = mrdfits(sciencefits, 1, /silent)

        self.nspec = n_elements(scienceall)
        speclist = mask_in+' '+strtrim(string(scienceall.objname), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
    endelse
end


pro sps_fit::writescience
    common mask_in, mask_in, copynum
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing to database ...'
    scienceall = *self.science
    sciencefits = self.directory+'sps_fit'+copynum+'.fits'
    mwrfits, scienceall, sciencefits, /create, /silent
    spawn, 'gzip -f '+sciencefits
    widget_control, widget_info(self.base,find_by_uname='status'), set_value='Ready.'
    print, 'saved file'
end


pro sps_fit::initialize_directory, directory=directory
    common mask_in, mask_in, copynum
    newdirectory = '/scr2/nichal/workspace4/sps_fit/data/choi/'+mask_in+'/'
    if ~file_test(newdirectory) then file_mkdir, newdirectory

    if strmid(directory, 0, 1, /reverse_offset) ne '/' then directory = directory + '/'
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

    countfiles:
    if mask_in eq 'NGC6528' then files = file_search(directory,'*.fits',count=c) else files = file_search(directory, '*.dat', count=c)

    sciencefits = newdirectory+'sps_fit'+copynum+'.fits.gz'
    if c eq 0 then begin
        files = file_search(directory+'*/*.dat', count=c)
    endif
    if c eq 0 and ~file_test(sciencefits) then begin
        message, 'Unknown mask.'
    endif

    self.directory = newdirectory
    if mask_in eq 'NGC6528' then self->getscience_ngc,files=files else $
       self->getscience, files=files
    self.i = 0
    science = (*self.science)[self.i]
    self->statusbox, science=science
    self->default_range
    self.lambdalim = [3300, 7000]
    self->newspec, increment=0    
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
    
end


; =============== INIT ================
function sps_fit::INIT, directory=directory, lowsn=lowsn
    common mask_in, mask_in, copynum
    common sps_spec, sps, spsz, spsage
    if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)

    base = widget_base(/row, title='sps_fit_choi'+copynum, uvalue=self, mbar=menu, tab_mode=0, units=1)
    ;------Top menu--------
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wdefaultrange = widget_button(tools_menu, value='Default Spectrum Settings', uname='default_range', uvalue='default_range')
    wdefault_cont = widget_button(tools_menu, value='Default Continuum Regions', uname='default_cont', uvalue='default_cont')
    wdefault_maskall = widget_button(tools_menu, value='Default Pixel Mask All', uname='default_maskall', uvalue='default_maskall')
    wdefault_goodspec = widget_button(tools_menu, value='Default Good Spectrum', uname='default_goodspec', uvalue='default_goodspec')
    wreprepare_all = widget_button(tools_menu, value='Reprepare All', uname='reprepare_all', uvalue='reprepare_all')
    wfit_all = widget_button(tools_menu, value='Fit All', uname='fit_all', uvalue='fit_all')
    wfitalpha_all = widget_button(tools_menu, value='Fit Alpha All', uname='fit_alpha_all', uvalue='fit_alpha_all')
    wfitomg_all = widget_button(tools_menu, value='Fit OMg All', uname='fit_omg_all', uvalue='fit_omg_all')
    wfitglorious_all = widget_button(tools_menu, value='Fit Glorious All', uname='fit_glorious_all', uvalue='fit_glorious_all')
    wcal_uncertainties_all = widget_button(tools_menu, value='Cal uncertainties All',uname='cal_uncertainties_all',uvalue='cal_uncertainties_all')
    wcal_uncertainties_alpha_all = widget_button(tools_menu, value='Cal uncertainties alpha All',uname='cal_uncertainties_alpha_all',uvalue='cal_uncertainties_alpha_all')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base

    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit', 'continuum division', 'sky line fit', 'rest frame'], /column, /exclusive, set_value=3, uname='mode', uvalue='mode', /no_release)
    wstep = widget_base(wleft, /row, /align_center)
    wbackward = widget_button(wstep, value='<---', uvalue='backward', uname='backward', tab_mode=1, xsize=75)
    wforward = widget_button(wstep, value='--->', uvalue='forward', uname='forward', tab_mode=1, xsize=75)
    wprepbase = widget_base(wleft, /row, /align_center)
    wfitalphaband = widget_button(wprepbase, value='Fit alpha band', uvalue='fitalphaband', uname='fitalphaband', tab_mode=1)
    wfit = widget_button(wprepbase, value='Fit', uvalue='fit', uname='fit', tab_mode=1, xsize=75)
    wfitalpha = widget_button(wprepbase, value='Fit alpha', uvalue='fitalpha', uname='fitalpha', tab_mode=1, xsize=85)
    wfitglorious = widget_button(wprepbase, value='Fit glorious', uvalue='fitglorious', uname='fitglorious', tab_mode=1, xsize=85)
    windicesbase = widget_base(wleft, /row, /align_center)
    wfitomg = widget_button(windicesbase,value='Fit O Mg',uvalue='fitomg',uname='fitomg',tab_mode=1,xsize=100)
    wdefault_mask = widget_button(windicesbase,value='Default Mask',uvalue='default_mask',uname='default_mask',tab_mode=1,xsize=100)
    wuncertbase = widget_base(wleft,/row,/align_center)
    wuncertainties = widget_button(wuncertbase,value='cal_uncertainties',uvalue='cal_uncertainties',uname='cal_uncertainties')
    wuncertainties = widget_button(wuncertbase,value='cal_uncertainties_alpha',uvalue='cal_uncertainties_alpha',uname='cal_uncertainties_alpha')
    wgoodbase = widget_base(wleft, /column, /align_center)
    wgood = cw_bgroup(wgoodbase, ['good spectrum','good fit'], /nonexclusive, set_value=[0,0], uname='good', uvalue='good')

    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=250)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=195)
    wmstarbase = widget_base(wcurobj, /align_center, /row)
    wmstarlabel = widget_label(wmstarbase, value='log M* = ', /align_right, uname='mstarlabel', xsize=95)
    wcurmstar = widget_label(wmstarbase, value='     ', /align_left, uname='curmstar', uvalue='curmstar', xsize=150)
    wagebase = widget_base(wcurobj, /align_center, /row)
    wagelabel = widget_label(wagebase, value='age = ', /align_right, uname='agelabel', xsize=95)
    wcurage = widget_label(wagebase, value='     ', /align_left, uname='curage', uvalue='curage', xsize=150)
    wageuncertbase = widget_base(wcurobj, /align_center, /row)
    wageuncertlabel = widget_label(wageuncertbase,value='    ',/align_right,uname='ageuncertlabel',xsize=95)
    wageuncert = widget_label(wageuncertbase,value='      ',/align_left, uname='curageuncert', uvalue='curageuncert', xsize=150)
    wagechoibase = widget_base(wcurobj, /align_center, /row)
    wagechoilabel = widget_label(wagechoibase, value='agechoi = ', /align_right, uname='agechoilabel', xsize=95)
    wcuragechoi = widget_label(wagechoibase, value='     ', /align_left, uname='curagechoi', uvalue='curagechoi', xsize=150)
    wfehbase = widget_base(wcurobj, /align_center, /row)
    wfehlabel = widget_label(wfehbase, value='[Fe/H] = ', /align_right, uname='fehlabel', xsize=95)
    wcurfeh = widget_label(wfehbase, value='     ', /align_left, uname='curfeh', uvalue='curfeh', xsize=150)
    wfehuncertbase = widget_base(wcurobj, /align_center, /row)
    wfehuncertlabel = widget_label(wfehuncertbase,value='    ',/align_right,uname='fehuncertlabel',xsize=95)
    wfehuncert = widget_label(wfehuncertbase,value='      ',/align_left, uname='curfehuncert', uvalue='curfehuncert', xsize=150)
    wfehchoibase = widget_base(wcurobj, /align_center, /row)
    wfehchoilabel = widget_label(wfehchoibase, value='[Fe/H] choi = ', /align_right, uname='fehchoilabel', xsize=95)
    wcurfehchoi = widget_label(wfehchoibase, value='     ', /align_left, uname='curfehchoi', uvalue='curfehchoi', xsize=150)
    wmgbase = widget_base(wcurobj, /align_center, /row)
    wmglabel = widget_label(wMgbase, value='[Mg/Fe] = ', /align_right, uname='mglabel', xsize=95)
    wcurMg = widget_label(wMgbase, value='     ', /align_left, uname='curmg', uvalue='curmg', xsize=150)
    wmguncertbase = widget_base(wcurobj, /align_center, /row)
    wmguncertlabel = widget_label(wmguncertbase,value='    ',/align_right,uname='mguncertlabel',xsize=95)
    wmguncert = widget_label(wmguncertbase,value='      ',/align_left, uname='curmguncert', uvalue='curmguncert', xsize=150)
    wmgfechoibase = widget_base(wcurobj, /align_center, /row)
    wmgfechoilabel = widget_label(wmgfechoibase, value='[Mg/Fe] choi = ', /align_right, uname='mgfechoilabel', xsize=95)
    wcurmgfechoi = widget_label(wmgfechoibase, value='     ', /align_left, uname='curmgfechoi', uvalue='curmgfechoi', xsize=150)
    wvdispbase = widget_base(wcurobj, /align_center, /row)
    wvdisplabel = widget_label(wvdispbase, value='sigma_v = ', /align_right, uname='vdisplabel', xsize=95)
    wcurvdisp = widget_label(wvdispbase, value='     ', /align_left, uname='curvdisp', uvalue='curvdisp', xsize=150)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='SNR = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)
    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='z = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)
    wzfitbase = widget_base(wcurobj, /align_center, /row)
    wzfitlabel = widget_label(wzfitbase, value='zfit = ', /align_right, uname='zfitlabel', xsize=95)
    wcurzfit = widget_label(wzfitbase, value='     ', /align_left, uname='curzfit', uvalue='curzfit', xsize=150)
    wchisqbase = widget_base(wcurobj, /align_center, /row)
    wchisqlabel = widget_label(wchisqbase, value='chisq = ', /align_right, uname='chisqlabel', xsize=95)
    wcurchisq = widget_label(wchisqbase, value='     ', /align_left, uname='curchisq', uvalue='curchisq', xsize=150)
    wnloopbase = widget_base(wcurobj, /align_center, /row)
    wnlooplabel = widget_label(wnloopbase, value='nloop = ', /align_right, uname='nlooplabel', xsize=95)
    wcurnloop = widget_label(wnloopbase, value='     ', /align_left, uname='curnloop', uvalue='curnloop', xsize=150)
    wnpolybase = widget_base(wcurobj, /align_center, /row)
    wnpolylabel = widget_label(wnpolybase, value='npoly = ', /align_right, uname='npolylabel', xsize=95)
    wcurnpoly = widget_label(wnpolybase, value='     ', /align_left, uname='curnpoly', uvalue='curnpoly', xsize=150)
    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wback = widget_button(wfile, value='Back', uvalue='back', uname='back', tab_mode=1)
    wfilelist = widget_combobox(wfile, uname='filelist', value='                 ', tab_mode=1, /dynamic_resize)
    wnext = widget_button(wfile, value='Next', uvalue='next', uname='next', tab_mode=1)
    wstatus = widget_text(wfile, xsize=108, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)

    wspec = widget_base(wright, /frame, /column)
    wspecplot = widget_draw(wspec, xsize=1600, ysize=400, uname='spec', /button_events, keyboard_events=1)
    w2dplot = widget_draw(wspec, xsize=1600, ysize=300, uname='2d', /tracking_events, /motion_events)

    wspeccontrol = widget_base(wright, /row, /align_center, tab_mode=1)
    wkeepoldfitcontrol = widget_base(wspeccontrol,/row,/frame)
    wkeepoldfit = cw_bgroup(wkeepoldfitcontrol, ['compare old chisq','dont compare(replace)'],/exclusive, set_value=0, uname='keepoldfit', uvalue='keepoldfit',row=1)
    wnpolycontrol = widget_base(wspeccontrol, /frame, /row)
    wnpolylabel = widget_label(wnpolycontrol, value='npoly=', /align_center, uname='npolylabel')
    wnpoly = widget_text(wnpolycontrol, xsize=8, /editable, uname='npoly', uvalue='npoly')
    wycontrol = widget_base(wspeccontrol, /frame, /row)
    wylow = widget_text(wycontrol, xsize=8, /editable, uname='ylow', uvalue='ylow')
    wylabel = widget_label(wycontrol, value=' < y < ', /align_center, uname='ylabel')
    wyhigh = widget_text(wycontrol, xsize=8, /editable, uname='yhigh', uvalue='yhigh')
    wlambdacontrol = widget_base(wspeccontrol, /frame, /row)
    wblue = widget_button(wlambdacontrol, value='<-', uname='blue', uvalue='blue', /align_center)
    wlambdalow = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdalow', uvalue='lambdalow')
    wlambdalabel = widget_label(wlambdacontrol, value=' < l < ', /align_center, uname='lambdalabel', font='-urw-standard symbols l-medium-r-normal--0-0-0-0-p-0-adobe-symbol')
    wlambdahigh = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdahigh', uvalue='lambdahigh')
    wred = widget_button(wlambdacontrol, value='->', uname='red', uvalue='red', /align_center)
    wnloopcontrol = widget_base(wspeccontrol,/frame,/row)
    wnloop = widget_text(wnloopcontrol,xsize=8,/editable,uname='maxnloop',uvalue='maxnloop')

    widget_control, base, /realize
    self.base = base
    xmanager, 'sps_fit', self.base, /no_block, cleanup='sps_fit_cleanup'

    readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
    wbands = [1,2,3,4,5]
    tellstart = tellstart[wbands]
    tellend = tellend[wbands]
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick, self.indstart, self.indend, self.indname
    self.linewaves = ptr_new([2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85])
    self.linenames = ptr_new(['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', '[NeIII]', 'HeI', 'H8', 'CaH', '[NeIII]', 'CaK', 'He', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', 'Mgb', 'Mgb', 'Mgb', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]'])
    self.linecolors = ptr_new(['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue'])
    self.tellstart = ptr_new(tellstart)
    self.tellend = ptr_new(tellend)
    self.tellthick = ptr_new([5, 2, 5, 2, 2])

    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    self.linestart = ptr_new(linestart)
    self.lineend = ptr_new(lineend)
    self.linetype = ptr_new(linetype)

    readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indbandstart,indbandend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,format='I,D,D,D,D,D,D,I,A',comment='#',/silent
    self.indstart = ptr_new(indbandstart)
    self.indend   = ptr_new(indbandend)
    self.indname  = ptr_new(indname)
   degendir = '/scr2/nichal/workspace4/sps_fit/sspdegen/sspdegen_choi/'
   degenfile = file_search(degendir+'age*_sspdegen_gridspec_choi.fits',count=cdegenfile)
   for i=0,cdegenfile-1 do begin
      degen_sps = mrdfits(degenfile(i),1,/silent)
      if i eq 0 then begin
           fits_info,degenfile(i),/silent,n_ext=nalpha
           nfeh = n_elements(degen_sps)
           nage = cdegenfile
           degen_fehgrid = fltarr(nfeh,nage,nalpha)
           degen_agegrid = fltarr(nfeh,nage,nalpha)
           degen_alphagrid = fltarr(nfeh,nage,nalpha)
           degen_file = strarr(nage)
           for j=0,nalpha-1 do begin
              degen_sps = mrdfits(degenfile(i),j+1,/silent)
              degen_alphagrid[*,*,j] = rebin(degen_sps.alpha,nfeh,nage)
           endfor
      endif
      filebase = file_basename(degenfile(i),'.fits')
      agei = fix(strmid(filebase,3,2))

      degen_file[agei] = degenfile(i)
      degen_fehgrid[*,agei,*] = rebin(degen_sps.feh,nfeh,nalpha)
      degen_agegrid[*,agei,*] = rebin(degen_sps.age,nfeh,nalpha)
    endfor
    self.degen_file = ptr_new(degen_file)
    self.degen_fehgrid = ptr_new(degen_fehgrid)
    self.degen_agegrid = ptr_new(degen_agegrid)
    self.degen_alphagrid = ptr_new(degen_alphagrid)
    common random, seed
    seed = systime(1)
    
    self.i = -1
    self.lowsn = keyword_set(lowsn) ? 1 : 0
    self->initialize_directory, directory=directory
    return, 1
end

pro sps_fit__define 
    state = {sps_fit, $
             base:0L, $
             directory:'', $
             science:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
             divylim:[-100d, 100d], $
             skyylim:[-100d, 100d], $
             linewaves:ptr_new(), $
             linenames:ptr_new(), $
             linecolors:ptr_new(), $
             tellstart:ptr_new(), $
             tellend:ptr_new(), $
             tellthick:ptr_new(), $
             linestart:ptr_new(), $
             lineend:ptr_new(), $
             linetype:ptr_new(), $
             indstart:ptr_new(), $
             indend:ptr_new(), $
             indname:ptr_new(), $
             degen_file:ptr_new(), $
             degen_fehgrid:ptr_new(),$
             degen_agegrid:ptr_new(),$
             degen_alphagrid:ptr_new(), $
             nspec:0L, $
             i:0L, $
             keystate:0, $
             lambda1:0d, $
             lowsn:0b}
end


pro science__define
    common npixcom, npix
    nsky = 200
    science = {science, $
               objname:'', $
               mask:'', $
               zrange:'', $
               massbin:'', $
               lambda:dblarr(npix), $
               spec:dblarr(npix), $
               ivar:dblarr(npix), $
               skyspec:dblarr(npix), $
               contmask:bytarr(npix), $
               continuum:dblarr(npix), $
               spscont:dblarr(npix), $
               contdiv:dblarr(npix), $
               contdivivar:dblarr(npix), $
               spsspec:dblarr(npix), $
               spsspecfull:dblarr(npix), $
               spsspecfullmg:dblarr(npix), $
               spsspecfullchoi:dblarr(npix), $
               spscontfull:dblarr(npix), $
               fitmask:bytarr(npix), $
               dlam:dblarr(npix), $
               skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
               skylinemask:lonarr(nsky), $
               goodsky:1B, $
               nloop:-999d, $
               zspec:-999d, $
               zfit:-999d, $
               zsource:0, $
               age:-999d, $
               ageerr:-999d, $
               ageupper:-999d, $
               agelower:-999d, $
               feh:-999d, $
               feherr:-999d, $
               fehupper:-999d, $
               fehlower:-999d,$
               agechoi:-999d, $
               agechoierr:-999d, $
               fehchoi:-999d, $
               fehchoierr:-999d, $
               MgFechoi:-999d, $
               MgFechoierr:-999d, $
               alphafe:-999d, $
               alphafeerr:-999d, $
               alphafeupper:-999d, $
               alphafelower:-999d, $
               mgfe:-999d, $
               mgfeupper:-999d, $
               mgfelower:-999d, $
               chisqmg:-999d, $
               logmstar:-999d, $
               vdisp:-999d, $
               vdisperr:-999d, $
               chisq:-999d, $
               sn:-999d, $
               snfit:-999d, $
               spec1dfile:'', $
               good:1B,$
               goodfit:1B,$
               npoly:0,$
               pfit:fltarr(13), $
               pfiterr:fltarr(13)}
end


pro sps_fit_choi,copyi=copyi
    common mask_in, mask_in, copynum
    common get_sps_alpha, element
    mask_in = 'ages'
    if ~keyword_set(copyi) then copyi=1
    copynum = strtrim(string(copyi,format='(I02)'),2)
    directory = '/scr2/nichal/workspace2/Choi14_spec/'+mask_in
    if ~file_test(directory) then message, 'Mask not found.'
    case copyi of
       1: element = ['Mg','O','Si','Ca','Ti'] 
       6: element = ['Mg','O','Si','Ca','Ti'] 
       3: element = ['Mg']
       5: element = ['Mg']
       7: element = ['Mg']
       8: element = ['Mg']
       9: element = ['Mg','O']
       10: element = ['Mg','N']
       11: element = ['Mg','N','Fe','O','C','N','Si','Ca','Ti']
       12: element = ['Mg','Fe']
       else: element=['nope']
    endcase
    n = obj_new('sps_fit', directory=directory)
end
