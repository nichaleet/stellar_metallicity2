; =================
pro sps_iterproc, funcname, p, iter, fnorm, functargs=functargs, parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    if iter gt 1 then print, contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
end


pro sps_fit::indices, science, noredraw=noredraw, nostatusbar=nostatusbar
    nmc = 1000
    cah = 3933.663
    cahindex = dblarr(nmc+1) - 1d9
    gbandindex = dblarr(nmc+1) - 1d9

    restlambda = science.lambda / (1d + science.zspec)
    flux = science.contdiv
    ivar = science.contdivivar
    w = where(finite(flux) and finite(ivar) and ivar gt 0, n)
    if n lt 10 then return
    restlambda = restlambda[w]
    flux = flux[w]
    ivar = ivar[w]
    wcah = where(restlambda gt cah-9 and restlambda lt cah+9, ccah)
    wcahside = where((restlambda gt 3900 and restlambda lt 3916) or (restlambda gt 4015 and restlambda lt 4029), ccahside)
    wgband = where(restlambda gt 4285 and restlambda lt 4318, cgband)
    wgbandside = where(restlambda gt 4228 and restlambda lt 4273, cgbandside)
    for i=0,nmc do begin
        if i eq 0 then fluxi = flux else fluxi = flux + (ivar)^(-0.5)*randomn(seed, n)
        if ccah gt 10 and ccahside gt 10 then cahindex[i] = weightedmean(1.0 - fluxi[wcah], (ivar[wcah])^(-0.5)) - weightedmean(1.0 - fluxi[wcahside], (ivar[wcahside])^(-0.5))
        if cgband gt 10 and cgbandside gt 10 then gbandindex[i] = weightedmean(1.0 - fluxi[wgband], (ivar[wgband])^(-0.5)) - weightedmean(1.0 - fluxi[wgbandside], (ivar[wgbandside])^(-0.5))
    endfor
    science.cah = cahindex[0]
    w = where(cahindex[1:nmc] gt -1d8)
    science.caherr = stddev(cahindex[w+1])
    science.gband = gbandindex[0]
    w = where(gbandindex[1:nmc] gt -1d8)
    science.gbanderr = stddev(gbandindex[w+1])
end


pro sps_fit::fit, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize, rest
    common toprint, agediff, zdiff

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting ...'
    savedata = 0 ;1

    znow = science.zspec
    reallambda = science.lambda
    nlambda = n_elements(reallambda)
    
    ;check if the blue chip was failed
    neg = where(reallambda lt 0., cneg)
    if cneg gt 0 then reallambda(neg) = reallambda(neg)+10000.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if min(science.dlam) gt 0.2 and max(science.dlam) lt 10.0 then begin
        dlam_all = science.dlam
    endif else begin
        specresfile = self.directory+'specres_poly.sav'
        if file_test(specresfile) then begin
            restore, specresfile
            dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
    endelse

    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
   
    pi[0].limits = [-0.5,0.1]
    pi[1].limits = [min(spsage),alog10(galage(znow,1000))<max(spsage)]
   ; pi[1].limits = [0,1.]
    pi[3].limits = [-0.05,0.05]+znow
    ;set the prior to velocity dispersion according to Faber-Jackson relation (Dutton2011)
    if science.logmstar gt 5. then begin
	logvdisp = 2.23+0.37*(science.logmstar-10.9)-0.19*alog10(0.5+0.5*(10.^science.logmstar/10.^10.9))
	pi[2].limits = [10.^(logvdisp-0.4),10.^(logvdisp+0.4)]
	print, 'velocity dispersion prior = ',pi[2].limits,' km/s'	
    endif else pi[2].limits = [40.,400.]; [30.,200.]
   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,4)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    ;pi[2].value = 106.
    ;pi[0].value = -0.2
    firstguess = pi.value    
    pi[0].limits = minmax(spsz) ;fix the limit back
    ;pi[0].limits = [-0.5,0.1]
    pi[1].limits = [min(spsage),alog10(galage(znow,1000))<max(spsage)]
    pi.step = double([0.1, 0.5, 25.0,0.0002])
    pi.parname = ['    Z', '  age', 'vdisp','redshift']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
    print, 'prior range:',pi.limits

    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)

    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d,-999d]
        perror = [-999d, -999d, -999d, -999d]
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
    maxnloop = 100
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF   ZDIFF  AGEDIFF'
    print, '--- ------- ----- ------- ---------  -------- ---- -----  ------'

;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.]
    besterror = [99.,99.,99.,99.]
    chisqarr = fltarr(maxnloop)
    valuearr = fltarr(4,maxnloop)
    errorarr = fltarr(4,maxnloop)
    
    ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop le maxnloop do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.telldivivar*(median(science.telldiv))^2
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc')

        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        vdispdiff = (pi[2].value-pars[2])/pi[2].value
        redshfdiff = (pi[3].value-pars[3])/pi[2].value
        ;print,'z,age,vdisp diff:',zdiff,agediff,vdispdiff,nloop
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_obs(reallambda, pars)                                
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

       ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont  
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.telldiv/median(science.telldiv)
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')
       
        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=150, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d, -999d]
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

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof

        ;;keeping the values
        if savedata eq 1 and nloop mod 10 eq 0 then begin
           if nloop eq 0 then $
              str={won:won,lambda:science.lambda,contdiv:science.contdiv,contdivivar:science.contdivivar,fakecont:science.spscont,spsbestfit:spsbestfitarr[*,0],chisq:curchisq,param:pars,perror:perror,loop:nloop} else $
                 str={won:won,lambda:science.lambda,contdiv:science.contdiv,contdivivar:science.contdivivar,fakecont:science.spscont,spsbestfit:spsbestfit,chisq:curchisq,param:pars,perror:perror,loop:nloop}
           if nloop eq 0 then loopstr = [str] else loopstr=[loopstr,str]
           save,loopstr,firstguess,filename=self.directory+'/obj'+strtrim(string(self.i+1, format='(I3)'), 2)+'sps_fit_loopinfo.sav'
        endif
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

    if savedata eq 1 then begin
       str = {chisq:chisqarr,param:valuearr,perror:errorarr,paraname:['Z','age','sigmav','redshift'],objname:strtrim(science.objname, 2),dof:dof}
       save,str,filename=self.directory+'/sps_fit_data_'+strtrim(string(self.i+1, format='(I3)'), 2)+'.sav'
    endif

    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'                            
    ;;check if the last chisq is the best chisq
    ;if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
    ;   print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
    ;   print,'The values used are:'
    ;   print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3)'
    ;   science.goodfit = 1.
    ;   spsbestfitarr = bestspsbestfitarr
    ;   pi.value = bestvalue
    ;   perror   = besterror
    ;   science.spscont = bestcont
    ;endif

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
    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       wcont = wcont[where(finite(science.telldiv[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       wcont = wcont[where(finite(science.telldiv[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fit_all
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        self->default_range
        science = scienceall[self.i]
        if science.good eq 0 then continue
	;if science.good eq 0 and science.goodfit eq 0 then continue
        self->fit, science
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


; =================
pro sps_fit_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)

    case (name) of
        'BUTT': obj->handle_button, ev
        'TEXT': obj->handle_text, ev
        'TRAC': obj->handle_tracking, ev
        'COMB': obj->handle_combobox, ev
        'DRAW': begin
            if ev.type eq 0 then obj->handle_draw_click, ev
            if ev.type eq 2 then obj->handle_motion, ev
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
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    
    case (uvalue) of
        'back': self->newspec, increment=-1
        'next': self->newspec, increment=1
        'default_range': begin
            self->default_range
            self->redraw
        end
        'indices': begin
            scienceall = *self.science
            science = scienceall[self.i]
            self->indices, science
            scienceall[self.i] = science
            ptr_free, self.science
            self.science = ptr_new(scienceall)
            self->statusbox
            ;self->redraw
         end
        'fit': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->fit, science, /noredraw
           scienceall[self.i] = science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
        end
        'fit[Mg/Fe]': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->maskmg,science
           self->fitmg, science, /noredraw
           scienceall[self.i] = science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
        end
        'fit_all': self->fit_all
        'default_cont': self->default_cont
        'default_mask': self->default_mask
        'default_maskall': self->default_maskall
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
    nmodes = 3
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 2 and increment gt 0: mode = 0
        mode eq 0 and increment lt 0: mode = 2
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
    scienceall[self.i].goodsky = good[0]
    scienceall[self.i].good = good[1]
    scienceall[self.i].goodfit = good[2]
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
    ;good = 0
    ;while(~good) do begin
        newi += increment
        if newi lt 0 then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the first spectrum.'
            good = 1
            return
        endif
        if newi gt self.nspec-1 then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the last spectrum.'
            good = 1
            return
        endif
    ;    good = (*self.science)[newi].good eq 1
    ;endwhile
    widget_control, widget_info(self.base, find_by_uname='filelist'), set_combobox_select=newi
    self.i = newi
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].telldiv,/nan)/median((*self.science)[self.i].telldiv)
    self->skyylim
    ;self->default_range
    
    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw

    common slit, bgood, rgood, bslit, rslit, bwidth, rwidth, bheight, rheight
    bgood = 0
    rgood = 0
    science = (*self.science)[self.i]
    if strmid(science.mask,0,4) eq '0024' then maskname = '0024' else maskname = science.mask
    if strmid(science.mask,0,5) eq '0024b' then maskname = science.mask
    bslitfile = file_dirname(science.spec1dfile, /mark_directory)+'slit.'+maskname+'.'+string(science.slit, format='(I03)')+'B.fits.gz'
    rslitfile = file_dirname(science.spec1dfile, /mark_directory)+'slit.'+maskname+'.'+string(science.slit, format='(I03)')+'R.fits.gz'
    widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
    wset, index
    loadct, 0, /silent
    tv, dblarr(1600, 300), 0, 0
    if file_test(bslitfile) then bslit = mrdfits(bslitfile, 1, /silent)
    if file_test(rslitfile) then rslit = mrdfits(rslitfile, 1, /silent)
    case 1 of
        file_test(bslitfile) and file_test(rslitfile): begin
            bflux = bslit.flux
            rflux = rslit.flux
            width = min([(size(bflux, /dimensions))[1], (size(rflux, /dimensions))[1]])
            flux = [bflux[*,0:width-1], rflux[*,0:width-1]]
        end
        file_test(bslitfile): flux = bslit.flux
        file_test(rslitfile): flux = rslit.flux
        else: return
    endcase
    med = median(flux)
    sig = stddev(flux)
    max_value = (med + (10 * sig)) < max(flux)
    min_value = (med - (2 * sig))  > min(flux)
    if (finite(min_value) EQ 0) then min_value = image_min
    if (finite(max_value) EQ 0) then max_value = image_max
    if (min_value GE max_value) then begin
        min_value = min_value - 1
        max_value = max_value + 1
    endif

    if file_test(bslitfile) then begin
        bflux = bytscl(bslit.flux, /nan, min=min_value, max=max_value, top=255) + 8
        bflux = -1.0*bflux + 255.0
        bheight = 1600./2047.*(size(bslit.dlambda, /dimensions))[1]
        if bheight gt 74 then begin
            bwidth = round(74./bheight*1600.)
            bheight = 74
        endif else begin
            bwidth = 1600
            bheight = round(bheight)
        endelse
        bflux1 = congrid(bflux[0:2047,*], bwidth, bheight, /interp)
        bflux2 = congrid(bflux[2048:4095,*], bwidth, bheight, /interp)
        tv, bflux1, 0, 225
        tv, bflux2, 0, 150
        bgood = 1
    endif else bgood = 0
    if file_test(rslitfile) then begin
        rflux = bytscl(rslit.flux, /nan, min=min_value, max=max_value, top=247) + 8
        rflux = -1.0*rflux + 255.0
        rheight = 1600./2047.*(size(rslit.dlambda, /dimensions))[1]
        if rheight gt 74 then begin
            rwidth = round(74./rheight*1600.)
            rheight = 74
        endif else begin
            rwidth = 1600
            rheight = round(rheight)
        endelse
        rflux1 = congrid(rflux[0:2047,*], rwidth, rheight, /interp)
        rflux2 = congrid(rflux[2048:4095,*], rwidth, rheight, /interp)
        tv, rflux1, 0, 75
        tv, rflux2, 0, 0
        rgood = 1
    endif else rgood = 0

    if bgood eq 1 and rgood eq 1 then begin
        n = n_elements(*self.tellstart)
        for i=0,n-1 do begin
            x1 = interpol(dindgen(4096), bslit.lambda0, (*self.tellstart)[i])
            offset1 = 0
            case 1 of
                x1 lt 0: x1 = -1
                x1 le 2047: offset1 = 225
                x1 le 4095: begin
                    offset1 = 150
                    x1 -= 2047
                end
                x1 gt 4095: begin
                    x1 = interpol(dindgen(4096), rslit.lambda0, (*self.tellstart)[i])
                    case 1 of
                        x1 lt 0: x1 = -1
                        x1 le 2047: offset1 = 75
                        x1 le 4095: begin
                            offset1 = 0
                            x1 -= 2047
                        end
                        x1 gt 4095: x1 = -1
                    endcase
                end
            endcase
            if x1 gt 0 then begin
                width1 = offset1 ge 150 ? bwidth : rwidth
                x1 *= width1 / 2047.
                height1 = offset1 ge 150 ? bheight : rheight
            endif else begin
                width1 = 0
                height1 = 0
            endelse

            x2 = interpol(dindgen(4096), bslit.lambda0, (*self.tellend)[i])
            offset2 = 0
            case 1 of
                x2 lt 0: begin
                    x2 = 2047
                    offset2 = offset1
                end
                x2 le 2047: offset2 = 225
                x2 le 4095: begin
                    offset2 = 150
                    x2 -= 2047
                end
                x2 gt 4095: begin
                    x2 = interpol(dindgen(4096), rslit.lambda0, (*self.tellend)[i])
                    case 1 of
                        x2 lt 0: begin
                            x2 = 2047
                            offset2 = offset1
                        end
                        x2 le 2047: offset2 = 75
                        x2 le 4095: begin
                            offset2 = 0
                            x2 -= 2047
                        end
                        x2 gt 4095: x2 = 2047
                    endcase
                end
            endcase
            if x2 gt 0 then begin
                width2 = offset2 ge 150 ? bwidth : rwidth
                x2 *= width2 / 2047.
                height2 = offset2 ge 150 ? bheight : rheight
            endif else begin
                width2 = 0
                height2 = 0
            endelse

            case 1 of
                x1 lt 0 and x2 lt 0: break
                x1 ge width1 and x2 ge width2: break
                offset1 eq offset2: begin
                    plots, [x1 > 0, x2 < width2], [0, 0]+height1-3+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [x1 > 0, x2 < width2], [0, 0]+2+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                end
                offset1 gt offset2: begin
                    plots, [x1 > 0, width2], [0, 0]+height1-3+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [x1 > 0, width2], [0, 0]+2+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [0, x2 < width2], [0, 0]+height2-3+offset2, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [0, x2 < width2], [0, 0]+2+offset2, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                end
                else:
            endcase
        endfor

        n = n_elements(*self.linewaves)
        for i=0,n-1 do begin
            x = interpol(dindgen(4096), bslit.lambda0, (*self.linewaves)[i] * (1d + science.zspec))
            case 1 of
                x lt 0: x = -1
                x le 2047: offset = 225
                x le 4095: begin
                    offset = 150
                    x -= 2047
                end
                x gt 4095: begin
                    x = interpol(dindgen(4096), rslit.lambda0, (*self.linewaves)[i] * (1d + science.zspec))
                    case 1 of
                        x lt 0: x = -1
                        x le 2047: offset = 75
                        x le 4095: begin
                            offset = 0
                            x -= 2047
                        end
                        x gt 4095: x = -1
                    endcase
                end
            endcase
            if x lt 0 then continue
            x *= (offset ge 150 ? bwidth : rwidth) / 2047.
            height = offset ge 150 ? bheight : rheight
            plots, [x, x], [0, 5]+offset, color=fsc_color((*self.linecolors)[i]), /device
            plots, [x, x], [height-6, height-1]+offset, color=fsc_color((*self.linecolors)[i]), /device
            xyouts, x+3, height-7+offset, (*self.linenames)[i], color=fsc_color((*self.linecolors)[i]), orientation=90, alignment=1, /device
        endfor
     endif
end


pro sps_fit::skyylim
    wsky = where((*self.science)[self.i].skylinemask ne -1, csky)
    if csky gt 0 then begin
        skyfit = (*self.science)[self.i].skyfit[wsky,*]
        yminmax = (*self.science)[self.i].dlam
        if skyfit[0,0] ne -1 then yminmax = [yminmax, skyfit[wsky,1]+skyfit[wsky,2], skyfit[wsky,1]-skyfit[wsky,2]]    
        self.skyylim = minmax(yminmax*2.35)    
    endif else begin
        self.skyylim = [0.9, 1.7]
    endelse
end


pro sps_fit::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.science)[self.i].telldiv,/nan)/median((*self.science)[self.i].telldiv)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self->skyylim
        self.lambdalim = (minmax((*self.science)[self.i].lambda / (1d + (*self.science)[self.i].zspec)) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d + (*self.science)[self.i].z)
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 1: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        mode eq 3: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.skyylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.skyylim[1], format='(D5.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
    endcase
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
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
    zl = mode lt 5 ? (*self.science)[self.i].z : 0.0
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


pro sps_fit::ylow, ylow, noredraw=noredraw
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 1: self.ylim[0] = ylow
        mode eq 3: self.skyylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 1: self.ylim[1] = yhigh
        mode eq 3: self.skyylim[1] = yhigh
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


; =============== TRACKING ==============
pro sps_fit::handle_tracking, ev
    uname = widget_info(ev.id, /uname)

    case (uname) of
        '2d': begin
            if ev.enter eq 0 then begin
                widget_control, widget_info(self.base, find_by_uname='obswave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='restwave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='yslit'), set_value='     '
            endif
        end
        else: 
    end
end


; ================ MOTION ===============
pro sps_fit::handle_motion, ev
    uname = widget_info(ev.id, /uname)

    case (uname) of
        '2d': begin
            common slit, bgood, rgood, bslit, rslit, bwidth, rwidth, bheight, rheight
            if bgood eq 0 or rgood eq 0 then return
            valid = 1
            case 1 of
                ev.y ge 225 and ev.y lt 225+bheight and ev.x lt bwidth: begin
                    slit = bslit
                    bheightorig = (size(bslit.dlambda, /dimensions))[1]
                    y = double(ev.y)-225.
                    if bheightorig gt 74 then y *= double(bheightorig) / 74.
                    x = double(ev.x) * 2047. / double(bwidth)
                end
                ev.y ge 150 and ev.y lt 150+bheight and ev.x lt bwidth: begin
                    slit = bslit
                    bheightorig = (size(bslit.dlambda, /dimensions))[1]
                    y = double(ev.y)-150.
                    if bheightorig gt 74 then y *= double(bheightorig) / 74.
                    x = (double(ev.x) * 2047. / double(bwidth)) + 2047.
                end
                ev.y ge 75 and ev.y lt 75+rheight and ev.x lt rwidth: begin
                    slit = rslit
                    rheightorig = (size(rslit.dlambda, /dimensions))[1]
                    y = double(ev.y)-75.
                    if rheightorig gt 74 then y *= double(rheightorig) / 74.
                    x = double(ev.x) * 2047. / double(rwidth)
                end
                ev.y ge 0 and ev.y lt rheight and ev.x lt rwidth: begin
                    slit = rslit
                    rheightorig = (size(rslit.dlambda, /dimensions))[1]
                    y = double(ev.y)
                    if rheightorig gt 74 then y *= double(rheightorig) / 74.
                    x = (double(ev.x) * 2047. / double(rwidth)) + 2047.
                end
                else: valid = 0
            endcase
            if ~valid then begin
                widget_control, widget_info(self.base, find_by_uname='obswave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='restwave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='yslit'), set_value='     '
            endif else begin
                lambda = slit.lambda0[x] + slit.dlambda[x, y]
                widget_control, widget_info(self.base, find_by_uname='obswave'), set_value='obs = '+string(lambda, format='(D7.1)')
                widget_control, widget_info(self.base, find_by_uname='restwave'), set_value='rest = '+string(lambda / (1.0 + (*self.science)[self.i].z), format='(D7.1)')
                widget_control, widget_info(self.base, find_by_uname='yslit'), set_value='y = '+string(y, format='(I3)')
            endelse
        end
        else: 
    end
end


; ============== COMBOBOX  =============
pro sps_fit::handle_combobox, ev
    ;widget_control, ev.id, get_uvalue=uvalue
    ;widget_control, ev.top, get_uvalue=obj
    self.i = ev.index
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].telldiv,/nan)/median((*self.science)[self.i].telldiv)
    self->skyylim
    self.keystate = 0
    self->redraw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; ============= DRAW CLICK =============
pro sps_fit::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 2 then click_coords /= 1d + (*self.science)[self.i].z
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
            zl = mode lt 2 ? (*self.science)[self.i].z : 0.0
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
            zl = mode lt 2 ? (*self.science)[self.i].z : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 0: ylim = self.ylim
                mode eq 1: ylim = self.skyylim
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
                mode le 1: begin
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


        1: begin        ;sky line fit
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

        2: begin        ;pixel mask
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
  for i=0,self.nspec-1 do begin
     self.i=i
     scienceall = *self.science
     science = scienceall[self.i]
     self->mask, science, nmc=0
     scienceall[self.i] = science
     ptr_free, self.science
     self.science = ptr_new(scienceall)
  endfor
  self.i = curi
  self->redraw
end

pro sps_fit::default_goodspec
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
       self.i=i
       science = scienceall[self.i]
       sn = science.sn
       z  = science.z
       if sn gt 6. and z gt 0.05 then science.good = 1 else science.good = 0
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

pro sps_fit::sn, science
    n = n_elements(science.lambda)
    wcont = where(science.contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.telldiv[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.telldiv[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])
end


; =============== MASK ================
pro sps_fit::mask, science, nomask=nomask, zfind=zfind, nozfind=nozfind, nmc=nmc
    specresfile = self.directory+'specres_poly.sav'
    globalres = 0
    if file_test(specresfile) then begin
        restore, specresfile
        dlam = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        globalres = 1
    endif

    if ~keyword_set(nomask) then begin

       ;1) remove emission lines
        linestart = *self.linestart
        lineend = *self.lineend
        linetype = *self.linetype
        wem = where(linetype eq 'e', cemlines)
        mask = bytarr(n_elements(science.lambda))+1
        w = where(science.contdivivar le 0 or ~finite(science.contdivivar), cw)
        if cw gt 0 then mask[w] = 0

        for i=0,cemlines-1 do begin
            w = where(science.lambda/(1d + science.zspec) ge linestart[wem[i]] and science.lambda/(1d + science.zspec) le lineend[wem[i]], c)
            if c gt 0 then mask[w] = 0
         endfor
       ;emission lines in Hb, Hg, Hd
       Hlines = [4863,4342,4103]
       Hwidth = [4,3,2]
       Hbreg = where(science.lambda/(1.+science.zspec) gt 4861.5 and  science.lambda/(1.+science.zspec) lt 4863.5,chbreg)
       ifem = where(science.contdiv(hbreg) gt 1.,cifem)
       if cifem eq chbreg then begin
          for i=0,n_elements(hlines)-1 do begin
             w = where(science.lambda/(1d + science.zspec) ge hlines[i]-Hwidth[i] and science.lambda/(1d + science.zspec) lt  hlines[i]+Hwidth[i], c)
             if c gt 0 then mask[w] = 0
          endfor
          science.goodsky = 0
       endif else science.goodsky = 1


        ;2) telluric region
        ;tellstart = [6864., 7591., 8938.]
        ;tellend = [6935., 7694., 9500.]
       ; tellstart = *self.tellstart
       ; tellend = *self.tellend
       ; tellend[4] = 9500. 
       ; for i=0,n_elements(tellstart)-1 do begin
       ;     w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
       ;     if c gt 0 then mask[w] = 0
       ;  endfor
        ;3) where it peaks greater than 3 sigmas
        ;w = where((science.contdiv-science.continuum) gt 3.*science.contdivivar^(-0.5), cw)
        ;if cw gt 0 then mask[w] = 0

        ;4) sky residuals
        ;first = deriv(science.lambda,science.contdiv)
        ;second = deriv(science.lambda,first)
        ;first  = first-mean(first)
        ;second = second-mean(second)
        ;thresh_first = 2.*stdev(first)
        ;thresh_second = 1.*stdev(second)
        ;skyres = where(abs(first) gt thresh_first and abs(second) gt thresh_second,cskyres)
        ;for i=0,cskyres-1 do begin
        ;   lambnow = science.lambda[skyres[i]]
        ;   w = where(science.lambda gt lambnow-2 and science.lambda lt lambnow+2,cw)
        ;   if cw gt 0 then mask[w]=0
        ;endfor
        ; 5) extreme noises
        w = where(science.contdiv gt 2. or science.contdiv lt -0.5, cw)
        nlambda = n_elements(science.lambda)
        if cw gt 0 then for i=0,cw-1 do mask[w[i]-2>0:w[i]+2<nlambda-1]=0

        ;6) where it's not finite
        w = where(finite(science.telldiv) eq 0 or finite(science.lambda) eq 0,c)
        if c gt 0 then mask[w]=0
        ;;

        ;7)where is is Mg region
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
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;8)wavelength range
        w = where(science.lambda/(1.+science.zspec) lt 3650. or science.lambda/(1.+science.zspec) gt 6000.,c)
        if c gt 0 then mask[w]=0
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        n = n_elements(science.lambda)
        nhalf = round(double(n) / 2.)
        mask[0:4] = 0
        mask[nhalf-4:nhalf+4] = 0
        mask[n-5:n-1] = 0
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
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        0: begin                ;continuum fit
           ;;making the highlight region
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
           ;; if cstart eq 0 or cend eq 0 then message, 'There are no continuum regions.'

           plot, science.lambda, science.telldiv/median(science.telldiv), xrange=self.lambdalim*(1d + science.zspec), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
           for i=0,cstart-1 do begin
              x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + science.zspec))) < (self.lambdalim[1] * (1d + science.zspec))
              y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
              polyfill, x, y, color=fsc_color('light cyan')
           endfor
           ;;finish highlighting masked regions
           ;;start ploting
           ;;plot the data and full model. Each normallized to around 1 independently
           goodlam = where(finite(science.telldiv))
           mediantelldiv = median(science.telldiv(goodlam))
           oplot, science.lambda, science.telldiv/mediantelldiv, color=fsc_color('black')
           oplot, science.lambda, science.continuum/mediantelldiv, color=fsc_color('green')
           if science.feh lt 3. then begin
              normfactor =median(science.spsspecfull)
              oplot, science.lambda,science.spsspecfull/normfactor,color=fsc_color('red')            
              oplot, science.lambda,science.spscontfull/normfactor,color=fsc_color('orange')            
           endif
           ;;make the axes labels
           plot, science.lambda, science.telldiv, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
           ;;highlight the telluric region and write line names
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
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


        1: begin        ;sky line fit
            wsky = where((*self.science)[self.i].skylinemask ne -1, csky)
            plot, science.lambda, 2.35*science.dlam, xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6sky line FWHM (!sA!r!u!9 %!6!n)!3', xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xrange=self.lambdalim * (1d + science.zspec), yrange=self.skyylim
            if csky gt 0 then begin
                fit = science.skyfit[wsky,*]
                if fit[0,0] ne -1 then begin
                    oploterror, fit[wsky,0], 2.35*fit[wsky,1], 2.35*fit[wsky,2], psym=1, color=fsc_color('black'), errcolor=fsc_color('black'), /nohat
                    wdel = where(science.skylinemask[wsky] eq 0, cdel)
                    if cdel gt 0 then plots, fit[wdel,0], 2.35*fit[wdel,1], color=fsc_color('red'), psym=7
                endif
            endif
            specresfile = self.directory+'specres_poly.sav'
            if file_test(specresfile) then begin
                restore, specresfile
                oplot, science.lambda, poly(science.lambda/1000 - 7.8, specres_poly), color=fsc_color('green')
            endif
        end
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        2: begin                ;rest frame
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
           plot, science.lambda/(1d + science.zspec), science.contdiv/science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
           
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
            nfull = n_elements(science.lambda)
            nhalf = fix(nfull/2)-1
            if science.wlfail_blue eq 1 then oplot, science.lambda[0:nhalf]/(1d +znow), science.contdiv[0:nhalf]/science.spscont[0:nhalf], color=fsc_color('darkgray') 
            if science.wlfail_red eq 1 then oplot, science.lambda[nhalf:nfull-1]/(1d +znow), science.contdiv[nhalf:nfull-1]/science.spscont[nhalf:nfull-1], color=fsc_color('darkgray') 

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
            plot, science.lambda/(1d + science.zspec), science.contdiv/science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6rest wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata, /noerase
            if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) then begin
               oplot, science.lambda/(1d + science.zspec), science.spsspecfull, color=fsc_color('red')
            endif
;	    oplot,science.lambda/(1.+science.zspec),science.contdiv*sqrt(science.contdivivar)/30.,color=fsc_color('pink')
         end
     endcase        
    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::statusbox, science=science
    if ~keyword_set(science) then science = (*self.science)[self.i]
    
    phot_color = science.phot_color
    case phot_color of
        'VK': begin
            color = science.v-science.k
            mag = science.k
            colorerr = sqrt((science.verr)^2. + (science.kerr)^2.)
            magerr = science.kerr
            noerror = science.verr le 0 or science.verr ge 10 or science.kerr le 0 or science.kerr ge 10
            clrlbl = 'V-K = '
            maglbl = 'K = '
        end
        'VJ': begin
            color = science.v-science.j
            mag = science.j
            colorerr = sqrt((science.verr)^2. + (science.jerr)^2.)
            magerr = science.jerr
            noerror = science.verr le 0 or science.verr ge 10 or science.jerr le 0 or science.jerr ge 10
            clrlbl = 'V-J = '
            maglbl = 'J = '
        end
        'VI': begin
            color = science.v-science.i
            mag = science.i
            colorerr = sqrt((science.verr)^2. + (science.ierr)^2.)
            magerr = science.ierr
            noerror = science.verr le 0 or science.verr ge 10 or science.ierr le 0 or science.ierr ge 10
            clrlbl = 'V-I = '
            maglbl = 'I = '
        end
        'VR': begin
            color = science.v-science.r
            mag = science.r
            colorerr = sqrt((science.verr)^2. + (science.rerr)^2.)
            magerr = science.rerr
            noerror = science.verr le 0 or science.verr ge 10 or science.rerr le 0 or science.rerr ge 10
            clrlbl = 'V-R = '
            maglbl = 'R = '
        end
        'BV': begin
            color = science.b-science.v
            mag = science.v
            colorerr = sqrt((science.berr)^2. + (science.verr)^2.)
            magerr = science.verr
            noerror = science.berr le 0 or science.berr ge 10 or science.verr le 0 or science.verr ge 10
            clrlbl = 'B-V = '
            maglbl = 'V = '
        end
        'BR': begin
            color = science.b-science.r
            mag = science.r
            colorerr = sqrt((science.berr)^2. + (science.rerr)^2.)
            magerr = science.rerr
            noerror = science.berr le 0 or science.berr ge 10 or science.rerr le 0 or science.rerr ge 10
            clrlbl = 'B-R = '
            maglbl = 'R = '
        end
        'RI': begin
            color = science.r-science.i
            mag = science.i
            colorerr = sqrt((science.rerr)^2. + (science.ierr)^2.)
            magerr = science.ierr
            noerror = science.rerr le 0 or science.rerr ge 10 or science.ierr le 0 or science.ierr ge 10
            clrlbl = 'R-I = '
            maglbl = 'I = '
        end
        'JK': begin
            color = science.j-science.k
            mag = science.k
            colorerr = sqrt((science.jerr)^2. + (science.kerr)^2.)
            magerr = science.kerr
            noerror = science.jerr le 0 or science.jerr ge 10 or science.kerr le 0 or science.kerr ge 10
            clrlbl = 'J-K = '
            maglbl = 'K = '
        end
        'F606WF814W': begin
            color = science.f606w-science.f814w
            mag = science.f814w
            colorerr = sqrt((science.f606werr)^2. + (science.f814werr)^2.)
            magerr = science.f814werr
            noerror = science.f606werr le 0 or science.f606werr ge 10 or science.f814werr le 0 or science.f814werr ge 10
            clrlbl = 'F606W-F814W = '
            maglbl = 'F814W = '
        end
        else: begin
            color = -999.0
            mag = -999.0
            colorerr = 0.0
            magerr = 0.0
            noerror = 1
            clrlbl = 'V-I = '
            maglbl = 'I = '
        end
    endcase

    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, science.good, science.goodfit]
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(science.objname, 2)+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='maglabel'), set_value=maglbl
    widget_control, widget_info(self.base, find_by_uname='collabel'), set_value=clrlbl
    widget_control, widget_info(self.base, find_by_uname='curmag'), set_value=mag gt 0 ? strcompress(string(mag, format='(D10.2)'), /rem)+(noerror ? '' : ' +/- '+strcompress(string(magerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curcol'), set_value=color gt -10 ? strcompress(string(color, format='(D10.2)'), /rem)+(noerror ? '' : ' +/- '+strcompress(string(colorerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(science.z, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzfit'), set_value=strcompress(string(science.zfit, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzquality'), set_value=strcompress(string(science.zquality, format='(D4.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=science.sn gt 0 ? strcompress(string(science.sn, format='(D10.1)'),/rem)+' ; '+strcompress(string(science.snfit, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curnloop'), set_value=science.nloop gt 0 ? strcompress(string(science.nloop, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curage'), set_value=science.age gt -100 ? strcompress(string(science.age, format='(D10.2)'), /rem)+(science.ageerr le 0 ? '' : ' +/- '+strcompress(string(science.ageerr, format='(D10.2)'), /rem))+' Gyr' : unknown
    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=science.logmstar gt 0 ? strcompress(string(science.logmstar, format='(D10.2)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfeh'), set_value=science.feh gt -100 ? strcompress(string(science.feh, format='(D10.2)'), /rem)+(science.feherr le 0 ? '' : ' +/- '+strcompress(string(science.feherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curcah'), set_value=science.caherr gt 0 ? strcompress(string(science.cah, format='(D10.2)'), /rem)+(science.caherr le 0 ? '' : ' +/- '+strcompress(string(science.caherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curgband'), set_value=science.gbanderr gt 0 ? strcompress(string(science.gband, format='(D10.2)'), /rem)+(science.gbanderr le 0 ? '' : ' +/- '+strcompress(string(science.gbanderr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curchisq'), set_value=science.chisq gt 0 ? strcompress(string(science.chisq, format='(D10.2)'), /rem) : unknown

    ;widget_control, widget_info(self.base, find_by_uname='curafe'), set_value=science.alphafe gt -100 ? strcompress(string(science.alphafe, format='(D10.2)'), /rem)+(science.alphafeerr le 0 ? '' : ' +/- '+strcompress(string(science.alphafeerr, format='(D10.2)'), /rem)) : unknown
    
    if science.vdisp gt 0 then widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=strcompress(string(science.vdisp, format='(D10.1)'), /rem)+(science.vdisperr le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr, format='(D10.1)'), /rem))+' km/s' else widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=unknown

    if science.vdisp_smm gt 0 then widget_control, widget_info(self.base, find_by_uname='curvdispsmm'), set_value=strcompress(string(science.vdisp_smm, format='(D10.1)'), /rem)+(science.vdisperr_smm le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr_smm, format='(D10.1)'), /rem))+' km/s (SMM)' else widget_control, widget_info(self.base, find_by_uname='curvdispsmm'), set_value=unknown
end


pro sps_fit::getscience, files=files
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'

    common mask_in, mask_in
    common mosfire, mosflag, deimosflag, lrisflag
   
    vdispfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451_sigmas.txt'
    readcol, vdispfile, vd_objname, vd_vdisp, vd_vdisperr, format='A,D,D', comment='#', /silent

    common npixcom, npix
    npix = 8192

    observatory, 'keck', obs
    sciencefits = self.directory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        c = n_elements(files)
        masks = strarr(c)
        slits = strarr(c)
        objnames = strarr(c)
        for i=0,c-1 do begin
            basefile = file_basename(files[i])
            extensions = strsplit(basefile, '.', /extract)
            masks[i] = extensions[1]
            slits[i] = extensions[2]
            objnames[i] = extensions[3]
        endfor
        w = where(strmatch(objnames, '*serendip*') eq 0 and strmatch(objnames, '*ACSstar*') eq 0, cw)
        if cw gt 0 then begin
            masks = masks[w]
            slits = slits[w]
            objnames = objnames[w]
            files = files[w]
        endif else begin
            message, 'No spectra found.'
        endelse

        nspec = n_elements(objnames)
        self.nspec = nspec
        speclist = masks+' '+strtrim(string(slits), 2)+' '+objnames
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1

        for i=0,nspec-1 do begin
            deimosflag = 1
            lrisflag = 0
            science = {science}
            science.skyfit = -1
            science.skylinemask = -1
            science.objname = objnames[i]
            science.mask = masks[i]
            science.slit = slits[i]

            data1 = mrdfits(files[i], 1, hdr, /silent) ;box car extraction
            fits_info,files[i],/silent,n_ext=n_ext
            phot = mrdfits(files[i], n_ext) ;photometries
            datafull = mrdfits(files[i], 2)
            if n_elements(datafull) gt 1 then begin; if the object is observed more than once 
                maxsn = max(datafull.sn,locmax)
                datafull = datafull[locmax]
            endif
            airmass = sxpar(hdr, 'AIRMASS')
            ras = sxpar(hdr, 'RA_OBJ')
            decs = sxpar(hdr, 'DEC_OBJ')
            decs_mask = sxpar(hdr, 'DEC')
            decssplit = strsplit(decs, ':', /extract)
            if stregex(decssplit[0], '\*', /boolean) eq 1 then begin
               decssplit_mask = strsplit(decs_mask, ':', /extract)
               decssplit[0] = decssplit_mask[0]
               decs = strjoin(decssplit, ':')
            endif
            if stregex(decssplit[1], '\*', /boolean) eq 1 then begin
               decssplit_mask = strsplit(decs_mask, ':', /extract)
               decssplit[0] = '-00'
               decssplit[1] = '00'
               decs = strjoin(decssplit, ':')
            endif
            get_coords, coords, instring=ras+'  '+decs
            science.ra = coords[0]*15.
            science.dec = coords[1]
            science.jdobs = double(sxpar(hdr, 'MJD-OBS')) + 2400000.5
            
            lambda = data1.lambda
            contdiv = data1.contdiv
            contdivivar = data1.contdivivar
            dlam = data1.dlam
	    npix = n_elements(science.lambda)
            if n_elements(lambda) gt npix then begin
               firstpix = where(lambda eq datafull.lambda[0],cfirstpix)
               if cfirstpix ne 1 then stop,'oops'
               ;take the blue side, this can be fixed for better handling in future
               lambda = lambda[firstpix:firstpix+npix-2]
               contdiv = contdiv[firstpix:firstpix+npix-2]
               contdivivar = contdivivar[firstpix:firstpix+npix-2]
               dlam = dlam[firstpix:firstpix+npix-2]
            endif  
            ;photometric stuff
            science.z = phot.z
            science.zquality = phot.zquality
            science.zspec = science.z
            science.zsource = phot.zsource
            science.b = phot.b_auto
            science.v = phot.v_auto
            science.r = phot.r_auto
            science.i = phot.i_auto
            science.j = phot.j_auto
            science.k = phot.k_auto
            science.f814w = phot.f814w_auto
            science.berr = phot.b_auto_err
            science.verr = phot.v_auto_err
            science.rerr = phot.r_auto_err
            science.ierr = phot.i_auto_err
            science.jerr = phot.j_auto_err
            science.kerr = phot.k_auto_err
            science.f814werr = phot.f814w_auto_err

            science.vdisp_smm = -999d
            science.vdisperr_smm = -999d
            if (size(vd_objname))[0] gt 0 then begin
                match, strtrim(vd_objname, 2), strtrim(science.objname, 2), w1, w2, count=cmatch
                if cmatch gt 0 then begin
                    science[w2].vdisp_smm = vd_vdisp[w1]
                    if strmid(mask_in, 0, 4) eq '0451' then begin
                        science[w2].vdisperr_smm = vd_vdisperr[w1]
                    endif else begin
                        science[w2].vdisperr_smm = -999d
                    endelse
                endif
            endif

            case 1 of
                science.b gt 10.0 and science.b lt 50.0 and science.v gt 10.0 and science.v lt 50.0: science.phot_color = 'BV'
                science.b gt 10.0 and science.b lt 50.0 and science.r gt 10.0 and science.r lt 50.0: science.phot_color = 'BR'
                science.v gt 10.0 and science.v lt 50.0 and science.i gt 10.0 and science.i lt 50.0: science.phot_color = 'VI'
                science.v gt 10.0 and science.v lt 50.0 and science.k gt 10.0 and science.k lt 50.0: science.phot_color = 'VK'
                science.v gt 10.0 and science.v lt 50.0 and science.j gt 10.0 and science.j lt 50.0: science.phot_color = 'VJ'
                science.j gt 10.0 and science.j lt 50.0 and science.k gt 10.0 and science.k lt 50.0: science.phot_color = 'JK'
                else: science.phot_color = 'BV'
            endcase

            science.spec1dfile = files[i]
            science.good = 1
            science.goodsky = 0
            w = where(science.age le 0, c)
            if c gt 0 then begin
                science[w].age = -999d
                science[w].ageerr = -999d
            endif
            science.feh = -999d
            science.feherr = -999d
            science.vdisp = -999d
            science.vdisperr = -999d

            self.i = i

            w = where(contdivivar gt 0 and finite(contdivivar), civar)
            if civar gt 10 then if w[civar-1] ne n_elements(contdivivar)-1 then contdivivar[w[civar-11:civar-1]] = 0

            n = n_elements(lambda)

            tell = mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent)
            wtell = n_elements(tell)-1
            tell = tell[wtell]
            ptr_free, self.tell
            self.tell = ptr_new(tell)

            t = (-1*ts_diff(lambda, 1))[0:n-2]
            wt = where(t le 0, ct)
            if ct gt 0 then begin
                message, 'Wavelength array for '+strtrim(objnames[i], 2)+' is not monotonic.  Omitting.', /info
                wgood[i] = 0
                continue
            endif

            science.lrisflag = lrisflag
            science.deimosflag = deimosflag

            science.lambda = lambda
            science.dlam = dlam
            science.contdiv = contdiv
            science.contdivivar = contdivivar

            science.contmask = datafull.contmask
            science.continuum = datafull.continuum
            science.skyfit = datafull.skyfit
            science.skylinemask = datafull.skylinemask
            science.goodsky = datafull.goodsky
            science.telldiv = contdiv*science.continuum
            science.telldivivar = contdivivar/(science.continuum)^2

            self->sn, science

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
        ;self->specres_mask, self.directory
        speclist = masks[wgood]+' '+strtrim(string(slits[wgood]), 2)+' '+objnames[wgood]
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=2
        ;self->fit_all
        self->writescience
     endif else begin ;if sps_fit.fits.gz exists or not
        scienceall = mrdfits(sciencefits, 1, /silent)

        self.nspec = n_elements(scienceall)
        speclist = scienceall.mask+' '+strtrim(string(scienceall.slit), 2)+' '+scienceall.objname
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        tell = (mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent))
        wtell = n_elements(tell)-1
        tell = tell[wtell]
        ptr_free, self.tell
        self.tell = ptr_new(tell)
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=2
    endelse
end


pro sps_fit::writescience
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing to database ...'
    scienceall = *self.science
    sciencefits = self.directory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits' : 'sps_fit.fits')
    mwrfits, scienceall, sciencefits, /create, /silent
    spawn, 'gzip -f '+sciencefits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::initialize_directory, directory=directory
    common mosfire, mosflag, deimosflag, lrisflag
    common mask_in, mask_in

    newdirectory = '/scr2/nichal/workspace2/sps_fit/data/'+mask_in+'/'
    if ~file_test(newdirectory) then file_mkdir, newdirectory

    if strmid(directory, 0, 1, /reverse_offset) ne '/' then directory = directory + '/'
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

    countfiles:
    files = file_search(directory, 'spec1d*.{fits,fits.gz}', count=c)
    sciencefits = newdirectory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    if c eq 0 then begin
       files = file_search(directory+'*/spec1d*.{fits,fits.gz}', count=c)
    endif
 
    if c eq 0 and ~file_test(sciencefits) then begin
        message, 'Unknown mask.'
    endif

    self.directory = newdirectory

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
    common sps_spec, sps, spsz, spsage
    if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)

    base = widget_base(/row, title='sps_fit', uvalue=self, mbar=menu, tab_mode=0, units=1)
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wdefaultrange = widget_button(tools_menu, value='Default Spectrum Settings', uname='default_range', uvalue='default_range')
    wdefault_cont = widget_button(tools_menu, value='Default Continuum Regions', uname='default_cont', uvalue='default_cont')
    wdefault_maskall = widget_button(tools_menu, value='Default Pixel Mask All', uname='default_maskall', uvalue='default_maskall')
    wdefault_goodspec = widget_button(tools_menu, value='Default Good Spectrum', uname='default_goodspec', uvalue='default_goodspec')
    wreprepare_all = widget_button(tools_menu, value='Fit All', uname='fit_all', uvalue='fit_all')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base

    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit','sky line fit', 'rest frame'], /column, /exclusive, set_value=2, uname='mode', uvalue='mode', /no_release)
    wstep = widget_base(wleft, /row, /align_center)
    wbackward = widget_button(wstep, value='<---', uvalue='backward', uname='backward', tab_mode=1, xsize=75)
    wforward = widget_button(wstep, value='--->', uvalue='forward', uname='forward', tab_mode=1, xsize=75)
    wprepbase = widget_base(wleft, /row, /align_center)
    wfit = widget_button(wprepbase, value='Fit', uvalue='fit', uname='fit', tab_mode=1, xsize=75)
    windicesbase = widget_base(wleft, /row, /align_center)
    windices = widget_button(windicesbase, value='Compute Indices', uvalue='indices', uname='indices', tab_mode=1, xsize=100)
    wdefault_mask = widget_button(windicesbase,value='Default Mask',uvalue='default_mask',uname='default_mask',tab_mode=1,xsize=100)
    wgoodbase = widget_base(wleft, /column, /align_center)
    wgood = cw_bgroup(wgoodbase, ['good sky', 'good spectrum','good fit'], /nonexclusive, set_value=[0, 0, 0], uname='good', uvalue='good')

    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wvbase = widget_base(wcurobj, /align_center, /row)
    wvlabel = widget_label(wvbase, value='B = ', /align_right, uname='maglabel', xsize=95)
    wcurmag = widget_label(wvbase, value='     ', /align_left, uname='curmag', uvalue='curmag', xsize=150)
    wbvbase = widget_base(wcurobj, /align_center, /row)
    wbvlabel = widget_label(wbvbase, value='B-V = ', /align_right, uname='collabel', xsize=95)
    wcurcol = widget_label(wbvbase, value='     ', /align_left, uname='curcol', uvalue='curcol', xsize=150)
    wagebase = widget_base(wcurobj, /align_center, /row)
    wagelabel = widget_label(wagebase, value='age = ', /align_right, uname='agelabel', xsize=95)
    wcurage = widget_label(wagebase, value='     ', /align_left, uname='curage', uvalue='curage', xsize=150)
    wmstarbase = widget_base(wcurobj, /align_center, /row)
    wmstarlabel = widget_label(wmstarbase, value='log M* = ', /align_right, uname='mstarlabel', xsize=95)
    wcurmstar = widget_label(wmstarbase, value='     ', /align_left, uname='curmstar', uvalue='curmstar', xsize=150)
    wfehbase = widget_base(wcurobj, /align_center, /row)
    wfehlabel = widget_label(wfehbase, value='[Fe/H] = ', /align_right, uname='fehlabel', xsize=95)
    wcurfeh = widget_label(wfehbase, value='     ', /align_left, uname='curfeh', uvalue='curfeh', xsize=150)
    wcahbase = widget_base(wcurobj, /align_center, /row)
    wcahlabel = widget_label(wcahbase, value='CaH = ', /align_right, uname='cahlabel', xsize=95)
    wcurcah = widget_label(wcahbase, value='     ', /align_left, uname='curcah', uvalue='curcah', xsize=150)
    wgbandbase = widget_base(wcurobj, /align_center, /row)
    wgbandlabel = widget_label(wgbandbase, value='Gband = ', /align_right, uname='gbandlabel', xsize=95)
    wcurgband = widget_label(wgbandbase, value='     ', /align_left, uname='curgband', uvalue='curgband', xsize=150)
    wvdispbase = widget_base(wcurobj, /align_center, /row)
    wvdisplabel = widget_label(wvdispbase, value='sigma_v = ', /align_right, uname='vdisplabel', xsize=95)
    wcurvdisp = widget_label(wvdispbase, value='     ', /align_left, uname='curvdisp', uvalue='curvdisp', xsize=150)

    wvdispsmmbase = widget_base(wcurobj, /align_center, /row)
    wvdispsmmlabel = widget_label(wvdispsmmbase, value='sigma_vsmm = ', /align_right, uname='vdispsmmlabel', xsize=95)
    wcurvdispsmm = widget_label(wvdispsmmbase, value='     ', /align_left, uname='curvdispsmm', uvalue='curvdispsmm', xsize=150)

    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='SNR = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)
    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='z = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)
    wzfitbase = widget_base(wcurobj, /align_center, /row)
    wzfitlabel = widget_label(wzfitbase, value='zfit = ', /align_right, uname='zfitlabel', xsize=95)
    wcurzfit = widget_label(wzfitbase, value='     ', /align_left, uname='curzfit', uvalue='curzfit', xsize=150)
    wzqualitybase = widget_base(wcurobj, /align_center, /row)
    wzqualitylabel = widget_label(wzqualitybase, value='zquality = ', /align_right, uname='zqualitylabel', xsize=95)
    wcurzquality = widget_label(wzqualitybase, value='     ', /align_left, uname='curzquality', uvalue='curzquality', xsize=150)
    wchisqbase = widget_base(wcurobj, /align_center, /row)
    wchisqlabel = widget_label(wchisqbase, value='chisq = ', /align_right, uname='chisqlabel', xsize=95)
    wcurchisq = widget_label(wchisqbase, value='     ', /align_left, uname='curchisq', uvalue='curchisq', xsize=150)
    wnloopbase = widget_base(wcurobj, /align_center, /row)
    wnlooplabel = widget_label(wnloopbase, value='nloop = ', /align_right, uname='nlooplabel', xsize=95)
    wcurnloop = widget_label(wnloopbase, value='     ', /align_left, uname='curnloop', uvalue='curnloop', xsize=150)

    wcursor = widget_base(wleft, /column, /align_right, tab_mode=0, /frame)
    wobswave = widget_label(wcursor, value='     ', /align_right, uname='obswave', uvalue='obswave', xsize=150)
    wrestwave = widget_label(wcursor, value='     ', /align_right, uname='restwave', uvalue='restwave', xsize=150)
    wyslit = widget_label(wcursor, value='     ', /align_right, uname='yslit', uvalue='yslit', xsize=150)

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

    widget_control, base, /realize
    self.base = base
    xmanager, 'sps_fit', self.base, /no_block, cleanup='sps_fit_cleanup'

    readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
    wbands = [1,2,3,4,5]
    tellstart = tellstart[wbands]
    tellend = tellend[wbands]
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick
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
             tell:ptr_new(), $
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
               slit:0L, $
               ra:-999d, $
               dec:-999d, $
               jdobs:-999d, $
               lambda:dblarr(npix), $
               telldiv:dblarr(npix), $
               telldivivar:dblarr(npix), $
               contmask:bytarr(npix), $
               continuum:dblarr(npix), $
               spscont:dblarr(npix), $
               contdiv:dblarr(npix), $
               contdivivar:dblarr(npix), $
               spsspec:dblarr(npix), $
               spsspecfull:dblarr(npix), $
               spscontfull:dblarr(npix), $
               fitmask:bytarr(npix), $
               dlam:dblarr(npix), $
               resscale:-999d, $
               skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
               skylinemask:lonarr(nsky), $
               goodsky:0, $
               airmass:-999d, $
               zobs:-999d, $
               zobserr:-999d, $
               zspec:-999d, $
               zfit:-999d, $
               z:-999d, $
               zquality:-999d, $
               zsource:0, $
               phot_color:'', $
               b:-999d, $
               v:-999d, $
               r:-999d, $
               i:-999d, $
               zmag:-999d, $
               j:-999d, $
               k:-999d, $
               f606w:-999d, $
               f814w:-999d, $
               berr:-999d, $
               verr:-999d, $
               rerr:-999d, $
               ierr:-999d, $
               jerr:-999d, $
               kerr:-999d, $
               f606werr:-999d, $
               f814werr:-999d, $
               age:-999d, $
               ageerr:-999d, $
               ageupper:-999d, $
               agelower:-999d, $
               feh:-999d, $
               feherr:-999d, $
               fehupper:-999d, $
               fehlower:-999d,$
               alphafe:-999d, $
	       alphafeerr:-999d, $
               alphafeupper:-999d, $
	       alphafelower:-999d, $
               logmstar:-999d, $
               b300:-999d, $
               vdisp:-999d, $
               vdisperr:-999d, $
               vdisp_smm:-999d, $
               vdisperr_smm:-999d, $
               cah:-999d, $
               caherr:-999d, $
               gband:-999d, $
               gbanderr:-999d, $
               chisq:-999d, $
               chisqmg:-999d, $
               chisqfe:-999d, $
               sn:-999d, $
               snfit:-999d, $
               spec1dfile:'', $
               good:0B,$
               goodfit:0B, $
               wlfail_blue:0B,$
               wlfail_red:0B,$
               nloop:0, $
               lrisflag: 0B, $
               deimosflag: 0B}
end


pro sps_fit_ms0451
    common mosfire, mosflag, deimosflag, lrisflag
    common mask_in, mask_in
    mask = 'combined_uniq_ms0451'
    mask_in = mask
    mosflag = 0
    deimosflag = 0
    lrisflag   = 0
    directory = '/scr2/nichal/keck/deimos/Cl0024MS0451/nicha_reduced/ms0451/'+mask
    if ~file_test(directory) then message, 'Mask not found.'
    n = obj_new('sps_fit', directory=directory)
end
