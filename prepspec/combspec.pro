pro combspec::newspec, increment=increment, noredraw=noredraw
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
    self.ylim = minmax((*self.spec)[self.i].indiv_spec,/nan)/median((*self.spec)[self.i].indiv_spec)

    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw
end


pro combspec::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.spec)[self.i].indiv_telldiv)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self.lambdalim = (minmax((*self.spec)[self.i].lambda / (1d + (*self.spec)[self.i].z)) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d + (*self.spec)[self.i].z)
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 0: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
    endcase
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 1 ? (*self.spec)[self.i].z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
end

pro combspec::indivcontinuum, spec
   common npixcom, npix,ndup
   common continuum, conttype
   nobs = spec.indiv_count
   for i=0,nobs-1 do begin
      contmask = spec.indiv_contmask[*,i]
      lambda = spec.indiv_lambda[*,i]
      telldivivar = spec.indiv_telldivivar[*,i]
      won = where(contmask eq 1, complement=woff, con)
      if con gt 100 then begin
         if conttype eq 'polycon' then begin
            degree = 7
            coeff =  poly_fit(lambda[won],spec[won],degree,measure_errors=1./sqrt(telldivivar[won]),status=status)
            if status ne 0 then stop,'stopped:(continuumnormalize_deimos.pro) blueside has some error'
            cont = poly(lambda,reform(coeff))
         endif else begin
            bkspace = 500.
            bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=telldivivar[won], bkspace=bkspace, upper=3, lower=3, /silent)
            if bkpt[0] eq -1 then stop,'cannot do continnum fit'
            cont = slatec_bvalu(lambda, bkpt, coeff)
         endelse
         check = where(finite(cont),ccheck)
         if ccheck lt 10 then stop
         contdiv = spec.indiv_telldiv[*,i]/cont
         contdivivar = spec.indiv_telldivivar[*,i]*cont^2
      endif else begin
         contdiv = dblarr(npix)
         cont = dblarr(npix)
         contdivivar = dblarr(npix)
         spec.indiv_good[i] = 0  
      endelse
      spec.indiv_continuum[*,i] = cont
      spec.indiv_contdivivar[*,i] = contdivivar
      spec.indiv_contdiv[*,i] = contdiv

      wcont = where(spec.indiv_contmask[3:npix-4,i] eq 1)+3
      wcont = wcont[where(finite(spec.indiv_telldiv[wcont,i]) and finite(spec.indiv_continuum[wcont,i]) and spec.indiv_continuum[wcont,i] ne 0)]
      dev = abs((spec.indiv_telldiv[wcont,i] - spec.indiv_continuum[wcont,i]) / spec.indiv_continuum[wcont,i])
      avgdev = mean(dev)
      w = where(dev lt 3.0*avgdev, c)
      if c gt 0 then spec.indiv_sn[i] = 1.0/mean(dev[w])
   endif
end

pro combspec::indivmask,spec
   nobs = spec.indiv_count
   for i=0,nobs-1 do begin
      lambda = spec.indiv_lambda[*,i]
      lambdarest = lambda/(1.+spec.z)
      spec = spec.indiv_telldiv[*,i]
      ivar = spec.indiv_ivar[*,i]
      telldivivar = spec.telldivivar[*,i]
   
      ;create continuum mask
      npix = n_elements(lambda)
      contmask = bytarr(npix)+1
      linestart = *self.linestart
      lineend = *self.lineend

      for i=0,n_elements(linestart)-1 do begin
         w = where(lambdarest ge linestart[i] and lambdarest le lineend[i], c)
         if c gt 0 then contmask[w] = 0
      endfor
   
      w = where(~finite(spec) or ~finite(lambda) or ~finite(ivar) or ~finite(telldivivar)  or ivar eq 0, c)
      if c gt 0 then contmask[w]=0
   
      tellmask = bytarr(n_elements(lambda))
      tellstart = [6864., 7591., 8938.]
      tellend = [6935., 7694., 9500.]
      for j=0,n_elements(tellstart)-1 do begin
          w = where(lambda ge tellstart[j] and lambda le tellend[j], c)
          if c gt 0 then begin
              if j eq 0 then contmask[w] = 0
              tellmask[w] = 1
          endif
      endfor
      nhalf = fix(npix/2)
      contmask[nhalf-3:nhalf+3]=0
      contmask[0:2] = 0
      contmask[npix-3:npix-1] = 0
      spec.indiv_contmask[*,i] = contmask
      spec.indiv_tellmask[*,i] = tellmask
   endfor
end

pro combspec::indivtelluric,spec
   nobs = spec.indiv_count
   for i=0,nobs-1 do begin
      aratio = spec.indiv_airmass[i]/((*self.tell).airmass)
      telllambda = (*self.tell)[0].lambda
      tellspec = ((*self.tell)[0].spec)^aratio
      tellivar = ((*self.tell)[0].ivar)*(((*self.tell)[0].spec)/(aratio*tellspec))^2.
      ivarmissing = 10d10
      w = where(tellivar ge 1d8, c) ;these are not in the telluric band
      if c gt 0 then tellivar[w] = ivarmissing
      f = where(finite(tellspec) and tellspec gt 0 and finite(tellivar) and tellivar gt 0 and tellivar lt 1d8) ;where the band is
      telllambda = telllambda[f]
      tellspec = tellspec[f]
      tellivarfull = tellivar
      tellivar = tellivar[f]

      tellspecnew = interpolate(tellspec, findex(telllambda, spec.indiv_lambda[*,i]), missing=1.)
      tellivarnew = interpolate(tellivarfull, findex(telllambda, spec.indiv_lambda[i,*]), missing=ivarmissing)

      spec.indiv_telldiv[*,i] = spec.indiv_spec[*,i] / tellspecnew
      spec.indiv_telldivivar[*,i] = (spec.indiv_ivar[*,i]^(-1.) + (spec.indiv_telldiv[*,i])^2.*tellivarnew^(-1.))^(-1.) * tellspecnew^2.
   endfor
end

function combspec::fitskylines,spec,iobs
    lambda = spec.indiv_lambda[*,iobs]
    skyspec = spec.skyspec[*,iobs]
    w = where(~finite(skyspec), c)
    if c gt 0 then skyspec[w] = 0
    skyspec = skyspec/max(skyspec)
    medskyspec = median(skyspec)

    deriv1skyspec = deriv(lambda, skyspec)
    deriv2skyspec = deriv(lambda, deriv1skyspec)

    thresh = 1.5
    nlines = 1000
    while nlines gt 200 do begin
        w = where(abs(deriv1skyspec) lt 0.2 and deriv2skyspec lt -0.01 and skyspec gt thresh*medskyspec)
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

        skyspecmax = max(skyspecfit, wmax)

        guess = [skyspecmax, lambdafit[wmax], 0, medskyspec]
        yfit = gaussfit(lambdafit, skyspecfit, a, estimates=guess, sigma=aerr, nterms=4, chisq=chisq)
        sigma[i] = abs(a[2])
        sigmaerr[i] = aerr[2]
        linelambda[i] = a[1]

        if chisq gt 1d-3 or sigmaerr[i] gt 0.8 then wgoodline[i] = 0
    endfor
    wgood = where(sigma gt 0 and wgoodline eq 1, c)
    if c gt 0 then begin
        linelambda = linelambda[wgood]
        sigma = sigma[wgood]
        sigmaerr = sigmaerr[wgood]
        return, [[linelambda], [sigma], [sigmaerr]]
    endif else return, [[-1], [-1], [-1]]
end

pro combspec::indivspecres, spec
   nobs = spec.indiv_count
   for i=0,nobs-1 do begin
      lambda = spec.indiv_lambda[*,i]
      wsky = where(spec.skylinemask[*,i] ne -1, cw)
      if cw eq 0 then begin
         fit = self->fitskylines(spec,i)
         w = where(2.35*fit[*,1] gt 0.8 and 2.35*fit[*,1] lt 7.0, cwprev)
         n = (size(fit))[1]
         if n gt 200 then message, 'Too many sky lines!'
         if n gt 0 then spec.indiv_skyfit[0:n-1,*,i] = fit
         if (n le 3) or (cwprev lt 3) then begin
            spec.indiv_dlam[*,i] = replicate(1.37/2.35, n_elements(lambda))
            spec.indiv_skylinemask[*,i] = lonarr(n_elements(spec.indiv_skylinemask[*,i]))-1
            spec.indiv_goodsky[i] = 0
            message, 'Unstable sky line fit.  Using FWHM = 1.37 A', /info
            return
         endif

         ;quadratic fit
         qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, $
                       measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)
         ;remove outliers and refit 4 times
         for j=0,4 do begin
            wnow = where(abs(2.35*fit[w,1] - yfit) lt 2.*2.35*fit[w,2], cw) ;good sigmas
            if cw eq cwprev then break
            cwprev = cw
            if cw lt 3 then begin
               spec.indiv_goodsky[i] = 0
               message, 'The spectral resolution fit is very poor.', /info
               break
            endif
            w = w[wnow]
            qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], $
                          chisq=chisq, /double, yfit=yfit)
         endfor
         n = (size(fit))[1]
         spec.skylinemask[*,i] = 0
         spec.skylinemask[w,i] = 1
         if n lt 200 then spec.skylinemask[n:n_elements(spec.skylinemask[*,i])-1,i] = -1
      endif else begin
         wsky = where(spec.skylinemask[*,i] eq 1, csky)
         if csky lt 3 then begin
            spec.indiv_dlam[*,i] = replicate(1.37/2.35, n_elements(lambda))
            spec.indiv_goodsky[i] = 0
            message, 'Too few arc lines for new fit.', /info
            return
         endif
         fit = science.skyfit[wsky,*,i]
         qf = poly_fit(fit[wsky,0]/1000.0 - 7.8, 2.35*fit[wsky,1], 2, measure_errors=2.35*fit[wsky,2],$
                       chisq=chisq, /double, yfit=yfit)
      endelse
      l = lambda / 1000. - 7.8
      dlam = poly(l, qf)
      dlam /= 2.35
      spec.indiv_dlam[*,i] = dlam
      spec.goodsky[i] = 1
   endfor
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

           plot, science.lambda, science.contdiv*science.spscontfull, xrange=self.lambdalim*(1d + science.zspec), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
           for i=0,cstart-1 do begin
              x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + science.zspec))) < (self.lambdalim[1] * (1d + science.zspec))
              y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
              polyfill, x, y, color=fsc_color('light cyan')
           endfor
           ;;finish highlighting masked regions
           ;;start ploting
           ;;plot the data(*spscontfull) and full model
           oplot, science.lambda, science.contdiv*science.spscontfull, color=fsc_color('black')
           if science.feh ne -999. then begin
              oplot, science.lambda,science.spsspecfull,color=fsc_color('red')
              oplot, science.lambda,science.spscontfull,color=fsc_color('orange')
           endif
           ;;make the axes labels
           plot, science.lambda, science.contdiv*science.spscontfull, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
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

        1: begin                ;rest frame
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
            if (science.feh gt -10 and science.age gt 0.0) or total(science.spsspec gt 0.0) then begin
               oplot, science.lambda/(1d + science.zspec), science.spsspecfull, color=fsc_color('red')
            endif
;           oplot,science.lambda/(1.+science.zspec),science.contdiv*sqrt(science.contdivivar)/30.,color=fsc_color('pink')
         end
     endcase
    zl = mode lt 2 ? (*self.science)[self.i].z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::statusbox, spec=spec
    common npixcom, npix,ndup
    if ~keyword_set(spec) then spec = (*self.spec)[self.i]
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='include'), set_value=spec.indiv_good
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(spec.objname, 2)+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='curcat'), set_value=strcompress(spec.specname)
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(spec.z, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=spec.logmstar gt 0 ? strcompress(string(spec.logmstar, format='(D10.2)'), /rem) : unknown
    inditable = replicate({mask:'',id:'',vhelio:-999,sn:-999,exptime:-999},ndup)
    inditable.mask = strcompress(spec.indiv_mask)
    inditable.id = strcompress(spec.indiv_objname)
    inditable.vhelio = spec.indiv_vhelio
    inditable.sn =  spec.indiv_sn
    inditable.exptime = spec.indiv_exptime
    widget_control,widget_info(self.base,find_by_uname='curinditable'),set_value = inditable
end


pro combspec::getspec, list=list, cat=cat, distance=distance, mass=mass
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common continuum, conttype
    common npixcom, npix,ndup

    npix = 8192

    observatory, 'keck', obs
    specfits = self.directory+'combspec_out.fits.gz'

    if ~file_test(specfits) then begin 
       nspec = n_elements(list)
       self.nspec = nspec

       specall = replicate({spec}, nspec)

       masks = strarr(nspec)
       slits = strarr(nspec)
       objnames = strarr(nspec)
       for i=0,nspec-1 do begin
          masks[i] = self.directory
          slits[i] = strtrim(string(i),2)
          objnames[i] = strtrim(cat[i].specname,2)
          spec.mask = masks[i]
          spec.slit = slits[i]
          spec.objname = objnames[i]
       endfor
       speclist = masks+' '+slits+' '+objnames
       widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

       for i=0,nspec-1 do begin
           spec = {spec}
           obj = list[i]
           
           struct_assign,cat[i],spec,/nozero ;copy cat over
           spec.mask = masks[i]
           spec.slit = slits[i]

           nobs = obj.count
           for iobs =0, nobs=1 do begin
               strblu = mrdfits(strtrim(obj.filename[iobs],2),1,hdrblu,/silent,status=status1)
               strred = mrdfits(strtrim(obj.filename[iobs],2),2,hdrred,/silent,status=status2)
               if status1 ne 0 or status2 ne 0 then continue
               if typename(strblu) eq 'INT' then continue
               spec.indiv_spec[*,iobs] = [strblu.spec,strred.spec]
               spec.indiv_lambda[*,iobs] = [strblu.lambda,strred.lambda]
               spec.indiv_ivar[*,iobs] = [strblu.ivar,strred.ivar]
               spec.indiv_skyspec[*,iobs] = [strblu.skyspec,strred.skyspec]
               spec.indiv_objpos[iobs] = strblu.objpos
               spec.indiv_fwhm[iobs] = strblu.fwhm
               spec.indiv_nsigma[iobs] = strblu.nsigma
               spec.indiv_R1[iobs] = strblu.R1
               spec.indiv_R2[iobs] = strblu.R2
               spec.indiv_ivarfudge[iobs] = strblu.ivarfudge
               spec.indiv_objname[iobs] = obj.objname[iobs]
               spec.indiv_rastring[iobs] = obj.rastring[iobs]
               spec.indiv_decstring[iobs] = obj.decstring[iobs]
               spec.indiv_ra[iobs] = obj.ra[iobs]
               spec.indiv_dec[iobs] = obj.dec[iobs]
               spec.indiv_mask[iobs] = obj.mask[iobs]
               spec.indiv_slit[iobs] = obj.slit[iobs]
               spec.indiv_filename[iobs] = obj.filename[iobs]
               spec.indiv_exptime[iobs] = obj.exptime[iobs]
               spec.indiv_count = obj.count[iobs]
               spec.indiv_airmass[iobs] = sxpar(hdrblu,'airmass')
               spec.indiv_good[iobs] = 1
           endfor    
           self->indivspecres, spec
           self->indivtelluric, spec
           self->indivmask,spec
           self->indivcontinuum,spec
           self->statusbox, spec=spec
           specall[i] = spec
      endfor
      self.i = 0
      ptr_free, self.spec
      self.spec = ptr_new(specall)
      self.nspec = cgood
      self->writespec
      speclist = masks[wgood]+' '+strtrim(string(slits[wgood]), 2)+' '+objnames[wgood]
      widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
      widget_control, widget_info(self.base, find_by_uname='mode'), set_value=1
      self->writespec
   endif else begin
      specall = mrdfits(specfits, 1, /silent)

      self.nspec = n_elements(specall)
      speclist = specall.mask+' '+strtrim(string(specall.slit), 2)+' '+specall.objname
      widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

      tell = (mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent))
      wtell = n_elements(tell)-1
      tell = tell[wtell]
      ptr_free, self.tell
      self.tell = ptr_new(tell)
      ptr_free, self.spec
      self.spec = ptr_new(specall)
      widget_control, widget_info(self.base, find_by_uname='mode'), set_value=1
   endelse
end

pro sps_fit::writespec
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing to database ...'
    specall = *self.spec
    specfits = self.directory+'combspec_out.fits.gz' 
    mwrfits, specall, specfits, /create, /silent
    spawn, 'gzip -f '+specfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro combspec::initialize_directory, directory=directory,rematch=rematch

   common continuum, conttype

   if keyword_set(rematch) then begin
      cat = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451master_specz.fits.gz',1)
      badcat = where(cat.ra eq 99., complement=goodcat)
      cat = cat(goodcat)
      list = mrdfits('matched_ms0451obj.fits',1)
      ;match object in the list and catalog 
      ras = list.ra[0]
      decs = list.dec[0]
      spherematch, cat.ra,cat.dec,ras,decs,2./3600.,w1,w2,distance
      cat = cat(w1)
      list = list(w2)
      ;only take the galaxies with z from 0.52 to 0.56
      goodz = where(cat.z gt 0.5 and cat.z lt 0.6,cgoodz)
      print,'total of',cgoodz, ' objects between z=0.5 and z=0.6'
      cat = cat(goodz)
      list = list(goodz)
      distance = distance(goodz)
      ;calculate mass
      mass = getms0451mass(cat)
      stop
      ;writefits
      mwrfits,list,'matched_ms0451obj_members.fits',/create,/silent
      mwrfits,cat,'matched_ms0451obj_members.fits',/silent
      mwrfits,distance,'matched_ms0451obj_members.fits',/silent
      mwrfits,mass,'matched_ms0451obj_members.fits',/silent
      makems0451regionfile,list,cat
   endif
   list = mrdfits('matched_ms0451obj_members.fits',1)
   cat  = mrdfits('matched_ms0451obj_members.fits',2)
   distance = mrdfits('matched_ms0451obj_members.fits',3)
   mass = mrdfits('matched_ms0451obj_members.fits',4)


   newdirectory = '/scr2/nichal/workspace4/prepspec/data/'+mask_in+'/'
   if ~file_test(newdirectory) then file_mkdir, newdirectory

   if strmid(directory, 0, 1, /reverse_offset) ne '/' then directory = directory + '/'
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

   countfiles:
   specfits = newdirectory+'combspec_out.fits.gz'

   self.directory = newdirectory

   self->getspec, list=list, cat=cat, distance=distance, mass=mass
   self.i = 0
   science = (*self.science)[self.i]
   self->statusbox, science=science
   self->default_range
   self.lambdalim = [3300, 7000]
   self->newspec, increment=0
   widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

; =============== INIT ================
function combspec::INIT, directory=directory
    common npixcom, npix,ndup

    base = widget_base(/row, title='combine_spec', uvalue=self, mbar=menu, tab_mode=0, units=1)
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wdefaultrange = widget_button(tools_menu, value='Default Spectrum Settings', uname='default_range', uvalue='default_range')
    wdefault_cont = widget_button(tools_menu, value='Default Continuum Regions', uname='default_cont', uvalue='default_cont')
    wdefault_maskall = widget_button(tools_menu, value='Default Pixel Mask All', uname='default_maskall', uvalue='default_maskall')
    wdefault_goodspec = widget_button(tools_menu, value='Default Good Spectrum', uname='default_goodspec', uvalue='default_goodspec')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base
    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit', 'rest frame'], /column, /exclusive, set_value=1, uname='mode', uvalue='mode', /no_release)
    wprepbase = widget_base(wleft, /row, /align_center)
    wfitcont = widget_button(wprepbase, value='Fit Continuum', uvalue='fitcont', uname='fitcont', tab_mode=1, xsize=100)
    wcombine = widget_button(wprepbase, value='Re-stack', uvalue='combine', uname='combine', tab_mode=1, xsize=100)
    wdefaultmaskbase = widget_base(wleft,/row,/align_center)
    wdefaultmask = widget_button(wdefaultmaskbase,value='Default Mask',uvalue='default_mask',uname='default_mask',tab_mode=1,xsize=100)
    wincludebase = widget_base(wleft, /column, /align_center,/frame)
    wincluderadio = cw_bgroup(wincludebase,['0','1','2','3','4','5','6'],/nonexclusive,set_value=[1,1,1,1,1,1,1],uname='include',uvalue='include')
    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wcatbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    wcatlabel = widget_label(wcatbase, value='cat ID:', /align_right, uname='catlabel', xsize=65)
    wcurcat = widget_label(wcatbase, value='     ', /align_left, uname='curcat', uvalue='curcat', xsize=180)

    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='redshift = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)

    wmstarbase = widget_base(wcurobj, /align_center, /row)
    wmstarlabel = widget_label(wmstarbase, value='log M* = ', /align_right, uname='mstarlabel', xsize=95)
    wcurmstar = widget_label(wmstarbase, value='     ', /align_left, uname='curmstar', uvalue='curmstar', xsize=150)
    windibase = widget_base(wcurobj,/align_left,/row,xsize=400)
    winditable = widget_table(windibase,value=replicate({mask:'',id:'',vhelio:-999,sn:-999,exptime:-999},7),/row_major,column_labels=['mask','id','vhelio','sn','exptime'],uname='curinditable',uvalue='curinditable')

    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wback = widget_button(wfile, value='Back', uvalue='back', uname='back', tab_mode=1)
    wfilelist = widget_combobox(wfile, uname='filelist', value='                 ', tab_mode=1, /dynamic_resize)
    wnext = widget_button(wfile, value='Next', uvalue='next', uname='next', tab_mode=1)
    wstatus = widget_text(wfile, xsize=108, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)

    wspec = widget_base(wright, /frame, /column)
    wspecplot = widget_draw(wspec, xsize=1400, ysize=400, uname='spec', /button_events, keyboard_events=1)
    wzoomplot = widget_draw(wspec, xsize=1400, ysize=300, uname='zoom')
    
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
    xmanager, 'combspec', self.base, /no_block, cleanup='combspec_cleanup'

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

    common random, seed
    seed = systime(1)

    self.i = -1
    self->initialize_directory, directory=directory
    return, 1
end


pro combspec__define
    state = {combspec, $
             base:0L, $
             directory:'', $
             spec:ptr_new(), $
             tell:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
             divylim:[-100d, 100d], $
             linewaves:ptr_new(), $
             linenames:ptr_new(), $
             linecolors:ptr_new(), $
             tellstart:ptr_new(), $
             tellend:ptr_new(), $
             tellthick:ptr_new(), $
             linestart:ptr_new(), $
             lineend:ptr_new(), $
             linetype:ptr_new(), $
             nspec:0L, $
             i:0L, $
             keystate:0, $
             lambda1:0d}
end

pro spec__define
   common npixcom, npix,ndup
   ndup = 7
   spec = {spec, $
           mask:'', $
           slits:'', $
           lambda:dblarr(npix), $
           dlam:dblarr(npix), $
           contdiv:dblarr(npix), $
           contdivivar:dblarr(npix), $
           sn:-999d, $
           exptime:-999d, $
           indiv_spec:dblarr(npix,ndup), $
           indiv_lambda:dblarr(npix,ndup), $
           indiv_ivar:dblarr(npix,ndup), $
           indiv_contmask:dblarr(npix,ndup), $
           indiv_telmask:dblarr(npix,ndup), $
           indiv_objpos:fltarr(ndup), $ 
           indiv_fwhm:fltarr(ndup), $
           indiv_nsigma:fltarr(ndup), $
           indiv_R1:lonarr(ndup), $
           indiv_R2:lonarr(ndup), $
           indiv_skyspec:dblarr(npix,ndup), $
           indiv_ivarfudge:fltarr(ndup),$
           indiv_dlam:dblarr(npix,ndup), $
           indiv_skyfit:dblarr(200,3,ndup), $
           indiv_skylinemask:lonarr(200,ndup), $
           indiv_airmass:dblarr(ndup), $
           indiv_telldiv:dblarr(npix,ndup), $
           indiv_telldivivar:dblarr(npix,ndup), $
           indiv_continuum:dblarr(npix,ndup), $
           indiv_contdiv:dblarr(npix,ndup), $
           indiv_contdivivar:dblarr(npix,ndup), $
           indiv_sn:fltarr(ndup),$
           indiv_objname:strarr(ndup),$
           indiv_rastring:strarr(ndup),$
           indiv_decstring:strarr(ndup),$
           indiv_ra:dblarr(ndup), $
           indiv_dec:dblarr(ndup), $
           indiv_mask:strarr(ndup),$
           indiv_slit:strarr(ndup),$
           indiv_filename:strarr(ndup),$
           indiv_exptime:dblarr(ndup), $
           indiv_good:bytarr(ndup), $
           indiv_count:-999L,$
           NUM:-999L,$   
           RA:-999d ,$   
           DEC:-999d ,$   
           SUBARU_RA:-999d ,$   
           SUBARU_DEC:-999d ,$   
           FUV_AUTO:-999d ,$   
           NUV_AUTO:-999d ,$   
           B_AUTO:-999d ,$   
           V_AUTO:-999d ,$   
           R_AUTO:-999d ,$   
           I_AUTO:-999d ,$   
           J_AUTO:-999d ,$   
           K_AUTO:-999d ,$   
           F814W_AUTO:-999d ,$   
           FUV_AUTO_ERR:-999d,$  
           NUV_AUTO_ERR:-999d,$  
           B_AUTO_ERR:-999d ,$  
           V_AUTO_ERR:-999d ,$  
           R_AUTO_ERR:-999d ,$  
           I_AUTO_ERR:-999d ,$  
           J_AUTO_ERR:-999d ,$  
           K_AUTO_ERR:-999d ,$  
           F814W_AUTO_ERR:-999d,$
           B_ISO:-999d ,$  
           V_ISO:-999d ,$  
           R_ISO:-999d ,$  
           I_ISO:-999d, $      
           J_ISO:-999d, $      
           K_ISO:-999d, $      
           B_ISO_ERR:-999d, $  
           V_ISO_ERR:-999d, $  
           R_ISO_ERR:-999d, $  
           I_ISO_ERR:-999d, $  
           J_ISO_ERR:-999d, $  
           K_ISO_ERR:-999d, $  
           HST_X:-999d, $      
           HST_Y:-999d, $      
           HST_POS:-999, $     
           HST_ID:-999, $      
           MORPH:-999, $       
           Z:-999d, $                                         
           ZQUALITY:-999d, $   
           ZSOURCE:-999d, $    
           AXIS_RATIO:-999d, $ 
           PA:-999d, $         
           SPECNAME:'', $               '
           MORPH_COMMENT:'', $                        '
           STAR:-999d, $       
           MASS_KCORRECT:-999d}

end

pro combspec,conttype=conttype
   common continuum, conttype
   n = obj_new('combspec',directory=directory)
end
