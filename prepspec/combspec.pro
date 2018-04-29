pro combspec::indivcontinuum, spec

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

pro combspec::getspec, list=list, cat=cat, distance=distance, mass=mass
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'

    common interpolation, interptype
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
       for i=0,c-1 do begin
          masks[i] = self.directory
          slits[i] = strtrim(string(i),2)
          objnames[i] = strtrim(cat[i].specname,2)
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
             endfor    
             self->indivspecres, spec
             self->indivtelluric, spec
             self->indivmask,spec
             self->indivcontinuum,spec




pro combspec::initialize_directory, directory=directory,rematch=rematch
   common interpolation, interptype

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
    winditable = widget_table(windibase,value=replicate({mask:'',id:'',vhelio:-999,sn:-999,exptime:-999},7),/row_major,column_labels=['mask','id','vhelio','sn','exptime'])

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

pro combspec,interptype=interptype
   common interpolation, interptype
   n = obj_new('combspec',directory=directory)
end
