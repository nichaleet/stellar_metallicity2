pro stackspec::combine,spec,noredraw=noredraw
   common npixcom, npix,ndup
   common smoothpar, smooth_veldisp

   wgood = where(spec.indiv_good eq 1, cgood)
   if cgood eq 1 then begin
       spec.lambda = spec.indiv_lambda[*,wgood[0]]*(1.+spec.indiv_vhelio[wgood[0]]/3.e5)
       spec.spec = spec.indiv_spec[*,wgood[0]]
       spec.ivar = spec.indiv_ivar[*,wgood[0]]
       spec.sn = spec.indiv_sn[wgood[0]]
       spec.dlam = spec.indiv_dlam[wgood[0]]
   endif

   if cgood gt 1 then begin ;stack continuum and make output structure (strout)
       sn = spec.indiv_sn[wgood]
       ;shift velocity to heliocentric frame
       lambdaarr = spec.indiv_lambda[*,wgood] ;[npix,cgood]
       for k=0,cgood-1 do begin
           lambda = (1.+spec.indiv_vhelio[wgood[k]]/3e5)*spec.indiv_lambda[*,wgood[k]]
           lambdaarr[*,k] = lambda  ;ok it's shifted
       endfor
       ;find the total wl range (lambdafull)
       maxsn = max(sn,ref)
       ;use the lambda of the spec with max sn as ref
       lambdafull = lambdaarr[*,ref]
       spec.lambda = lambdafull
       fdlam = (smooth_veldisp/3.e5)*lambdafull/2.35
       spec.dlam = fdlam

       ;deal with the mask
       combmaskarr = bytarr(npix,cgood)+1
       for k=0,cgood-1 do begin
          woffmask = where(spec.indiv_combmask[*,wgood[k]] eq 0,coffmask)
          if coffmask gt 0 then begin
             loc = value_locate(lambdafull,lambdaarr[woffmask,k])
             loc = loc(where(loc ge 0))
             combmaskarr[loc,k] = 0
          endif
       endfor

       ;get data in and interpolate
       specarr=fltarr(npix,cgood)+1./0.
       ivararr = fltarr(npix,cgood)+1./0.
       skyarr = fltarr(npix,cgood)+1./0.
       framearr = bytarr(npix,cgood)+1
       for k=0,cgood-1 do begin
          loc = value_locate(lambdaarr[*,k],lambdafull) ;refvec,value
          woffrange = where(loc eq -1 or loc eq npix-1,coffrange,complement=winrange,$
                             ncomplement=cinrange)
          idlam = interpol(spec.indiv_dlam[*,wgood[k]],lambdaarr[*,k],lambdafull)
          specout = smooth_gauss_wrapper(lambdaarr[*,k],spec.indiv_spec[*,wgood[k]],$
                        lambdafull,sqrt(fdlam^2-idlam^2),ivar1=spec.indiv_ivar[*,wgood[k]],$
                        ivar2=ivarout)
          if coffrange gt 0 then begin
             specout[woffrange] = !values.f_infinity
             ivarout[woffrange] = 0
          endif 
          spec.indiv_spec_smoothed[*,wgood[k]] = specout
          spec.indiv_ivar_smoothed[*,wgood[k]] = ivarout 
          specarr[*,k] = specout
          ivararr[*,k] = ivarout

          framearr[woffrange,k] = 0
          woff = where(combmaskarr[*,k] eq 0,coff)
          if coff gt 0 then framearr[woff,k]=0
       endfor

       ;weighted average
       nframearr = total(framearr,2) ;size of npix
       wnodup = where(nframearr eq 1, cwnodup)
       if cwnodup gt 0 then begin
           for ii=0,cwnodup-1 do begin
              wframe = where(framearr[wnodup[ii],*] eq 1, cwframe)
              if cwframe ne 1 then stop,'oops something is wrong'
              spec.spec[wnodup[ii]] = specarr[wnodup[ii],wframe]
              spec.ivar[wnodup[ii]] = ivararr[wnodup[ii],wframe]
           endfor
       endif
       wdup = where(nframearr gt 1, cwdup)
       for ii=0,cwdup-1 do begin
           wnan = where(~finite(ivararr[wdup[ii],*]) or (ivararr[wdup[ii],*] eq 0) $
                   or ~finite(specarr[wdup[ii],*]) or (framearr[wdup[ii],*] eq 0),cnan)
           dspec = 1./sqrt(abs(ivararr[wdup[ii],*]))
           if cnan gt 0 then dspec(wnan) = 1./0.
           spec.spec[wdup[ii]]=wmean(specarr[wdup[ii],*],dspec,error=specerr,/nan)
           spec.ivar[wdup[ii]]=1./specerr^2
       endfor
       ;signal to noise
       contpara = poly_fit(spec.lambda,spec.spec,6,yfit=cont)
       ptr_free,self.cont
       self.cont = ptr_new(cont)

       ;calculate signal to noise
       dev = abs((spec.spec-cont)/cont)
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then begin
          sn_perpix = 1.0/mean(dev[w])
          spec.sn = sn_perpix
       endif else spec.sn = 0.
   endif
   self->statusbox, spec=spec

   if ~keyword_set(noredraw) then begin
      self->redraw
   endif
end

pro stackspec::indivsn,spec
   common npixcom, npix,ndup
   for i=0, ndup-1 do begin
       contpara = poly_fit(spec.indiv_lambda[*,i],spec.indiv_spec[*,i],6,yfit=cont)

       ;calculate signal to noise
       dev = abs((spec.indiv_spec[*,i]-cont)/cont)
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then sn_perpix = 1.0/mean(dev[w])
       spec.indiv_sn[i] = sn_perpix
   endfor
end

; =================
pro stackspec_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)
    case (name) of
	'BUTT': obj->handle_button, ev
	'TEXT': obj->handle_text, ev
	'TRAC': obj->handle_tracking, ev
	'DRAW': begin
	    if ev.type eq 0 then obj->handle_draw_click, ev
	    if ev.type eq 2 then obj->handle_motion, ev
	    if ev.type eq 5 then obj->handle_draw_key, ev
	end
	'DONE': widget_control, ev.top, /destroy
	else: begin
	    case (value) of
		'include': begin 
		    obj->toggle_good
		    obj->redraw
		end
                'iplotcont': begin
                    obj->statusbox
                    obj->redraw
                end
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


pro stackspec::handle_button, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    case (uvalue) of
        'back': self->newspec, increment=-1
        'next': self->newspec, increment=1
	'default_range': begin
	    self->default_range
	    self->redraw
	end
	'combine': begin
	   spec = *self.spec
	   self->combine, spec, /noredraw
	   ptr_free, self.spec
	   self.spec = ptr_new(spec)
	   self->redraw
	end
	'default_mask': self->default_mask
	'default_vhelio':self->default_vhelio
        'backward': self->step, -1
        'forward': self->step, 1
	'save': begin
	    spec = *self.spec
	    self->writespec,spec
	end
	'savespec':begin
	    spec = *self.spec
	    self->writespec,spec
	end
	'exit': begin
	    spec = *self.spec
	    self->writespec,spec
	    widget_control, ev.top, /destroy
	end
	else:
    endcase
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro stackspec::step, increment
    common npixcom, npix,ndup
    curi = self.indivi
    curi += increment
    if curi ge ndup then curi=0
    if curi lt 0 then curi = ndup-1

    self.indivi = curi

    self.keystate = 0
    self->statusbox
    self->redraw
end


pro stackspec::toggle_good
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
    spec = *self.spec
    widget_control, widget_info(self.base, find_by_uname='include'), get_value=good
    spec.indiv_good = good
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
    
end

pro stackspec::newspec, increment=increment, noredraw=noredraw
   common npixcom, npix,ndup
    if ~keyword_set(increment) then increment=0
    newi= self.i
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
    self.indivi = 0
    self->readspec
    spec = *self.spec
    ndup =  spec.indiv_count
    newindivlist = strcompress(indgen(ndup)+1,/rem)+'/'+strcompress(ndup, /rem)+' '+spec.indiv_objname 
    widget_control, widget_info(self.base, find_by_uname='indivfilelist'), set_value=newindivlist
    self->statusbox    
    self.ylim = [0,max((*self.spec).spec,/nan)]

    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw

end


pro stackspec::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.spec).indiv_spec,/nan)
        self.ylim[1] *= 1.1
        self.lambdalim = (minmax((*self.spec).lambda) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938.
    endif
    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0], format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1], format='(D7.1)'), /rem)
end

; ============== TEXT BOXES =============
pro stackspec::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
        'ylow': self->ylow, val
        'yhigh': self->yhigh, val
        'lambdalow': self->lambdalow, val
        'lambdahigh': self->lambdahigh, val
        'vhelio':self->editvhelio, val
        'z':self->editz,val
        else:
    end
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro stackspec::editvhelio,vhelio,noredraw=noredraw
   wplotcon = self.indivi
   spec = *self.spec
   spec.indiv_vhelio[wplotcon]=vhelio
   ptr_free, self.spec
   self.spec = ptr_new(spec)
   if ~keyword_set(noredraw) then self->redraw
end

pro stackspec::lambdalow, lambdalow, noredraw=noredraw
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->redraw
end


pro stackspec::lambdahigh, lambdahigh, noredraw=noredraw
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->redraw
end


pro stackspec::ylow, ylow, noredraw=noredraw
    self.ylim[0] = ylow
    if ~keyword_set(noredraw) then self->redraw
end


pro stackspec::yhigh, yhigh
    self.ylim[1] = yhigh
    self->redraw
end



; ============= DRAW CLICK =============
pro stackspec::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    case ev.modifiers of
        0: begin ;normal click
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
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0], format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1], format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        1: begin  ;shift and click
            if ev.press ne 1 then begin
                widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                return
            endif
            lrange = self.lambdalim[1] - self.lambdalim[0]
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0], format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1], format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin  ;control and click
            ylim = self.ylim
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
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(g8.1)'), /rem)
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
pro stackspec::handle_draw_key, ev
   common npixcom, npix,ndup 
   if ev.press ne 1 then return
   key = string(ev.ch)
   coords = convert_coord(ev.x, ev.y, /device, /to_data)
   spec= *self.spec
   wplotcon = self.indivi

   case key of
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
       's': begin
           case self.keystate of
               0: begin
                   self->writespec, spec
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
                   spec = *self.spec
                   w = where(spec.indiv_lambda[*,wplotcon] gt lambda1 and spec.indiv_lambda[*,wplotcon] lt lambda2, c)
                   if c eq 0 then begin
                       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                   endif else begin
                       spec.indiv_combmask[w,wplotcon] = 0
                       ptr_free, self.spec
                       self.spec = ptr_new(spec)
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
                   spec = *self.spec
                   w = where(spec.indiv_lambda[*,wplotcon] gt lambda1 and spec.indiv_lambda[*,wplotcon] lt lambda2, c)
                   if c eq 0 then begin
                       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                   endif else begin
                       spec.indiv_combmask[w,wplotcon] = 1
                       ptr_free, self.spec
                       self.spec = ptr_new(spec)
                       self->redraw
                   endelse
               end
               else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
           endcase
       end
       else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
   endcase
end


; ============== DESTROY ===============
pro stackspec_cleanup, ev
    widget_control, ev, get_uvalue=obj
    obj_destroy, obj
end


pro stackspec::cleanup
    ptr_free, self.spec
end

pro stackspec::default_mask
    wplotcon = self.indivi
    spec = *self.spec
    spec.indiv_combmask[*,wplotcon] = 1
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    self->redraw
end

pro stackspec::default_vhelio
    wplotcon = self.indivi
    spec = *self.spec
    spec.indiv_vhelio[wplotcon] = 0. 
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    self->redraw
    self->statusbox
end

; ============== REDRAW ===============
pro stackspec::redraw
   common npixcom, npix,ndup
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Redrawing ...'
   spec = *self.spec
   indivcolors = *self.indivcolors
   indivbgcolors = *self.indivbgcolors
   wplotcon = self.indivi

   ;plot zoom
    widget_control, widget_info(self.base, find_by_uname='zoom'), get_value=index
    wset, index
    plot,spec.lambda,spec.spec,xrange=[3700,3760],title='OII',position=[0.05,0.08,0.32,0.9],/normal, background=fsc_color('white'), color=fsc_color('black')
    oplot,spec.indiv_lambda[*,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5),spec.indiv_spec[*,wplotcon],color=fsc_color((*self.indivcolors)[0])
    oplot,spec.lambda,spec.spec,color=fsc_color('black'),thick=2

    plot,spec.lambda,spec.spec,xrange=[3900,4000],title='CaII',position=[0.37,0.08,0.64,0.9],/noerase,/normal, background=fsc_color('white'), color=fsc_color('black')
    oplot,spec.indiv_lambda[*,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5),spec.indiv_spec[*,wplotcon],color=fsc_color((*self.indivcolors)[0])
    oplot,spec.lambda,spec.spec,color=fsc_color('black'),thick=2

    plot,spec.lambda,spec.spec,xrange=[4830,4890],title='Hb',position=[0.69,0.08,0.96,0.9],/noerase,/normal, background=fsc_color('white'), color=fsc_color('black')
    oplot,spec.indiv_lambda[*,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5),spec.indiv_spec[*,wplotcon],color=fsc_color((*self.indivcolors)[0])
    oplot,spec.lambda,spec.spec,color=fsc_color('black'),thick=2

   ;plot main
   self->default_range, /update
   widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
   wset, index

   plot, spec.lambda, spec.spec, xrange=self.lambdalim, yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), /nodata
   ;;making the highlight region
   if wplotcon lt spec.indiv_count then begin
      t = round(-1*ts_diff(spec.indiv_combmask[*,wplotcon], 1))
      wstart = where(t eq 1, cstart)+1
      wend = where(t eq -1, cend)
      if spec.indiv_combmask[0,wplotcon] eq 1 then begin
         if cstart eq 0 then begin
            wstart = 0
         endif else begin
            wstart = [0, wstart]
           endelse
         cstart += 1
      endif
      if spec.indiv_combmask[n_elements(t)-1,wplotcon] eq 1 then begin
         if cend eq 0 then begin
            wend = n_elements(t)-1
         endif else begin
            wend = [wend, n_elements(t)-1]
         endelse
         cend += 1
      endif
      if cstart ne cend then message, 'There are a different number of starting and ending continuum wavelengths.'
      ;; if cstart eq 0 or cend eq 0 then message, 'There are no continuum regions.'

      for i=0,cstart-1 do begin
         multfactor = (1.+spec.indiv_vhelio[wplotcon]/3e5)
         x = ([spec.indiv_lambda[wstart[i],wplotcon]*multfactor,$
               spec.indiv_lambda[wstart[i],wplotcon]*multfactor, $
               spec.indiv_lambda[wend[i],wplotcon]*multfactor, $
               spec.indiv_lambda[wend[i],wplotcon]*multfactor] > (self.lambdalim[0])) < (self.lambdalim[1])
         dy = 0.05*(self.ylim[1]-self.ylim[0])
         y = [self.ylim[0]+dy, self.ylim[1]-dy, self.ylim[1]-dy, self.ylim[0]+dy]
         polyfill, x, y, color=fsc_color(indivbgcolors[0],/brewer)
      endfor
      ;;finish highlighting masked regions
   endif
   ;plotting
   oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
   oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')

   ;plot spectrum
   oplot,spec.indiv_lambda[*,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5),spec.indiv_spec[*,wplotcon],color=fsc_color(indivcolors[0])
   oplot,spec.lambda,spec.indiv_spec_smoothed[*,wplotcon],color=fsc_color(indivcolors[1])
   
   if total(spec.indiv_good) gt 1 then thick=2 else thick=1 
   oplot, spec.lambda, spec.spec, color=fsc_color('black'),thick=thick

   ;;highlighting the telluric bands and label line names
   n = n_elements(*self.linewaves)
   for i=0,n-1 do begin
      if (*self.linewaves)[i] le !X.CRANGE[0] or (*self.linewaves)[i] ge !X.CRANGE[1] then continue
      oplot, [(*self.linewaves)[i], (*self.linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
      xyouts, (*self.linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
   endfor

   widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0], format='(D7.1)'), /rem)
   widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1], format='(D7.1)'), /rem)

   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro stackspec::statusbox, spec=spec
    common npixcom, npix,ndup
    if ~keyword_set(spec) then spec = *self.spec
    unknown = '???'
    curi = self.indivi
    ndup = spec.indiv_count
    widget_control, widget_info(self.base, find_by_uname='include'), set_value=~spec.indiv_good[curi]
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(spec.objname, 2)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=strcompress(string(spec.sn, format='(D7.2)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='indivfilelist'), set_combobox_select=curi
    inditable = replicate({id:'',vhelio:-999,sn:-999,z:-999},5)
    fi = curi+4 < (ndup-1)
    nshow = fi-curi
    inditable[0:nshow].id = strcompress(spec.indiv_objname[curi:fi],/rem)
    inditable[0:nshow].vhelio = spec.indiv_vhelio[curi:fi]
    inditable[0:nshow].sn =  spec.indiv_sn[curi:fi]
    inditable[0:nshow].z = spec.indiv_z[curi:fi]
    widget_control,widget_info(self.base,find_by_uname='curinditable'),set_value = inditable
 
    widget_control, widget_info(self.base, find_by_uname='vhelio'), set_value=strcompress(string(spec.indiv_vhelio[curi], format='(D7.2)'), /rem)
end

pro stackspec::getspec, redo=redo
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common npixcom, npix,ndup
    common file, inputfile, massrange, outfits
    nspec = n_Elements(outfits)
    checkfile = file_test(outfits)
    self.nspec = nspec
    speclist = file_basename(outfits,'.sav')
    widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
    npix = 8192
    if (total(checkfile) ne nspec) or keyword_set(redo) then begin 
       ;read all the spec files
       data = mrdfits(inputfile,1)
       ;calculate sn per angstrom
       snperpix = data.sn
       snperang = snperpix
       for i=0,n_elements(snperpix)-1 do begin
          dlam = median(-1.*ts_diff(data[i].lambda,1)/(1.+data[i].z)) ;rest lambda
          snperang[i] = snperpix[i]/sqrt(dlam)  
       endfor
       ;get data for each mass range
       for i=0,n_elements(massrange)-2 do begin
          sel = where(data.logmstar ge massrange[i] and data.logmstar lt massrange[i+1] $
                      and snperang gt 3 and snperang le 15 and data.oiiew gt -5.,ndup)
          if ndup lt 3 then stop,'Warning: there are less than 3 galaxies in this mass range'
          print, outfits[i],ndup

          spec = {mask:'', $
           objname:'', $
           indiv_count:0L, $
           sn:-999d, $
           fileref:'', $
           lambda:dblarr(npix), $
           spec:dblarr(npix), $
           ivar:dblarr(npix), $
           dlam:dblarr(npix), $
           indiv_spec:dblarr(npix,ndup), $
           indiv_lambda:dblarr(npix,ndup), $
           indiv_ivar:dblarr(npix,ndup), $
           indiv_spec_smoothed:dblarr(npix,ndup), $
           indiv_ivar_smoothed:dblarr(npix,ndup), $
           indiv_combmask:bytarr(npix,ndup), $
           indiv_dlam:dblarr(npix,ndup), $
           indiv_sn:fltarr(ndup),$
           indiv_objname:strarr(ndup),$
           indiv_filename:strarr(ndup),$
           indiv_z:dblarr(ndup), $
           indiv_vhelio:dblarr(ndup), $
           indiv_good:bytarr(ndup)}

          objname = speclist[i] 
          spec.mask = 'spline'
          spec.fileref = inputfile
          spec.objname = objname
          spec.indiv_count = float(ndup)
          for iobs =0, ndup-1 do begin
             datanow = data[sel[iobs]]
             znow = datanow.zspec 
             if znow le 0. then znow = datanow.z
             spec.indiv_z[iobs] = znow
             spec.indiv_spec[*,iobs] = datanow.contdiv
             spec.indiv_lambda[*,iobs] = datanow.lambda/(1.+znow)
             spec.indiv_ivar[*,iobs] = datanow.contdivivar
             spec.indiv_dlam[*,iobs] = datanow.dlam/(1.+znow)
             spec.indiv_combmask[*,iobs]=1
             spec.indiv_objname[iobs] = datanow.objname
             spec.indiv_sn[iobs] = snperang[sel[iobs]]
             spec.indiv_vhelio[iobs] = 0. 
             spec.indiv_good[iobs] = 1
           endfor
         self.i = i
         self->combine,spec,/noredraw
         self->statusbox, spec=spec
         self->writespec,spec
      endfor
   endif 
   self.i = 0
end


pro stackspec::writespec,spec
    common file, inputfile, massrange, outfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing spectrum to database ...'
    specfits = outfits[self.i]
    ;mwrfits, spec, specfits, /create, /silent
    save,spec,filename=specfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro stackspec::readspec
    common file, inputfile, massrange, outfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading spectrum from database ...'
    specfits = outfits[self.i]
    if file_test(specfits) eq 0 then stop,'no file found'
    ;spec = mrdfits(specfits,1,/silent)
    restore, specfits
    ptr_free,self.spec
    self.spec = ptr_new(spec)    
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro stackspec::initialize_directory,directory=directory, redo=redo
   common file, inputfile, massrange, outfits

   newdirectory = '/scr2/nichal/workspace4/prepspec_ms0451/data/stacked/'
   if ~file_test(newdirectory) then file_mkdir, newdirectory

   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'
   
   nfiles = n_elements(massrange)-1
   outfits = strarr(nfiles)
   for i=0,nfiles-1 do outfits[i] = newdirectory+'stackedms0451_mass'+$
                                     strtrim(string(massrange[i],format='(f5.2)'),2)+'_'+$
                                     strtrim(string(massrange[i+1],format='(f5.2)'),2)+'.sav'
   self.directory = newdirectory

   self->getspec
   self.i = 0
   self->readspec
   spec = *self.spec
   self->statusbox, spec=spec
   self->default_range
   self.lambdalim = [3300, 7000]
   self->newspec
   widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

; =============== INIT ================
function stackspec::INIT,directory=directory
    common file, inputfile, massrange, outfits
    common npixcom, npix,ndup
    ndup = 15
    base = widget_base(/row, title='stack_spec_ms0451', uvalue=self, mbar=menu, tab_mode=0, units=1)
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wcombineall = widget_button(tools_menu, value = 'combine all spectra',uname='combine_all',uvalue='combine_all')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base
    ; ------ LEFT -------
    wprepbase = widget_base(wleft, /row, /align_center)
    wcombine = widget_button(wprepbase, value='Re-stack', uvalue='combine', uname='combine', tab_mode=1, xsize=100)
    wsavespec = widget_button(wprepbase,value='Save spec', uvalue='savespec',uname='savespec',tab_mode=1,xsize=100)
    windivbase = widget_base(wleft,/row,/align_center)
    wdefaultmaskbase = widget_base(windivbase,/row,/align_center)
    wdefaultmask = widget_button(wdefaultmaskbase,value='Default Mask',uvalue='default_mask',uname='default_mask',tab_mode=1,xsize=120)    
    wdefaultvheliobase = widget_base(windivbase,/row,/align_center)
    wdefaultvhelio = widget_button(wdefaultvheliobase,value='Reset Vhelio',uvalue='default_vhelio',uname='default_vhelio',tab_mode=1,xsize=120)

    windivfile = widget_base(wleft, /row, /align_center)
    windivfilelist = widget_combobox(windivfile, uname='indivfilelist', value='            ', tab_mode=1, /dynamic_resize)
    wstep = widget_base(wleft, /row, /align_center)
    wbackward = widget_button(wstep, value='<---', uvalue='backward', uname='backward', tab_mode=1, xsize=75)
    wforward = widget_button(wstep, value='--->', uvalue='forward', uname='forward', tab_mode=1, xsize=75)
    wincludebase = widget_base(wleft, /row, /align_center,/frame)
    wincluderadio = cw_bgroup(wincludebase,["include","discard"],/exclusive,set_value=0,uname='include',uvalue='include',label_top='include',/row)
    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='final s/n = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)


    windibase = widget_base(wcurobj,/align_left,/row,xsize=400)
    winditable = widget_table(windibase,value=replicate({id:'',vhelio:-999,sn:-999,z:-999},5),/row_major,column_labels=['id','vhelio','sn','redshift'],uname='curinditable',uvalue='curinditable')

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
    wvheliocontrol = widget_base(wspeccontrol,/frame,/row)
    wvheliolabel = widget_label(wvheliocontrol,value='vehelio:',/align_center,uname='vheliolabel',font='-urw-standard symbols l-medium-r-normal--0-0-0-0-p-0-adobe-symbol')
    wvhelio = widget_text(wvheliocontrol,xsize=8,/editable,uname='vhelio',uvalue='vhelio')

    widget_control, base, /realize
    self.base = base
    xmanager, 'stackspec', self.base, /no_block, cleanup='stackspec_cleanup'

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
    self.indivcolors = ptr_new(['blueviolet','blue','seagreen','springgreen','gold','chocolate','firebrick','saddlebrown'])
    self.indivbgcolors = ptr_new(['pbg1','blu1','ryb5','grn1','wt4','org1','red2','tg4'])
    common random, seed
    seed = systime(1)

    self->initialize_directory,directory=directory
    return, 1
end


pro stackspec__define
    state = {stackspec, $
             base:0L, $
             directory:'', $
             specinfo:ptr_new(), $
             spec:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
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
             keystate:0, $
             lambda1:0d, $
             indivcolors:ptr_new(),$
             indivbgcolors:ptr_new(),$
             cont:ptr_new(),$
             i:0L,$
             indivi:0L}

end

pro stackspec_ms0451
   common file, inputfile, massrange, outfits
   common npixcom, npix,ndup
   common smoothpar, smooth_veldisp
   inputfile = '/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit01.fits.gz'
   smooth_veldisp = 350 ;km/s
   if file_test(inputfile) eq 0 then stop,'Stopped: cannot find input file'
   massrange=[9.,10.,10.5,10.7,10.82,11.0,11.5]
   n = obj_new('stackspec',directory=directory)

end
