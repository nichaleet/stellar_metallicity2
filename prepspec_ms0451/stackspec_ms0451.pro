pro stackspec::combine,spec,noredraw=noredraw
   common npixcom, npix,ndup
   wgood = where(spec.indiv_good eq 1, cgood)
   if cgood eq 1 then begin
       spec.lambda = spec.indiv_lambda[*,wgood[0]]*(1.+spec.indiv_vhelio[wgood[0]]/3.e5)
       spec.spec = spec.indiv_spec[*,wgood[0]]
       spec.ivar = spec.indiv_ivar[*,wgood[0]]
       spec.sky = spec.sky[*,wgood[0]]
       spec.exptime = spec.indiv_exptime[wgood[0]]
       spec.sn = spec.indiv_sn[wgood[0]]
   endif

   if cgood gt 1 then begin ;stack continuum and make output structure (strout)
       spec.exptime = total(spec.indiv_exptime[wgood])
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
          woffrange = where(loc eq -1 or loc eq npix-1,coffrange,complement = winrange)
          specarr[winrange,k]=interpol(spec.indiv_spec[*,wgood[k]],lambdaarr[*,k],lambdafull[winrange])
          ivararr[winrange,k]=interpol(spec.indiv_ivar[*,wgood[k]],lambdaarr[*,k],lambdafull[winrange])
          skyarr[winrange,k]=interpol(spec.indiv_sky[*,wgood[k]],lambdaarr[*,k],lambdafull[winrange])
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
              spec.sky[wnodup[ii]] = skyarr[wnodup[ii],wframe]
           endfor
       endif
       wdup = where(nframearr gt 1, cwdup)
       for ii=0,cwdup-1 do begin
           wnan = where(~finite(ivararr[wdup[ii],*]) or ~finite(specarr[wdup[ii],*]) or (framearr[wdup[ii],*] eq 0),cnan)
           dspec = 1./sqrt(abs(ivararr[wdup[ii],*]))
           if cnan gt 0 then dspec(wnan) = 1./0.
           err = abs(dspec/specarr[wdup[ii],*])
           spec.sky[wdup[ii]] = wmean(skyarr[wdup[ii],*],abs(err*skyarr[wdup[ii],*]),/nan)
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
       if c gt 0 then sn_perpix = 1.0/mean(dev[w])
       spec.sn = sn_perpix
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
	'save': begin
	    spec = *self.spec
	    self->writespec,spec
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



pro stackspec::toggle_good
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
    spec = *self.spec
    widget_control, widget_info(self.base, find_by_uname='include'), get_value=good
    spec.indiv_good = good
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
    
end

pro stackspec::newspec, noredraw=noredraw

    self->readspec
    self->statusbox    
    self.ylim = minmax((*self.spec).indiv_spec,/nan)

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

pro stackspec::editz,z
    self.z=z
    self->redraw
end

pro stackspec::editvhelio,vhelio,noredraw=noredraw
   widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
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
   widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
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
    widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
    spec = *self.spec
    spec.indiv_combmask[*,wplotcon] = 1
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    self->redraw
end

pro stackspec::default_vhelio
    widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
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
   widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon       
   znow = self.z

   ;plot zoom
    widget_control, widget_info(self.base, find_by_uname='zoom'), get_value=index
    wset, index
    z = self.z
    plot,spec.lambda/(1.+z),spec.spec,xrange=[3700,3760],title='OII',position=[0.05,0.08,0.32,0.9],/normal, background=fsc_color('white'), color=fsc_color('black')
    for iplot=spec.indiv_count-1,0,-1 do oplot,spec.indiv_lambda[*,iplot]/(1.+z)*(1.+spec.indiv_vhelio[iplot]/3e5),spec.indiv_spec[*,iplot],color=fsc_color((*self.indivcolors)[iplot])
    oplot,spec.lambda/(1.+z),spec.spec,color=fsc_color('black'),thick=2

    plot,spec.lambda/(1.+z),spec.spec,xrange=[3900,4000],title='CaII',position=[0.37,0.08,0.64,0.9],/noerase,/normal, background=fsc_color('white'), color=fsc_color('black')
    for iplot=spec.indiv_count-1,0,-1 do oplot,spec.indiv_lambda[*,iplot]/(1.+z)*(1.+spec.indiv_vhelio[iplot]/3e5),spec.indiv_spec[*,iplot],color=fsc_color((*self.indivcolors)[iplot])
    oplot,spec.lambda/(1.+z),spec.spec,color=fsc_color('black'),thick=2

    plot,spec.lambda/(1.+z),spec.spec,xrange=[4830,4890],title='Hb',position=[0.69,0.08,0.96,0.9],/noerase,/normal, background=fsc_color('white'), color=fsc_color('black')
    for iplot=spec.indiv_count-1,0,-1 do oplot,spec.indiv_lambda[*,iplot]/(1.+z)*(1.+spec.indiv_vhelio[iplot]/3e5),spec.indiv_spec[*,iplot],color=fsc_color((*self.indivcolors)[iplot])
    oplot,spec.lambda/(1.+z),spec.spec,color=fsc_color('black'),thick=2

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
         x = ([spec.indiv_lambda[wstart[i],wplotcon]*multfactor, spec.indiv_lambda[wstart[i],wplotcon]*multfactor, $
               spec.indiv_lambda[wend[i],wplotcon]*multfactor, $
               spec.indiv_lambda[wend[i],wplotcon]*multfactor] > (self.lambdalim[0])) < (self.lambdalim[1])
         dy = 0.05*(self.ylim[1]-self.ylim[0])
         y = [self.ylim[0]+dy, self.ylim[1]-dy, self.ylim[1]-dy, self.ylim[0]+dy]
         polyfill, x, y, color=fsc_color(indivbgcolors[wplotcon],/brewer)
      endfor
      ;;finish highlighting masked regions
   endif
   ;plotting
   oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
   oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')

   for idup=0,ndup-1 do begin
      if spec.indiv_good[idup] eq 1 then begin
         oplot,spec.indiv_lambda[*,idup]*(1.+spec.indiv_vhelio[idup]/3e5),spec.indiv_spec[*,idup],color=fsc_color(indivcolors(idup))
      endif
   endfor
   
   if total(spec.indiv_good) gt 1 then thick=2 else thick=1 
   oplot, spec.lambda, spec.spec, color=fsc_color('black'),thick=thick

   ;;highlighting the telluric bands and label line names
   n = n_elements(*self.tellstart)
   for i=0,n-1 do begin
      oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
   endfor
   n = n_elements(*self.linewaves)
   for i=0,n-1 do begin
      z = self.z
      if (*self.linewaves)[i]*(1.+z) le !X.CRANGE[0] or (*self.linewaves)[i]*(1.+z) ge !X.CRANGE[1] then continue
      oplot, [(*self.linewaves)[i]*(1.+z), (*self.linewaves)[i]*(1.+z)], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
      xyouts, (*self.linewaves)[i]*(1.+z)+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
   endfor



   widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0], format='(D7.1)'), /rem)
   widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1], format='(D7.1)'), /rem)

   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro stackspec::statusbox, spec=spec
    common npixcom, npix,ndup
    if ~keyword_set(spec) then spec = *self.spec
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='include'), set_value=spec.indiv_good
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(spec.objname, 2)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=strcompress(string(spec.sn, format='(D7.2)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='curexptime'), set_value=strcompress(string(spec.exptime, format='(I8)'), /rem)

    inditable = replicate({id:'',vhelio:-999,sn:-999,exptime:-999},ndup)
    inditable.id = strcompress(spec.indiv_objname)
    inditable.vhelio = spec.indiv_vhelio
    inditable.sn =  spec.indiv_sn
    inditable.exptime = spec.indiv_exptime
    widget_control,widget_info(self.base,find_by_uname='curinditable'),set_value = inditable
 
    widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
    widget_control, widget_info(self.base, find_by_uname='vhelio'), set_value=strcompress(string(spec.indiv_vhelio[wplotcon], format='(D7.2)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='z'), set_value=strcompress(string(self.z, format='(D5.3)'), /rem)
end


pro stackspec::getspec, list=list, redo=redo
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common npixcom, npix,ndup
    common mask_in, mask_in, listfile

    npix = 2834
    objname = (strsplit(listfile,'ms',/extract,/regex))[0]
    specfits = self.directory+objname+'.fits'
    datadir = '/scr2/nichal/workspace_lowmass/dbsp_data/'
    if ~file_test(specfits) or keyword_set(redo) then begin 
       ndup = n_elements(list)
       spec = {spec}
       spec.mask = mask_in
       spec.objname = objname
       spec.indiv_count = ndup
       for iobs =0, ndup-1 do begin
           filenow = datadir+mask_in+'/'+list[iobs]
           if file_test(filenow) eq 0 then begin
               stop, 'cannot find file '+filenow
               continue
           endif
           array = readfits(filenow,hdr,/silent)
           spec.indiv_spec[*,iobs] = array[*,0,0]
           spec.indiv_lambda[*,iobs] = (findgen(npix)-sxpar(hdr,'CRPIX1')+1)*sxpar(hdr,'CD1_1')+sxpar(hdr,'CRVAL1')
           spec.indiv_ivar[*,iobs] = 1./(array[*,0,3])^2
           spec.indiv_sky[*,iobs] = array[*,0,2]
           spec.indiv_combmask[*,iobs]=1
           spec.indiv_airmass[iobs] = float(sxpar(hdr,'airmass'))
           spec.indiv_jd[iobs] = double(sxpar(hdr,'jd'))
           radec = stringradec2deci(sxpar(hdr,'ra'), sxpar(hdr,'dec'))
           spec.indiv_ra[iobs] = radec[0] 
           spec.indiv_dec[iobs] = radec[1]
           spec.indiv_objname[iobs] = sxpar(hdr,'OBJECT')
           spec.indiv_filename[iobs] = filenow
           spec.indiv_exptime[iobs] = float(sxpar(hdr,'exptime'))
           spec.indiv_vhelio[iobs] = 0. 
           spec.indiv_good[iobs] = 1
      endfor
      self->indivsn,spec
      self->combine,spec,/noredraw
      self->statusbox, spec=spec
      self->writespec,spec
      ptr_free, self.spec
      self->readspec
   endif else begin
      spec = mrdfits(specfits, 1, /silent)
      ptr_free,self.spec
      self.spec = ptr_new(spec)
      self.i = 0
      self->readspec
   endelse
end


pro stackspec::writespec,spec
    common mask_in, mask_in, listfile
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing spectrum to database ...'
    specfits = self.directory+(strsplit(listfile,'ms',/extract,/regex))[0]+'.fits'

    mwrfits, spec, specfits, /create, /silent
    spawn, 'gzip -f '+specfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro stackspec::readspec
    common mask_in, mask_in, listfile
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading spectrum from database ...'
    specfits = self.directory+(strsplit(listfile,'ms',/extract,/regex))[0]+'.fits.gz'
    if file_test(specfits) eq 0 then stop,'no file found'
    spec = mrdfits(specfits,1,/silent)
    ptr_free,self.spec
    self.spec = ptr_new(spec)    
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'

end

pro stackspec::initialize_directory, directory=directory,redo=redo

   common mask_in, mask_in, listfile
   datadir = '/scr2/nichal/workspace_lowmass/dbsp_data/'
   readcol,datadir+mask_in+'/'+listfile,list,format='A'

   newdirectory = '/scr2/nichal/workspace_lowmass/dbsp_data/reduced_data/'+mask_in+'/'
   if ~file_test(newdirectory) then file_mkdir, newdirectory

   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

   specfits = newdirectory+(strsplit(listfile,'ms',/extract,/regex))[0]+'.fits'
   self.directory = newdirectory

   self->getspec, list=list, redo=redo
   self->readspec
   spec = *self.spec
   self->statusbox, spec=spec
   self->default_range
   self.lambdalim = [3300, 7000]
   self->newspec
   widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

; =============== INIT ================
function stackspec::INIT, inputfile=inputfile
    common file, inputspecfile, massrange
    common npixcom, npix,ndup

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
    wincludebase = widget_base(wleft, /row, /align_center,/frame)
    wincluderadio = cw_bgroup(wincludebase,strtrim(string(indgen(ndup)),2),/column,/nonexclusive,set_value=bytarr(ndup)+1,uname='include',uvalue='include',label_top='include')
    wplotconradio = cw_bgroup(wincludebase,strtrim(string(indgen(ndup)),2),/column,/exclusive,set_value=0,uname='iplotcont',uvalue='iplotcont',label_top='plot')
    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='s/n = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)

    wexptimebase = widget_base(wcurobj, /align_center, /row)
    wexptimelabel = widget_label(wexptimebase, value='Exp time = ', /align_right, uname='exptimelabel', xsize=95)
    wcurexptime = widget_label(wexptimebase, value='     ', /align_left, uname='curexptime', uvalue='curexptime', xsize=150)

    windibase = widget_base(wcurobj,/align_left,/row,xsize=400)
    winditable = widget_table(windibase,value=replicate({id:'',redshift:-999,mass:-999,sn:-999},ndup),/row_major,column_labels=['id','redshift','mass','sn'],uname='curinditable',uvalue='curinditable')

    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
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

    wzcontrol = widget_base(wspeccontrol,/frame,/row)
    wzlabel = widget_label(wzcontrol,value='z:',/align_center,uname='zlabel',font='-urw-standard symbols l-medium-r-normal--0-0-0-0-p-0-adobe-symbol')
    wz = widget_text(wzcontrol,xsize=8,/editable,uname='z',uvalue='z')

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

    self->initialize_directory, directory=directory
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
             z:0d, $
             cont:ptr_new()}
end

pro spec__define
   common npixcom, npix,ndup
   spec = {spec, $
           mask:'', $
           objname:'', $
           lambda:dblarr(npix), $
           spec:dblarr(npix), $
           ivar:dblarr(npix), $
           sky:dblarr(npix), $
           sn:-999d, $
           exptime:-999d, $
           indiv_spec:dblarr(npix,ndup), $
           indiv_lambda:dblarr(npix,ndup), $
           indiv_ivar:dblarr(npix,ndup), $
           indiv_sky:bytarr(npix,ndup), $
           indiv_combmask:bytarr(npix,ndup), $
           indiv_airmass:fltarr(ndup), $ 
           indiv_jd:fltarr(ndup), $
           indiv_ra:fltarr(ndup), $
           indiv_dec:dblarr(ndup), $
           indiv_sn:fltarr(ndup),$
           indiv_objname:strarr(ndup),$
           indiv_filename:strarr(ndup),$
           indiv_exptime:dblarr(ndup), $
           indiv_vhelio:dblarr(ndup), $
           indiv_good:bytarr(ndup), $
           indiv_count:-999L}

end

pro stackspec_ms0451
   common file, inputspecfile, massrange
   common npixcom, npix,ndup
   inputfile = '/scr2/nichal/workspace4/sps_fit/data/spline/sps_fit.fits.gz'
   massrange = [10.,10.5,11.,11.5]   
   n = obj_new('stackspec',inputfile=inputfile)

end
