pro combspec::combine,spec,noredraw=nosdaredraw
   common npixcom, npix,npixcomb,ndup
   wgood = where(spec.indiv_good eq 1, cgood)
   if cgood eq 1 then begin
       goodpix = where(spec.indiv_lambda[*,wgood[0]] ne 0.,npixcomb)
       ipix = min(goodpix)
       fpix = max(goodpix)
       struct_replace_field,spec,'lambda',dblarr(npixcomb)
       struct_replace_field,spec,'dlam',dblarr(npixcomb)
       struct_replace_field,spec,'contdiv',dblarr(npixcomb)
       struct_replace_field,spec,'contdivivar',dblarr(npixcomb)
       spec.lambda = spec.indiv_lambda[ipix:fpix,wgood[0]]*(1.+spec.indiv_vhelio[wgood[0]]/3.e5)
       spec.dlam = spec.indiv_dlam[ipix:fpix,wgood[0]]
       spec.contdiv = spec.indiv_contdiv[ipix:fpix,wgood[0]]
       spec.contdivivar = spec.indiv_contdivivar[ipix:fpix,wgood[0]]
       spec.sn = spec.indiv_sn[wgood[0]]
   endif

   if cgood gt 1 then begin ;stack continuum and make output structure (strout)
       sn = spec.indiv_sn[wgood] ;per angstrom

       ;First, find the total wl range (lambdafull) 
       ;in the overlapped region use the pixel template of the spectrum with higher signal to noise
       ;Also, create a dlam array, which is the highest dlam at a given wavelength
       snorder = reverse(sort(sn)) ;order by high to low
       wgap=[]
       for k=0,cgood-1 do begin
          know = snorder[k] ;k-now
          lambnow = (1.+spec.indiv_vhelio[wgood[know]]/3e5)*spec.indiv_lambda[*,wgood[know]] ;shifted to helio centric
          goodpix = where(lambnow ne 0.,cgoodpix)
          lambnow = lambnow(goodpix)
          if k eq 0 then lambdafull = lambnow
          if k ge 1 then begin
            if lambnow[0] lt min(lambdafull) then begin
               waddpix = where(lambnow lt min(lambdafull),cwaddpix)
               lambdafull = [lambnow(waddpix),lambdafull]
               ;note if there is a gap and edit the existing gaps
               if n_elements(wgap) gt 0 then wgap = wgap+cwaddpix
               if cwaddpix eq cgoodpix then wgap = [wgap,cwaddpix-1]
            endif
            if lambnow[cgoodpix-1] gt max(lambdafull) then begin
               waddpix = where(lambnow gt max(lambdafull),cwaddpix)
               lambdafull = [lambdafull,lambnow(waddpix)]
               if cwaddpix eq cgoodpix then wgap = [wgap,n_elements(lambdafull)-cwaddpix-1]
            endif
            ;check the gap
            nwgap = n_elements(wgap)   
            if nwgap gt 0 then begin
               rmgap = []
               for gg=0,nwgap-1 do begin
                 waddpix = where(lambnow gt lambdafull[wgap[gg]] and lambnow lt lambdafull[wgap[gg]+1],cwaddpix)
                 if cwaddpix gt 0 then begin
                    nlambdafull = n_elements(lambdafull)
                    lambdafull = [lambdafull[0:wgap[gg]],lambnow(waddpix),lambdafull[wgap[gg]+1:nlambdafull-1]]
                    if max(waddpix) eq cgoodpix-1 then wgap = [wgap,wgap[gg]+cwaddpix]
                    if min(waddpix) eq 0 then wgap = [wgap,wgap[gg]] 
                    rmgap = [rmgap,gg]
                 endif
               endfor
               if n_elements(rmgap) gt 0 then begin
                  if n_Elements(rmgap) eq n_Elements(wgap) then wgap=[] else remove,rmgap,wgap
               endif
            endif
          endif
       endfor
       npixcomb = n_Elements(lambdafull)
       ;change the dimension of the output in the spec structure
       struct_replace_field,spec,'lambda',dblarr(npixcomb)
       struct_replace_field,spec,'dlam',dblarr(npixcomb)
       struct_replace_field,spec,'contdiv',dblarr(npixcomb)
       struct_replace_field,spec,'contdivivar',dblarr(npixcomb)
 
       spec.lambda = lambdafull

       ;get data in smooth and interpolate
       combmaskarr = bytarr(npixcomb,cgood)
       contdivarr = dblarr(npixcomb,cgood)
       contdivivararr = dblarr(npixcomb,cgood)
       dlamarr = dblarr(npixcomb,cgood)
       for k=0,cgood-1 do begin
          lambnow = (1.+spec.indiv_vhelio[wgood[k]]/3e5)*spec.indiv_lambda[*,wgood[k]] ;shifted to helio centric
          goodpix = where(lambnow ne 0.,cgoodpix)
          lambnow = lambnow(goodpix)
          specnow = spec.indiv_contdiv[goodpix,wgood[k]]
          ivarnow = spec.indiv_contdivivar[goodpix,wgood[k]]
          dlamnow = spec.indiv_dlam[goodpix,wgood[k]]
          masknow = spec.indiv_combmask[goodpix,wgood[k]]
          wonlamb = where(lambdafull ge min(lambnow) and lambdafull le max(lambnow),conlamb,$
                          complement=wofflamb,ncomplement=cofflamb)
          combmaskarr[*,k] = byte(interpol(masknow,lambnow,lambdafull))
          combmaskarr[wofflamb,k] = 0
          contdivarr[*,k] = interpol(specnow,lambnow,lambdafull)          
          contdivivararr[*,k] = interpol(ivarnow,lambnow,lambdafull)
          dlamarr[*,k] = interpol(dlamnow,lambnow,lambdafull)
          if cofflamb gt 0 then begin 
             combmaskarr[wofflamb,k] = 0
             contdivarr[wofflamb,k]=1./0.
             contdivivararr[wofflamb,k]=0
             dlamarr[wofflamb,k]=1./0.
          endif
       endfor
       wzero = where(contdivarr eq 0,czero)
       if czero gt 0 then combmaskarr(wzero) = 0

       ;weighted average
       ncombmaskarr = total(combmaskarr,2) ;size of ncombpix
       wnodup = where(ncombmaskarr eq 1, cwnodup)
       if cwnodup gt 0 then begin
           for ii=0,cwnodup-1 do begin
              wframe = where(combmaskarr[wnodup[ii],*] eq 1, cwframe)
              if cwframe ne 1 then stop,'oops something is wrong'
              spec.dlam[wnodup[ii]] = dlamarr[wnodup[ii],wframe]
              spec.contdiv[wnodup[ii]] = contdivarr[wnodup[ii],wframe]
              spec.contdivivar[wnodup[ii]] = contdivivararr[wnodup[ii],wframe]
           endfor
       endif
       wdup = where(ncombmaskarr gt 1, cwdup)
       for ii=0,cwdup-1 do begin
           wnan = where(~finite(contdivivararr[wdup[ii],*]) or (combmaskarr[wdup[ii],*] eq 0),cnan)
           dcontdiv = 1./sqrt(abs(contdivivararr[wdup[ii],*]))
           if cnan gt 0 then dcontdiv(wnan) = 1./0.
           err = abs(dcontdiv/contdivarr[wdup[ii],*])
           spec.dlam[wdup[ii]] = wmean(dlamarr[wdup[ii],*],abs(err*dlamarr[wdup[ii],*]),/nan)
           spec.contdiv[wdup[ii]]=wmean(contdivarr[wdup[ii],*],dcontdiv,error=contdiverr,/nan)
           spec.contdivivar[wdup[ii]]=1./contdiverr^2
       endfor
       ;signal to noise
       ;create mask
       contmask = bytarr(npixcomb)+1
       lambdarest = lambdafull/(1.+spec.z)
       for i=0,n_elements(*self.linestart)-1 do begin
          w = where(lambdarest ge (*self.linestart)[i] and lambdarest le (*self.lineend)[i], c)
          if c gt 0 then contmask[w] = 0
       endfor
       w = where(~finite(spec.contdiv) or ~finite(spec.contdivivar), c)
       if c gt 0 then contmask[w]=0
       ;calculate the deviation
       wcont = where(contmask eq 1)
       dev = abs(spec.contdiv[wcont] - 1.)
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then spec.sn = 1.0/mean(dev[w])/sqrt(median(spec.dlam))
   endif
   self->statusbox, spec=spec

   if ~keyword_set(noredraw) then begin
      self->redraw
   endif

end

pro combspec::combine_all
    specinfo = *self.specinfo
    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        self->default_range
        self->readspecindiv
        spec = *self.spec
        self->combine, spec
        ptr_free, self.spec
        self.spec = ptr_new(spec)
        self->writespecindiv,spec
    endfor
    self.i = curi
    self->readspecindiv
    spec = *self.spec
    self->statusbox, spec=spec
    self->redraw
end



; =================
pro combspec_event, ev
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


pro combspec::handle_button, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    case (uvalue) of
	'back': self->newspec, increment=-1
	'next': self->newspec, increment=1
	'default_range': begin
	    self->default_range
	    self->redraw
	end
	'fitcont': begin
	   spec = *self.spec
	   self->indivcontinuum,spec
           self->statusbox, spec=spec
	   ptr_free, self.spec
	   self.spec = ptr_new(spec)
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
	'combine_all': self->combine_all
	'backward': self->step, -1
	'forward': self->step, 1
	'save': begin
	    spec = *self.spec
	    self->writespecinfo
	    self->writespecindiv,spec
	end
	'savespec':begin
	    spec = *self.spec
	    self->writespecindiv,spec
	end
	'exit': begin
	    spec = *self.spec
	    self->writespecindiv,spec
	    self->writespecinfo
	    widget_control, ev.top, /destroy
	end
	else:
    endcase
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


pro combspec::step, increment ;for backward and forward button
    nmodes = 2
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
	mode eq nmodes-1 and increment gt 0: mode = 0
	mode eq 0 and increment lt 0: mode = nmodes-1
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


pro combspec::toggle_good
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
    spec = *self.spec
    widget_control, widget_info(self.base, find_by_uname='include'), get_value=good
    spec.indiv_good = good
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
    
end

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
    self->readspecindiv
    self->statusbox    
    self.ylim = minmax((*self.spec).indiv_spec,/nan)

    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw
    widget_control, widget_info(self.base, find_by_uname='zoom'), get_value=index
    wset, index
    spec = *self.spec
    z = spec.z
    plot,spec.lambda/(1.+z),spec.contdiv,xrange=[3700,3760],title='OII',position=[0.05,0.05,0.32,0.95],/normal, background=fsc_color('white'), color=fsc_color('black')
    for iplot=spec.indiv_count-1,0,-1 do oplot,spec.indiv_lambda[*,iplot]/(1.+z)*(1.+spec.indiv_vhelio[iplot]/3e5),spec.indiv_contdiv[*,iplot],color=fsc_color((*self.indivcolors)[iplot])
    oplot,spec.lambda/(1.+z),spec.contdiv,color=fsc_color('black'),thick=2

    plot,spec.lambda/(1.+z),spec.contdiv,xrange=[3900,4000],title='CaII',position=[0.37,0.05,0.64,0.95],/noerase,/normal, background=fsc_color('white'), color=fsc_color('black')
    for iplot=spec.indiv_count-1,0,-1 do oplot,spec.indiv_lambda[*,iplot]/(1.+z)*(1.+spec.indiv_vhelio[iplot]/3e5),spec.indiv_contdiv[*,iplot],color=fsc_color((*self.indivcolors)[iplot])
    oplot,spec.lambda/(1.+z),spec.contdiv,color=fsc_color('black'),thick=2

    plot,spec.lambda/(1.+z),spec.contdiv,xrange=[4830,4890],title='Hb',position=[0.69,0.05,0.96,0.95],/noerase,/normal, background=fsc_color('white'), color=fsc_color('black')
    for iplot=spec.indiv_count-1,0,-1 do oplot,spec.indiv_lambda[*,iplot]/(1.+z)*(1.+spec.indiv_vhelio[iplot]/3e5),spec.indiv_contdiv[*,iplot],color=fsc_color((*self.indivcolors)[iplot])
    oplot,spec.lambda/(1.+z),spec.contdiv,color=fsc_color('black'),thick=2


end


pro combspec::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.spec).indiv_spec,/nan)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 5]
        self.lambdalim = (minmax((*self.spec).lambda / (1d + (*self.spec).z)) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d +(*self.spec).z)
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
    zl = mode lt 1 ? (*self.spec).z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
end

; ============== TEXT BOXES =============
pro combspec::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
        'ylow': self->ylow, val
        'yhigh': self->yhigh, val
        'lambdalow': self->lambdalow, val
        'lambdahigh': self->lambdahigh, val
        'vhelio':self->editvhelio, val
        else:
    end
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro combspec::editvhelio,vhelio,noredraw=noredraw
   widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
   spec = *self.spec
   spec.indiv_vhelio[wplotcon]=vhelio
   ptr_free, self.spec
   self.spec = ptr_new(spec)
   self->statusbox, spec=spec
   if ~keyword_set(noredraw) then self->redraw
end

pro combspec::lambdalow, lambdalow, noredraw=noredraw
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->redraw
end


pro combspec::lambdahigh, lambdahigh, noredraw=noredraw
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->redraw
end


pro combspec::ylow, ylow, noredraw=noredraw
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 0: self.ylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->redraw
end


pro combspec::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 0: self.ylim[1] = yhigh
        else: self.divylim[1] = yhigh
    endcase
    self->redraw
end


; ============== COMBOBOX  =============
pro combspec::handle_combobox, ev
    self.i = ev.index
    self->statusbox
    self.ylim = minmax(((*self.spec).indiv_contdiv)*((*self.spec).indiv_continuum),/nan)
    self.keystate = 0
    self->redraw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end



; ============= DRAW CLICK =============
pro combspec::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 1 then click_coords /= 1d + (*self.spec).z
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
            zl = mode lt 2 ? (*self.spec).z : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
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
            zl = mode lt 1 ? (*self.spec).z : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin  ;control and click
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 0: ylim = self.ylim
                mode eq 1: ylim = self.divylim
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
pro combspec::handle_draw_key, ev
   common npixcom, npix,npixcomb,ndup 
   if ev.press ne 1 then return
   key = string(ev.ch)
   coords = convert_coord(ev.x, ev.y, /device, /to_data)
   widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
   spec= *self.spec
   widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
   case mode of
      0: begin        ;continuum fit
          case key of
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
              's': begin
                  case self.keystate of
                      0: begin
                          self->writespecindiv, spec
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
                              spec.indiv_contmask[w,wplotcon] = 0
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
                              spec.indiv_contmask[w,wplotcon] = 1
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

      1: begin        ;pixel mask
          case key of
              '0': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[~spec.indiv_good[0],spec.indiv_good[1:ndup-1]]
                          self->toggle_good
                      end
                      else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                  endcase
              end

              '1': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[spec.indiv_good[0],~spec.indiv_good[1],spec.indiv_good[2:ndup-1]]
                          self->toggle_good
                      end
                      else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                  endcase
              end

              '2': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[spec.indiv_good[0:1],~spec.indiv_good[2],spec.indiv_good[3:ndup-1]]
                          self->toggle_good
                      end
                      else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                  endcase
              end

              '3': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[spec.indiv_good[0:2],~spec.indiv_good[3],spec.indiv_good[4:ndup-1]]
                          self->toggle_good
                      end
                      else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                  endcase
              end

              '4': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[spec.indiv_good[0:3],~spec.indiv_good[4],spec.indiv_good[5:ndup-1]]
                          self->toggle_good
                      end
                      else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                  endcase
              end

              '5': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[spec.indiv_good[0:4],~spec.indiv_good[5],spec.indiv_good[6]]
                          self->toggle_good
                      end
                      else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                  endcase
              end

              '6': begin
                  case self.keystate of
                      0: begin
                          widget_control, widget_info(self.base, find_by_uname='good'), set_value=[spec.indiv_good[0:5],~spec.indiv_good[6]]
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
              's': begin
                  case self.keystate of
                      0: begin
                          self->writespecindiv, spec
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
                          multfactor = 1./(1.+spec.z)*(1.+spec.indiv_vhelio[wplotcon]/3e5)
                          w = where(spec.indiv_lambda[*,wplotcon]*multfactor gt lambda1 and spec.indiv_lambda[*,wplotcon]*multfactor lt lambda2, c)
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
                          multfactor = 1./(1.+spec.z)*(1.+spec.indiv_vhelio[wplotcon]/3e5)
                          w = where(spec.indiv_lambda[*,wplotcon]*multfactor gt lambda1 and spec.indiv_lambda[*,wplotcon]*multfactor lt lambda2, c)
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
              else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
          endcase
      end
   endcase
end


; ============== DESTROY ===============
pro combspec_cleanup, ev
    widget_control, ev, get_uvalue=obj
    obj_destroy, obj
end


pro combspec::cleanup
    ptr_free, self.spec
end

pro combspec::default_mask
    spec = *self.spec
    self->indivmask, spec
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    self->redraw
end

pro combspec::default_vhelio
    widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
    spec = *self.spec
    spec.indiv_vhelio[wplotcon] = helio_deimos(spec.indiv_ra[wplotcon],spec.indiv_dec[wplotcon],2000,jd=spec.indiv_jd[wplotcon])
    ptr_free, self.spec
    self.spec = ptr_new(spec)
    self->redraw
    self->statusbox
end

pro combspec::indivcontinuum, spec
   common npixcom, npix,npixcomb,ndup
   common continuum, conttype
   nobs = spec.indiv_count
   for i=0,nobs-1 do begin
         ipix=0
         fpix = spec.indiv_npix[i]-1
         contmask = spec.indiv_contmask[ipix:fpix,i]
         contmask[0:15]=0
         contmask[fpix-ipix-15:fpix-ipix]=0
         lambda = spec.indiv_lambda[ipix:fpix,i]
         ivar = spec.indiv_ivar[ipix:fpix,i]
         spectrum = spec.indiv_spec[ipix:fpix,i]
         won = where(contmask eq 1, complement=woff, con)
         if con gt 100 then begin
            case conttype of 
               'polycon': begin
                  degree = 7
                  coeff =  poly_fit(lambda[won],spectrum[won],degree,$
                                    measure_errors=1./sqrt(ivar[won]),status=status)
                  if status ne 0 then stop,'stopped:(continuumnormalize_deimos.pro) has some error'
                  cont = poly(lambda,reform(coeff))
                end
               'spline': begin
                  bkspace =  fix(400./median(abs(ts_diff(lambda,1)))) ;approximately every 400 A
                  bkpt = slatec_splinefit(lambda[won], spectrum[won], coeff, invvar=ivar[won], $
                                          bkspace=bkspace, upper=3, lower=3, /silent)
               
                  if bkpt[0] eq -1 then stop,'cannot do spline continnum fit'
                  cont = slatec_bvalu(lambda, bkpt, coeff)
                end
                else: stop,'no continuum type found'
            endcase
            check = where(finite(cont),ccheck)
            if ccheck lt 10 then stop
            contdiv = spec.indiv_spec[ipix:fpix,i]/cont
            contdivivar = spec.indiv_ivar[ipix:fpix,i]*cont^2
         endif else begin
            contdiv = dblarr(fpix-ipix+1)
            cont = dblarr(fpix-ipix+1)
            contdivivar = dblarr(fpix-ipix+1)
            spec.indiv_good[i] = 0  
         endelse
         spec.indiv_continuum[ipix:fpix,i] = cont
         spec.indiv_contdivivar[ipix:fpix,i] = contdivivar
         spec.indiv_contdiv[ipix:fpix,i] = contdiv
      wcont = where(spec.indiv_contmask[3:npix-4,i] eq 1)+3
      wcont = wcont[where(finite(spec.indiv_spec[wcont,i]) and finite(spec.indiv_continuum[wcont,i]) and spec.indiv_continuum[wcont,i] ne 0 and spec.indiv_lambda[wcont,i] gt 6700.)]
      dev = abs((spec.indiv_spec[wcont,i] - spec.indiv_continuum[wcont,i]) / spec.indiv_continuum[wcont,i])
      avgdev = mean(dev)
      w = where(dev lt 3.0*avgdev, c)
      if c gt 0 then spec.indiv_sn[i] = 1.0/mean(dev[w])/sqrt(median(spec.indiv_dlam[wcont,i]))
   endfor
end

pro combspec::indivmask,spec
   nobs = spec.indiv_count
   for i=0,nobs-1 do begin
      lambda = spec.indiv_lambda[*,i]
      lambdarest = lambda/(1.+spec.z)
      indivspec = spec.indiv_spec[*,i]
      ivar = spec.indiv_ivar[*,i]
      indiv_npix = spec.indiv_npix[i] 
      ;create continuum mask
      npix = n_elements(lambda)
      contmask = bytarr(npix)+1
      linestart = *self.linestart
      lineend = *self.lineend

      for j=0,n_elements(linestart)-1 do begin
         w = where(lambdarest ge linestart[j] and lambdarest le lineend[j], c)
         if c gt 0 then contmask[w] = 0
      endfor
   
      w = where(~finite(indivspec) or ~finite(lambda) or ~finite(ivar) or $
                ivar eq 0 or lambda eq 0 or indivspec eq 0, c)
      if c gt 0 then contmask[w]=0
   
      tellmask = bytarr(n_elements(lambda))
      tellstart = [6864., 7591., 8938.]
      tellend = [6935., 7694., 9500.]
      for j=0,n_elements(tellstart)-1 do begin
          w = where(lambda ge tellstart[j] and lambda le tellend[j], c)
          if c gt 0 then begin
              if j eq 0 then contmask[w] = 0
          endif
      endfor
      contmask[0:2] = 0
      contmask[npix-3:npix-1] = 0
      contmask[indiv_npix-3:indiv_npix-1]=0
      spec.indiv_contmask[*,i] = contmask
   endfor
end



; ============== REDRAW ===============
pro combspec::redraw
   common npixcom, npix,npixcomb,ndup
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Redrawing ...'
   self->default_range, /update
   widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
   widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
   wset, index
   spec = *self.spec
   indivcolors = *self.indivcolors
   indivbgcolors = *self.indivbgcolors
   widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
   case mode of
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       0: begin ;continuum fit
          if wplotcon lt spec.indiv_count then begin
             ipix = 0 ;min(where(spec.indiv_lambda[*,wplotcon] ne 0.))
             fpix = spec.indiv_npix[wplotcon]-1
             plot, spec.indiv_lambda[ipix:fpix,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5), spec.indiv_contdiv[ipix:fpix,wplotcon]*spec.indiv_continuum[ipix:fpix,wplotcon], xrange=self.lambdalim*(1d + spec.z), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
         
             ;;making the highlight region
             t = round(-1*ts_diff(spec.indiv_contmask[*,wplotcon], 1))
             wstart = where(t eq 1, cstart)+1
             wend = where(t eq -1, cend)
             if spec.indiv_contmask[0,wplotcon] eq 1 then begin
                if cstart eq 0 then begin
                   wstart = 0
                endif else begin
                   wstart = [0, wstart]
                  endelse
                cstart += 1
             endif
             if spec.indiv_contmask[n_elements(t)-1,wplotcon] eq 1 then begin
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
                x = ([spec.indiv_lambda[wstart[i],wplotcon], spec.indiv_lambda[wstart[i],wplotcon], spec.indiv_lambda[wend[i],wplotcon], $
                      spec.indiv_lambda[wend[i],wplotcon]] > (self.lambdalim[0] * (1d + spec.z))) < (self.lambdalim[1] * (1d + spec.z))
                y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
                polyfill, x, y, color=fsc_color(indivbgcolors[wplotcon],/brewer)
             endfor
             ;;finish highlighting masked regions

             ;;start ploting
             for idup=0,ndup-1 do begin 
                if spec.indiv_good[idup] eq 1 then begin
                   ipix = min(where(spec.indiv_lambda[*,idup] ne 0.))
                   fpix = spec.indiv_npix[idup]-1
                   oplot, spec.indiv_lambda[ipix:fpix,idup]*(1.+spec.indiv_vhelio[idup]/3e5), spec.indiv_contdiv[ipix:fpix,idup]*spec.indiv_continuum[ipix:fpix,idup], color=fsc_color(indivcolors(idup))
                endif
             endfor
             ipix = min(where(spec.indiv_lambda[*,wplotcon] ne 0.))
             fpix = spec.indiv_npix[wplotcon]-1
             oplot, spec.indiv_lambda[ipix:fpix,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5), spec.indiv_contdiv[ipix:fpix,wplotcon]*spec.indiv_continuum[ipix:fpix,wplotcon], color=fsc_color((*self.indivcolors)[wplotcon])
             oplot,spec.indiv_lambda[ipix:fpix,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5),spec.indiv_continuum[ipix:fpix,wplotcon], color=fsc_color('black')
             oplot,spec.indiv_lambda[ipix:fpix,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5),1./sqrt(spec.indiv_ivar[ipix:fpix,wplotcon]), color=fsc_color('gray') 
             ;;make the axes labels
             plot, spec.indiv_lambda[ipix:fpix,wplotcon]*(1.+spec.indiv_vhelio[wplotcon]/3e5), spec.indiv_contdiv[ipix:fpix,wplotcon]*spec.indiv_continuum[ipix:fpix,wplotcon], xrange=self.lambdalim * (1d + spec.z), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
             ;;highlight the telluric region and write line names
             n = n_elements(*self.tellstart)
             for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
             endfor
             n = n_elements(*self.linewaves)
             for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + spec.z) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + spec.z) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + spec.z), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + spec.z))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
             endfor
          endif
       end
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       1: begin                ;rest frame
          
          znow = spec.z
          plot, spec.lambda/(1d + znow), spec.contdiv, xrange=self.lambdalim, yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), /nodata
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
                multfactor = 1./(1d +znow)*(1.+spec.indiv_vhelio[wplotcon]/3e5)
                x = ([spec.indiv_lambda[wstart[i],wplotcon]*multfactor, spec.indiv_lambda[wstart[i],wplotcon]*multfactor, $
                      spec.indiv_lambda[wend[i],wplotcon]*multfactor, $
                      spec.indiv_lambda[wend[i],wplotcon]*multfactor] > (self.lambdalim[0])) < (self.lambdalim[1])
                y = [self.divylim[0], self.divylim[1], self.divylim[1], self.divylim[0]]
                polyfill, x, y, color=fsc_color(indivbgcolors[wplotcon],/brewer)
             endfor
             ;;finish highlighting masked regions
          endif
          ;plotting
          oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
          oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')

          for idup=0,ndup-1 do begin
             if spec.indiv_good[idup] eq 1 then begin
                ipix = min(where(spec.indiv_lambda[*,idup] ne 0.))
                fpix = spec.indiv_npix[idup]-1
                oplot,spec.indiv_lambda[ipix:fpix,idup]/(1d +znow)*(1.+spec.indiv_vhelio[idup]/3e5),spec.indiv_contdiv[ipix:fpix,idup],color=fsc_color(indivcolors(idup))
             endif
          endfor
         
          if total(spec.indiv_good) gt 1 then thick=1 else thick=1 
          oplot, spec.lambda/(1d +znow), spec.contdiv, color=fsc_color('black'),thick=thick

          ;;highlighting the telluric bands and label line names
          n = n_elements(*self.tellstart)
          for i=0,n-1 do begin
             oplot, [(*self.tellstart)[i], (*self.tellend)[i]] / (1d + znow), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
          endfor
          n = n_elements(*self.linewaves)
          for i=0,n-1 do begin
             if (*self.linewaves)[i] le !X.CRANGE[0] or (*self.linewaves)[i] ge !X.CRANGE[1] then continue
             oplot, [(*self.linewaves)[i], (*self.linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
             xyouts, (*self.linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
          endfor
        end
    endcase
   zl = mode lt 2 ? (*self.spec).z : 0.0
   widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
   widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro combspec::statusbox, spec=spec
    common npixcom, npix,npixcomb,ndup
    if ~keyword_set(spec) then spec = *self.spec
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='include'), set_value=spec.indiv_good
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(spec.objname, 2)+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(spec.z, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=strcompress(string(spec.sn, format='(D7.2)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=spec.LOGMSTAR_SED_BL gt 0 ? strcompress(string(spec.logmstar_sed_bl, format='(D10.2)'), /rem) : unknown

    widget_control, widget_info(self.base, find_by_uname='curoiiew'), set_value=strcompress(string(spec.EW_OII, format='(I8)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curf814w'), set_value=spec.f814wmag gt 0 ? strcompress(string(spec.f814wmag, format='(D8.2)'), /rem) : unknown


    inditable = replicate({instru:'',mask:'',vhelio:-999,sn:-999},ndup)
    inditable.instru = strcompress(spec.indiv_instru)
    inditable.mask = strcompress(spec.indiv_mask)
    inditable.vhelio = spec.indiv_vhelio
    inditable.sn =  spec.indiv_sn
    widget_control,widget_info(self.base,find_by_uname='curinditable'),set_value = inditable
 
    widget_control, widget_info(self.base, find_by_uname='iplotcont'), get_value=wplotcon
    widget_control, widget_info(self.base, find_by_uname='vhelio'), set_value=strcompress(string(spec.indiv_vhelio[wplotcon], format='(D7.2)'), /rem)
end


pro combspec::getspec, list=list
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common continuum, conttype
    common npixcom, npix,npixcomb,ndup

    npix = 10625 ;size of deimos old
    npixcomb = npix+1373 ;size of deimos old+mosfire
    observatory, 'keck', obs
    specfits = self.directory+'combspec_out.fits.gz'
    if ~file_test(specfits) then begin 
       nspec = n_elements(list)
       self.nspec = nspec
       specinfo = replicate({conttype:'',objnum:'',objname:''},nspec)
       conttypes = replicate(conttype,nspec)
       objnums = strarr(nspec)
       objnames = strarr(nspec)
       for i=0,nspec-1 do begin
          spec = {spec}
          obj = list[i]
          objnums[i] = strtrim(string(i),2)
          objnames[i] = strtrim(obj.phot_id,2)
          spec.conttype = conttypes[i]
          spec.objnum = objnums[i]
          spec.objname = objnames[i]
          specinfo[i].conttype = conttypes[i]
          specinfo[i].objnum = objnums[i]
          specinfo[i].objname = objnames[i]
          struct_assign,obj,spec,/nozero ;copy cat over
          spec.indiv_count = obj.ndup
          nobs = obj.ndup
          for iobs =0, nobs-1 do begin
              instru = obj.instrus[iobs]
              case strtrim(instru,2) of 
                 'deimos_old': begin
                    ;the structure will have 3 tags: spec, lambda, ivar of 10625 pixels
                    ;lambda is in rest-frame
                    str = mrdfits(strtrim(obj.files[iobs],2),1,hdr,/silent,status=status)
                    goodpix = where(str.ivar ne 0,cgoodpix)
                    ipix = min(goodpix)
                    fpix = max(goodpix)
                    indiv_npix = fpix-ipix+1
                    spec.indiv_npix[iobs] = indiv_npix
                    spec.indiv_spec[0:indiv_npix-1,iobs]=str.spec[ipix:fpix]
                    spec.indiv_lambda[0:indiv_npix-1,iobs]=str.lambda[ipix:fpix]*(obj.z+1.) ;change to observed lambda
                    spec.indiv_ivar[0:indiv_npix-1,iobs] = str.ivar[ipix:fpix]
                    spec.indiv_combmask[0:indiv_npix-1,iobs]=1
                    spec.indiv_filename[iobs] = obj.files[iobs]
                    spec.indiv_mask[iobs] = obj.masks[iobs]
                    spec.indiv_slit[iobs] = obj.slits[iobs]
                    spec.indiv_instru[iobs] = obj.instrus[iobs]
                    spec.indiv_ra[iobs] = obj.ra[iobs]
                    spec.indiv_dec[iobs] = obj.dec[iobs]
                    spec.indiv_dlam[0:indiv_npix-1,iobs] = 1.7 ;DEIMOS 1200 line grating, center at 7500, 1" slit
                    spec.indiv_good[iobs] = 1
                    end
                 'lris_old': begin
                    str = mrdfits(strtrim(obj.files[iobs],2),1,hdr,/silent,status=status)
                    indiv_npix = n_Elements(str.spec)
                    spec.indiv_npix[iobs] = indiv_npix
                    spec.indiv_spec[0:indiv_npix-1,iobs]=str.spec
                    spec.indiv_lambda[0:indiv_npix-1,iobs]=str.lambda*(obj.z+1.) ;change to observed lambda
                    spec.indiv_ivar[0:indiv_npix-1,iobs] = str.ivar
                    spec.indiv_combmask[0:indiv_npix-1,iobs]=1
                    spec.indiv_filename[iobs] = obj.files[iobs]
                    spec.indiv_mask[iobs] = obj.masks[iobs]
                    spec.indiv_slit[iobs] = obj.slits[iobs]
                    spec.indiv_instru[iobs] = obj.instrus[iobs]
                    spec.indiv_ra[iobs] = obj.ra[iobs]
                    spec.indiv_dec[iobs] = obj.dec[iobs]
                    spec.indiv_dlam[0:indiv_npix-1,iobs] = 9.18 ;LRIS 300 line grating, 1" slit
                    spec.indiv_good[iobs] = 1
                    end
                 'lris': begin
                    str = readfits(strtrim(obj.files[iobs],2),hdr)
                    indiv_npix = n_elements(str)
                    spec.indiv_spec[0:indiv_npix-1,iobs]=str
                    cdelt = sxpar(hdr,'cdelt1')
                    crval = sxpar(hdr,'crval1')
                    spec.indiv_lambda[0:indiv_npix-1,iobs]= dindgen(indiv_npix)*cdelt+crval
                    spec.indiv_npix[iobs] = indiv_npix
                    spec.indiv_combmask[0:indiv_npix-1,iobs]=1
                    spec.indiv_filename[iobs] = obj.files[iobs]
                    spec.indiv_instru[iobs] = obj.instrus[iobs]
                    radec = stringradec2deci(sxpar(hdr,'ra'),sxpar(hdr,'dec'))
                    spec.indiv_ra[iobs] = radec[0]
                    spec.indiv_dec[iobs] = radec[1]
                    spec.indiv_dlam[0:indiv_npix-1,iobs] = 6.9 ;LRIS400/8500 grating, 1" slit
                    spec.indiv_good[iobs] = 1
                    spec.indiv_npix[iobs] = indiv_npix
                    spec.indiv_jd = sxpar(hdr,'mjd-obs')+2400000.5 
                    basefile =  file_basename(obj.files[iobs])
                    extensions = strsplit(basefile,'.',/extract)
                    spec.indiv_mask[iobs] = extensions[0]
                    extensions[1]= 'rsig'
                    sigmafile = strtrim(file_dirname(obj.files[iobs])+'/'+strjoin(extensions,'.'),2)
                    if file_test(sigmafile) eq 0 then stop, 'cannot find file '+sigmafile
                    sigstr = readfits(sigmafile,hdr)
                    if sxpar(hdr,'crval1') ne crval or sxpar(hdr,'cdelt1') ne cdelt then $
                       stop,'check LRIS sigma and spec file:',obj.files[iobs]
                    spec.indiv_ivar[0:indiv_npix-1,iobs] = 1./sigstr^2
                    end
                 'mosfire': begin
                    str = mrdfits(strtrim(obj.files[iobs],2),1)
                    hdr = str.hdr
                    indiv_npix = n_Elements(str.lambda)
                    spec.indiv_npix[iobs] = indiv_npix
                    spec.indiv_spec[0:indiv_npix-1,iobs]=str.spec
                    spec.indiv_lambda[0:indiv_npix-1,iobs]=str.lambda
                    spec.indiv_ivar[0:indiv_npix-1,iobs] = 1./str.errspec^2
                    spec.indiv_combmask[0:indiv_npix-1,iobs]=1
                    spec.indiv_filename[iobs] = obj.files[iobs]
                    spec.indiv_instru[iobs] = obj.instrus[iobs]
                    spec.indiv_ra[iobs] = sxpar(hdr,'targra')
                    spec.indiv_dec[iobs]= sxpar(hdr,'targdec')
                    spec.indiv_dlam[0:indiv_npix-1,iobs] = 3.06 ;MOSFIRE Y BAND
                    spec.indiv_good[iobs] = 1
                    spec.indiv_npix[iobs] = indiv_npix
                    spec.indiv_jd = sxpar(hdr,'mjd-obs')+2400000.5 
                    basefile =  file_basename(obj.files[iobs])
                    extensions = strsplit(basefile,'.',/extract)
                    spec.indiv_mask[iobs] = extensions[1]
                    end
                 else: stop,'no instrument type found for '+obj.instrus[iobs]
              endcase
              wbadpix = where(spec.indiv_ivar[*,iobs] eq 0,cbadpix)
              if cbadpix gt 0 then spec.indiv_combmask[wbadpix,iobs] = 0
          endfor 
          self.i = i
          tell = mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent)
          wtell = n_elements(tell)-1
          tell = tell[wtell]
          ptr_free, self.tell
          self.tell = ptr_new(tell)
   
          self->indivmask,spec
          self->indivcontinuum,spec
          self->statusbox, spec=spec
          self->writespecindiv,spec
      endfor
      speclist = conttypes+' '+objnums+' '+objnames
      widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
      self.i = 0
      ptr_free, self.specinfo
      self.specinfo = ptr_new(specinfo)
      self.nspec = nspec
      self->writespecinfo
      speclist = conttypes+' '+strtrim(string(objnums), 2)+' '+objnames
      widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
      self->readspecindiv
      widget_control, widget_info(self.base, find_by_uname='mode'), set_value=1
   endif else begin
      specinfo = mrdfits(specfits, 1, /silent)
      ptr_free,self.specinfo
      self.specinfo = ptr_new(specinfo)

      self.nspec = n_elements(specinfo)
      speclist = specinfo.conttype+' '+strtrim(string(specinfo.objnum), 2)+' '+specinfo.objname
      widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

      tell = (mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent))
      wtell = n_elements(tell)-1
      tell = tell[wtell]
      ptr_free, self.tell
      self.tell = ptr_new(tell)

      self.i = 0
      self->readspecindiv

      widget_control, widget_info(self.base, find_by_uname='mode'), set_value=1
   endelse
end

pro combspec::writespecinfo
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing info to database ...'
    specinfo = *self.specinfo
    specfits = self.directory+'combspec_out.fits' 
    mwrfits, specinfo, specfits, /create, /silent
    spawn, 'gzip -f '+specfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro combspec::writespecindiv,spec
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing spectrum to database ...'
    specfits = self.directory+strcompress(spec.conttype)+'_'+spec.objnum+'_'+spec.objname+'.fits'
    mwrfits, spec, specfits, /create, /silent
    spawn, 'gzip -f '+specfits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end

pro combspec::readspecindiv
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading spectrum from database ...'
    i=self.i
    specinfo = *self.specinfo
    specfits = self.directory+strtrim(specinfo[i].conttype,2)+'_'+strtrim(specinfo[i].objnum,2)+'_'+strtrim(specinfo[i].objname,2)+'.fits.gz'
    if file_test(specfits) eq 0 then stop,'no file found'
    spec = mrdfits(specfits,1,/silent)
    ptr_free,self.spec
    self.spec = ptr_new(spec)    
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'

end

pro combspec::initialize_directory, directory=directory,rematch=rematch
   common continuum, conttype
   list = mrdfits('matched_cl1604obj_members.fits',1)

   newdirectory = '/scr2/nichal/workspace4/prepspec_cl1604/data/'+conttype+'/'
   if ~file_test(newdirectory) then file_mkdir, newdirectory

   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

   specfits = newdirectory+'combspec_out.fits.gz'
   self.directory = newdirectory

   self->getspec, list=list
   self.i = 0
   self->readspecindiv
   spec = *self.spec
   self->statusbox, spec=spec
   self->default_range
   self.lambdalim = [3300, 7000]
   self->newspec, increment=0
   widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

; =============== INIT ================
function combspec::INIT, directory=directory
    common npixcom, npix,npixcomb,ndup

    base = widget_base(/row, title='combine_spec_cl1604', uvalue=self, mbar=menu, tab_mode=0, units=1)
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wcombineall = widget_button(tools_menu, value = 'combine all spectra',uname='combine_all',uvalue='combine_all')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base
    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit', 'rest frame'], /column, /exclusive, set_value=1, uname='mode', uvalue='mode', /no_release)
    wprepbase = widget_base(wleft, /row, /align_center)
    wfitcont = widget_button(wprepbase, value='Fit Continuum', uvalue='fitcont', uname='fitcont', tab_mode=1, xsize=100)
    wcombine = widget_button(wprepbase, value='Re-stack', uvalue='combine', uname='combine', tab_mode=1, xsize=100)
    wsavespec = widget_button(wprepbase,value='Save spec', uvalue='savespec',uname='savespec',tab_mode=1,xsize=100)
    windivbase = widget_base(wleft,/row,/align_center)
    wdefaultmaskbase = widget_base(windivbase,/row,/align_center)
    wdefaultmask = widget_button(wdefaultmaskbase,value='Default Cont Mask',uvalue='default_mask',uname='default_mask',tab_mode=1,xsize=120)    
    wdefaultvheliobase = widget_base(windivbase,/row,/align_center)
    wdefaultvhelio = widget_button(wdefaultvheliobase,value='Reset Vhelio',uvalue='default_vhelio',uname='default_vhelio',tab_mode=1,xsize=120)
    wincludebase = widget_base(wleft, /row, /align_center,/frame)
    wincluderadio = cw_bgroup(wincludebase,['0','1','2','3'],/column,/nonexclusive,set_value=[1,1,1,1],uname='include',uvalue='include',label_top='include')
    wplotconradio = cw_bgroup(wincludebase,['0','1','2','3'],/column,/exclusive,set_value=0,uname='iplotcont',uvalue='iplotcont',label_top='plot')
    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='s/n = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)

    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='redshift = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)

    wmstarbase = widget_base(wcurobj, /align_center, /row)
    wmstarlabel = widget_label(wmstarbase, value='log M* = ', /align_right, uname='mstarlabel', xsize=95)
    wcurmstar = widget_label(wmstarbase, value='     ', /align_left, uname='curmstar', uvalue='curmstar', xsize=150)

    woiiewbase = widget_base(wcurobj, /align_center, /row)
    woiiewlabel = widget_label(woiiewbase, value='OII EW = ', /align_right, uname='oiiewlabel', xsize=95)
    wcuroiiew = widget_label(woiiewbase, value='     ', /align_left, uname='curoiiew', uvalue='curoiiew', xsize=150)

    wf814wbase = widget_base(wcurobj, /align_center, /row)
    wf814wlabel = widget_label(wf814wbase, value='F814W = ', /align_right, uname='f814wlabel', xsize=95)
    wcurf814w = widget_label(wf814wbase, value='     ', /align_left, uname='curf814w', uvalue='curf814w', xsize=150)


    windibase = widget_base(wcurobj,/align_left,/row,xsize=400)
    winditable = widget_table(windibase,value=replicate({instru:'',mask:'',vhelio:-999,sn:-999},5),/row_major,column_labels=['instru','mask','vhelio','sn'],uname='curinditable',uvalue='curinditable')

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
    xmanager, 'combspec', self.base, /no_block, cleanup='combspec_cleanup'

    readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
    wbands = [1,2,3,4,5]
    tellstart = tellstart[wbands]
    tellend = tellend[wbands]
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick
    self.linewaves = ptr_new([2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85,4384,4455,4531,5015,5270,5335,5406,5709,5782])
    self.linenames = ptr_new(['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', '[NeIII]', 'HeI', 'H8', 'CaH', '[NeIII]', 'CaK', 'He', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', 'Mgb', 'Mgb', 'Mgb', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]','Fe4384','Fe4455','Fe4531','Fe5015','Fe5270','Fe5335','Fe5406','Fe5709','Fe5782'])
    self.linecolors = ptr_new(['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue','red','red','red','red','red','red','red','red','red'])
    self.tellstart = ptr_new(tellstart)
    self.tellend = ptr_new(tellend)
    self.tellthick = ptr_new([5, 2, 5, 2, 2])

    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    self.linestart = ptr_new(linestart)
    self.lineend = ptr_new(lineend)
    self.linetype = ptr_new(linetype)
    self.indivcolors = ptr_new(['blueviolet','blue','seagreen','springgreen','gold','chocolate','firebrick'])
    self.indivbgcolors = ptr_new(['pbg1','blu1','ryb5','grn1','wt4','org1','red2'])
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
             specinfo:ptr_new(), $
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
             lambda1:0d, $
             indivcolors:ptr_new(),$
             indivbgcolors:ptr_new()}
end

pro spec__define
   common npixcom, npix,npixcomb,ndup
   spec = {spec, $
           conttype:'', $
           objnum:'', $
           objname:'', $
           lambda:dblarr(npixcomb), $
           dlam:dblarr(npixcomb), $
           contdiv:dblarr(npixcomb), $
           contdivivar:dblarr(npixcomb), $
           sn:-999d, $
           indiv_spec:dblarr(npix,ndup), $
           indiv_lambda:dblarr(npix,ndup), $
           indiv_ivar:dblarr(npix,ndup), $
           indiv_contmask:bytarr(npix,ndup), $
           indiv_combmask:bytarr(npix,ndup), $
           indiv_dlam:dblarr(npix,ndup), $
           indiv_continuum:dblarr(npix,ndup), $
           indiv_contdiv:dblarr(npix,ndup), $
           indiv_contdivivar:dblarr(npix,ndup), $
           indiv_sn:fltarr(ndup),$
           indiv_ra:dblarr(ndup), $
           indiv_dec:dblarr(ndup), $
           indiv_mask:strarr(ndup),$
           indiv_slit:strarr(ndup),$
           indiv_filename:strarr(ndup),$
           indiv_jd:dblarr(ndup), $
           indiv_vhelio:dblarr(ndup), $
           indiv_good:bytarr(ndup), $
           indiv_count:-999L,$
           indiv_instru:strarr(ndup),$
           indiv_npix:lonarr(ndup),$
           PHOT_ID:'',$      
           MASK:'',$         
           SLIT:'',$         
           RA:-999d,$        
           DEC:-999d,$       
           RA_LFC:-999d,$    
           DEC_LFC:-999d,$   
           RMAG:-999d,$      
           IMAG:-999d,$      
           ZMAG:-999d,$      
           Z:-999d,$         
           ZERR:-999d,$      
           ZQUALITY:-999,$   
           OLDID:'',$        
           PHOT_FLAG:'',$    
           RA_ACS:-999d,$    
           DEC_ACS:-999d,$   
           ACS_ID:'',$       
           F606WMAG:-999d,$  
           F814WMAG:-999d,$  
           COMMENT:'',$                          
           D4000:-999d,$     
           D4000ERR:-999d,$  
           KMAG_AB:-999d,$   
           M_K:-999d,$       
           M_KERR:-999d,$    
           L_K:-999d,$       
           L_KERR:-999d,$                                          
           MSTAR_K:-999d,$   
           MSTAR_KERR:-999d,$
           LOGMSTAR_K:-999d,$
           LOGMSTAR_KERR:-999d,$   
           MSTAR_SED_BL:-999d,$    
           MSTAR_SED_BLERR:-999d,$ 
           LOGMSTAR_SED_BL:-999d,$ 
           LOGMSTAR_SED_BLERR:-999d,$ 
           CONTEXT:-999,$         
           CONTEXTFLAG:-999,$     
           JMAG:-999d,$           
           JMAGERR:-999d,$        
           SFR_24M:-999d,$        
           EW_OIIERR:-999d,$      
           SFR_OII:-999d,$        
           EW_OII:-999d,$         
           LOGMSTAR_SED_PFW:-999d,$ 
           AGE_SED:-999d,$       
           LOGAGE_SED:-999d,$    
           CAH:-999d,$           
           CAHERR:-999d,$        
           GBAND:-999d,$         
           GBANDERR:-999d,$      
           MASS_KCORRECT:-999d}

end

pro combspec_cl1604
   common continuum, conttype
   common npixcom, npix,npixcomb,ndup
   ndup = 4
   conttype = 'spline'
   n = obj_new('combspec',directory=directory)
end
