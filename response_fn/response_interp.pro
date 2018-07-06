function response_interp,element,abund,zh,age,silent=silent
common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
;Function to return response function of the specified element with model lambda
;It's very similar to add_response.f90 in alf. Some parts are from get_model.f90
;It's also an equivalent of sps_interp.pro
;INPUT
;element is a string array: choices are Fe,Mg,O,C,N,Si,Ca (Oxygen can't be negative)
;abund is a float array: same size as element - values of each element enhancement
;zh is the metallicity of base spectrum in log scale
;age is the age of base spectrum in gyr
;OUTPUT
;structure of {lambda:lambda,rspec:response_function}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;get input values
   logage = alog10(age)
   ;read structure 
   if (size(rsp_str))[1] eq 0 then begin
      rsp_str = mrdfits('/scr2/nichal/workspace4/response_fn/atlas_ssp_abund_krp.fits',1)
      logzgrid_rsp = rsp_str[uniq(rsp_str.zmet, sort(rsp_str.zmet))].zmet                 ;monotonic increasing
      agegrid_rsp = rsp_str[uniq(rsp_str.agegyr, sort(rsp_str.agegyr))].agegyr
      ;following setup.f90 in alf, don't use the models with logz=0.2
      wbad = where(logzgrid_rsp eq 0.2, cbad,complement=wgood)
      if cbad gt 0 then logzgrid_rsp=logzgrid_rsp(wgood)
   endif

   logagegrid_rsp = alog10(agegrid_rsp) ; "
   nage = n_elements(logagegrid_rsp)
   nzmet = n_Elements(logzgrid_rsp)

   ;set up interpolation
   vr = max([min([value_locate(logagegrid_rsp,logage),nage-2]),0])  ;interpolate age in log scale
   dr = (logage-logagegrid_rsp[vr])/(logagegrid_rsp[vr+1]-logagegrid_rsp[vr])
   dr = max([min([1.,dr]),0.])

   vm = max([min([value_locate(logzgrid_rsp,zh),nzmet-2]),0])
   dm = (zh-logzgrid_rsp)/(logzgrid_rsp[vm+1]-logzgrid_rsp[vm])
   dm = max([min([1.,dm]),0.])

   nlam = n_elements(rsp_str[0].lam)
   combr = fltarr(nlam)+1.
   validr = bytarr(n_elements(element))+1
   for ie = 0,n_Elements(element)-1 do begin
      range = 0.3
      ;get the plus and minus for specific elements
      case element[ie] of  
        'Mg': begin
               plusarr = rsp_str.mgp
               minusarr = rsp_str.mgm
              end
        'Fe': begin
               plusarr = rsp_str.fep
               minusarr = rsp_str.fem
              end
        'O': begin
              plusarr = rsp_str.asfep
              if abund[ie] lt 0. then begin
                  if ~keyword_set(silent) then print, 'Warning: oxygen cant be negative value. Skipped'
                  validr[ie] = 0
              endif else minusarr = rsp_str.asfep
              end
        'C': begin
              plusarr = rsp_str.cp
              minusarr = rsp_str.cm
              range = 0.15
             end
        'N': begin
              plusarr = rsp_str.np
              minusarr = rsp_str.nm
             end
        'Si': begin
               plusarr = rsp_str.sip
               minusarr = rsp_str.sim
              end
        'Ca': begin
               plusarr = rsp_str.cap
               minusarr = rsp_str.cam
              end
        'Ti': begin
               plusarr = rsp_str.tip
               minusarr = rsp_str.tim
              end
      else: begin
              print, 'element '+element[ie]+' not found, ignored'
              validr[ie] = 0
            end
      endcase

      if validr[ie] eq 1 then begin   
         ;get response fn at vr,vr+1,vm,vm+1 in to array of [nlam,2,2] size
         plus = fltarr(nlam,2,2)
         minus = fltarr(nlam,2,2)
         solar = fltarr(nlam,2,2)
      
         for ir=0,1 do for im=0,1 do begin
            loc = where(rsp_str.agegyr eq agegrid_rsp[vr+ir] and rsp_str.zmet eq logzgrid_rsp[vm+im],cloc)
            if cloc ne 1 then stop,'something is wrong with response grid'
            plus[*,ir,im] = plusarr[*,loc]
            minus[*,ir,im] = minusarr[*,loc]
            solar[*,ir,im] = rsp_str[loc].solar
         endfor

         if abund[ie] gt 0. then begin
            ;interpolate the response fn of element=+0.3 dex to the correct
            ;age and metallicity
            tmpr = dr*dm*plus[*,1,1]/solar[*,1,1] + $
                   (1.-dr)*dm*plus[*,0,1]/solar[*,0,1] + $
                   dr*(1.-dm)*plus[*,1,0]/solar[*,1,0] + $
                   (1.-dr)*(1.-dm)*plus[*,0,0]/solar[*,0,0]
            ;interpolate to the correct element abundace
            tmpr = (1+(tmpr-1.)*abund[ie]/range)
         endif else begin
            tmpr = dr*dm*minus[*,1,1]/solar[*,1,1] + $
                   (1.-dr)*dm*minus[*,0,1]/solar[*,0,1] + $
                   dr*(1.-dm)*minus[*,1,0]/solar[*,1,0] + $
                   (1.-dr)*(1.-dm)*minus[*,0,0]/solar[*,0,0]
            tmpr = (1.+(tmpr-1.)*ABS(abund[ie])/range)
         endelse
         combr = combr*tmpr
      endif
   endfor
   if total(validr) eq 0 then begin
      print,'No valid element for adding respons fn.'
      print, 'Return flat response function...'
   endif
   return, {lambda:rsp_str[loc].lam,rspec:combr}
end
