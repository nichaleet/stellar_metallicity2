function response_interp,pos,range,dr,dm,solar,plus,minus
;very similar to add_response.f90 in alf
   if pos gt 0. then begin
      ;interpolate the response fn of element=+0.3 dex to the correct
      ;age and metallicity
      tmpr = dr*dm*plus[*,1,1]/solar[*,1,1] + $
             (1.-dr)*dm*plus[*,0,1]/solar[*,0,1] + $
             dr*(1.-dm)*plus[*,1,0]/solar[*,1,0] + $
             (1.-dr)*(1.-dm)*plus[*,0,0]/solar[*,0,0]
      ;interpolate to the correct element abundace
      tmpr = (1+(tmpr-1.)*pos/range)
   endif else begin
          tmpr = dr*dm*minus[*,1,1]/solar(*,1,1] + $
             (1.-dr)*dm*minus[*,0,1]/solar[*,0,1] + $
             dr*(1.-dm)*minus[*,1,0]/solar[*,1,0] + $
             (1.-dr)*(1.-dm)*minus[*,0,0]/solar[*,0,0]

      tmpr = (1.+(tmpr-1.)*ABS(pos)/range)
   endelse
   return, tmpr
end

function add_response_fn,element,xin,yin,a
;Function to return response function of the specified element
;element is string array: choices are Fe,Mg,O,C,N,Si,Ca (Oxygen can't be negative)
;a is array of 5 elements: [Zmet,age,vdisp, redshift,abund]
;xin is observed wavelength
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;get input values
   age = a[1]
   logage = alog10(a[1])
   zh = a[0]
   vdisp = a[2]
   redshift = a[3]
   abund = a[4]
   ;read structure 
   rstr = mrdfits('/scr2/nichal/workspace4/response_fn/atlas_ssp_abund_krp.fits',1)
   logzgrid = rstr[uniq(rstr.zmet, sort(rstr.zmet))].zmet                 ;monotonic increasing
   agegrid = rstr[uniq(rstr.agegyr, sort(rstr.agegyr))].agegyr
   logagegrid = alog10(agegrid) ; "
   ;following setup.f90 in alf, don't use the models with logz=0.2
   wbad = where(logzgrid eq 0.2, cbad,complement=wgood)
   if cbad gt 0 then logzgrid=logzgrid(wgood)
   nage = n_elements(logagegrid)
   nzmet = n_Elements(logzgrid)

   ;set up interpolation
   vr = max([min([value_locate(logagegrid,logage),nage-2]),0])  ;interpolate age in log scale
   dr = (logage-logagegrid[vr])/(logagegrid[vr+1]-logagegrid[vr])
   dr = max([min([1.,dr]),0.])

   vm = max([min([value_locate(logzgrid,zh),nzmet-2]),0])
   dm = (zh-logzgrid)/(logzgrid[vm+1]-logzgrid[vm])
   dm = max([min([1.,dm]),0.])

   range = 0.3
   ;get the plus and minus for specific elements
   case element of 
     'Mg': begin
           plusarr = rstr.mgp
           minusarr = rstr.mgm
           end
     'Fe': begin
           plusarr = rstr.fep
           minusarr = rstr.fem
           end
     'O': begin
           plusarr = rstr.asfep
           minusarr = (rstr.asfep*0.)+1.
           end
     'C': begin
           plusarr = rstr.cp
           minusarr = rstr.cm
           range = 0.15
           end
     'N': begin
           plusarr = rstr.np
           minusarr = rstr.nm
           end
     'Si': begin
           plusarr = rstr.sip
           minusarr = rstr.sim
           end
     'Ca': begin
           plusarr = rstr.cap
           minusarr = rstr.cam
           end
   else: print, 'element '+element+' not found'
   endcase

   ;get response fn at vr,vr+1,vm,vm+1 in to array of [nlam,2,2] size
   nlam = n_elements(rstr[0].lam)
   plus = fltarr(nlam,2,2)
   minus = fltarr(nlam,2,2)
   solar = fltarr(nlam,2,2)

   for ir=0,1 do for im=0,1 do begin
      loc = where(rstr.agegyr eq agegrid[vr+ir] and rstr.zmet eq logzgrid[vm+im],cloc)
      if cloc ne 1 then stop,'something is wrong with response grid'
      plus[*,ir,im] = plusarr[*,loc]
      minus[*,ir,im] = minuxarr[*,loc]
      solar[*,ir,im] = rstr[loc].solar
   endfor

   response_interp,abund,range,dr,dm,solar,plus,minus
 

end
