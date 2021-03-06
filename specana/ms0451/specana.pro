pro specana,retrieve=retrieve,redooii=redooii
   dir = ['combined_uniq_ms0451_polynomial','combined_uniq_ms0451_spline']
   headword = ['polynomial','spline']
   strname = dir+'.info.fits'
   oiifile = 'ms0451_oiiew.fits'
   if ~keyword_set(redooii) then oiistr = mrdfits(oiifile,1)
   if keyword_set(retrieve) then begin
      path = '/scr2/nichal/workspace4/prepspec/'
      for idir =0,1 do begin
         files = file_search(path+dir[idir]+'/*.fits', count=cspec)
         mass = fltarr(cspec)
         sn = fltarr(cspec) ;boxcar
         sngauss = fltarr(cspec)
         f814w = fltarr(cspec)
         redshift = fltarr(cspec)
         exptime = fltarr(cspec)   
         basename = strarr(cspec)
         oiiew = fltarr(cspec)
         oiiewerr =fltarr(cspec)
         for i=0,cspec-1 do begin
            spec = mrdfits(files[i],1)
            fits_open,files[i],fitsinfo
            cat = mrdfits(files[i],fitsinfo.nextend)
            free_lun,fitsinfo.unit
            mass[i] = cat.mass_kcorrect
            sn[i] = spec.sn
            f814w[i] = cat.f814w_auto
            redshift[i] = cat.z
            exptime[i] = spec.exptime
            basename[i] = file_basename(files[i])
            if keyword_Set(redooii) then begin 
               oiiew[i] = oii_ew(spec.contdiv,spec.lambda,spec.contdivivar,$
                              cat.z,widtherr=widtherr,znew=znew)
               oiiewerr[i] = widtherr
            endif else begin
               if strmatch(oiistr.basename[i],'*'+basename[i]+'*') ne 1 then stop,'ERROR:oii match'
               oiiew[i] = oiistr.oiiew[i]
               oiiewerr[i] = oiistr.oiiewerr[i]
            endelse
            if fitsinfo.nextend gt 4 then begin
               spec2 = mrdfits(files[i],3)
               sngauss[i] = spec2.sn
            endif
         endfor
         strsave = {mass:mass,sn:sn,sngauss:sngauss,f814w:f814w,redshift:redshift,$
                    exptime:exptime,basename:basename,dirname:dir[idir],$
                    headword:headword[idir],oiiew:oiiew,oiiewerr:oiiewerr}
         mwrfits,strsave,strname[idir],/silent,/create
      endfor
   endif
   
   for idir=0,1 do begin
      str = mrdfits(strname[idir],1)
      if idir eq 0 then strall = str else strall=[strall,str]
         set_plot,'ps'
      !p.multi = [0,2,1]
      !p.font = 0
      sunsym = sunsymbol()
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      psname=headword[idir]+'_SN_info.eps'
      device, filename = psname,xsize = 20,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
          plot,str.mass,str.sn,psym=cgsymcat(14),xtitle='log mass (M'+sunsym+')',ytitle='SN per pix',$
               xrange=[9,12],/nodata
          good = where(str.oiiew ge -5. and str.redshift gt 0.52 and str.redshift lt 0.56)
          bad = where(str.oiiew lt -5. and str.redshift gt 0.52 and str.redshift lt 0.56)
          oplot,str.mass(bad),str.sn(bad),psym=cgsymcat(15),color=fsc_color('blue')
          oplot,str.mass(good),str.sn(good),psym=cgsymcat(14),color=fsc_color('red')
          plot,str.sn,str.sngauss,xtitle='boxcar SN',ytitle='gauss SN',psym=cgsymcat(14)
          oplot,!x.crange,!x.crange,color=fsc_color('red')
      device,/close   
   endfor

   !p.multi = [0,1,1]
   psname=headword[0]+'_vs_'+headword[1]+'_SNcompare.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plot,strall[0].sn,strall[1].sn,psym=cgsymcat(14),xtitle='SN '+strall[0].headword,$
              ytitle='SN '+strall[1].headword
       oplot,!x.crange,!x.crange,color=fsc_color('red')
   device,/close
end
