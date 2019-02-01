pro compare_feonly_mgfe
;compare Fe only to MgFe 
   mgfe_ms = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
   feonly_ms = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit02.fits.gz',1)
   if total(mgfe_ms.ra-feonly_ms.ra) ne 0 then stop,'matching problem. Stopped.'
   good = where(mgfe_ms.good eq 1 and mgfe_ms.oiiew gt -5 and mgfe_ms.snfit gt 6.47 and $
                mgfe_ms.logmstar gt 6.,cgalms)
   mgfe_ms = mgfe_ms(good)
   feonly_ms = feonly_ms(good)
   print,'number of ms0451 galaxies:',cgalms
   ;calculate deltafeh using montecarlo, assuming two sided gaussian 
   ntime = 1000
   deltafeh_ms = fltarr(cgalms)
   deltafeherr_ms = fltarr(cgalms)

   for i=0,cgalms-1 do begin
      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_mgfe = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_mgfe(locpos) = gauss(locpos)*(mgfe_ms[i].fehupper-mgfe_ms[i].feh)+mgfe_ms[i].feh
      rndm_feh_mgfe(locneg) = gauss(locneg)*(mgfe_ms[i].feh-mgfe_ms[i].fehlower)+mgfe_ms[i].feh 

      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_feonly = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_feonly(locpos) = gauss(locpos)*(feonly_ms[i].fehupper-feonly_ms[i].feh)+feonly_ms[i].feh
      rndm_feh_feonly(locneg) = gauss(locneg)*(feonly_ms[i].feh-feonly_ms[i].fehlower)+feonly_ms[i].feh
   
      deltaarr = rndm_feh_feonly-rndm_feh_mgfe
      deltafeh_ms[i] =  median(deltaarr)
      deltafeherr_ms[i] = stdev(deltaarr)
   endfor 

   mgfe_cl = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit03.fits.gz',1)
   feonly_cl = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit02.fits.gz',1)
   if total(mgfe_cl.ra-feonly_cl.ra) ne 0 then stop,'matching problem. Stopped.'
   good = where(mgfe_cl.good eq 1 and mgfe_cl.oiiew gt -5 and mgfe_cl.snfit gt 6.83 and $
                mgfe_cl.logmstar gt 6.,cgalcl)
   mgfe_cl = mgfe_cl(good)
   feonly_cl = feonly_cl(good)
   print,'number of cl0024 galaxies:',cgalcl
  ;calculate deltafeh using montecarlo, assuming two sided gaussian 
   ntime = 1000
   deltafeh_cl = fltarr(cgalcl)
   deltafeherr_cl = fltarr(cgalcl)

   for i=0,cgalcl-1 do begin
      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_mgfe = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_mgfe(locpos) = gauss(locpos)*(mgfe_cl[i].fehupper-mgfe_cl[i].feh)+mgfe_cl[i].feh
      rndm_feh_mgfe(locneg) = gauss(locneg)*(mgfe_cl[i].feh-mgfe_cl[i].fehlower)+mgfe_cl[i].feh

      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_feonly = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_feonly(locpos) = gauss(locpos)*(feonly_cl[i].fehupper-feonly_cl[i].feh)+feonly_cl[i].feh
      rndm_feh_feonly(locneg) = gauss(locneg)*(feonly_cl[i].feh-feonly_cl[i].fehlower)+feonly_cl[i].feh

      deltaarr = rndm_feh_feonly-rndm_feh_mgfe
      deltafeh_cl[i] =  median(deltaarr)
      deltafeherr_cl[i] = stdev(deltaarr)
   endfor

   set_plot,'ps'
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()

   psname = 'compare_feonly_mgfe_mass_feh.eps'
   deltasym = '!9' + String("104B) + '!X'
   deltasym = '!9' + String("104B) + '!X'
   device, filename = psname,xsize = 30,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,2,1]
      xrange=[9.5,12]
      yrange=[-1.,0.3]
      plot,mgfe_ms.logmstar,mgfe_ms.feh,/nodata,xrange=xrange,xstyle=1,$
           yrange=yrange,ystyle=1,title='MS0451'
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

      mgfe_mssymsize = 1.5/(max(mgfe_ms.snfit)-min(mgfe_ms.snfit))*mgfe_ms.snfit+0.7
      feonly_mssymsize = 0.8/(max(feonly_ms.snfit)-min(feonly_ms.snfit))*feonly_ms.snfit+0.7
      for i=0,cgalms-1 do begin
        oplot,[mgfe_ms[i].logmstar],[mgfe_ms[i].feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=mgfe_mssymsize[i]
        oplot,[feonly_ms[i].logmstar],[feonly_ms[i].feh],psym=cgsymcat(16),color=fsc_Color('slateblue'),symsize=feonly_mssymsize[i]
      endfor
      al_legend,['mg fe','fe only'],/fill,psym=[46,16],colors=['tomato','slateblue'],box=0,/bottom,/right

      plot,mgfe_cl.logmstar,mgfe_cl.feh,/nodata,xrange=xrange,xstyle=1,$
           yrange=yrange,ystyle=1,title='CL0024'
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

      mgfe_clsymsize = 1.5/(max(mgfe_cl.snfit)-min(mgfe_cl.snfit))*mgfe_cl.snfit+0.7
      feonly_clsymsize = 0.8/(max(feonly_cl.snfit)-min(feonly_cl.snfit))*feonly_cl.snfit+0.7
      for i=0,cgalcl-1 do begin
        oplot,[mgfe_cl[i].logmstar],[mgfe_cl[i].feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=mgfe_clsymsize[i]
        oplot,[feonly_cl[i].logmstar],[feonly_cl[i].feh],psym=cgsymcat(16),color=fsc_Color('slateblue'),symsize=feonly_clsymsize[i]
      endfor
   device,/close

   !p.multi = [0,1,1]
   psname = 'compare_oldnew_deltafeh_deltaage.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plot,feonly_ms.feh-mgfe_ms.feh,alog10(feonly_ms.age/mgfe_ms.age),psym=cgsymcat(46),$
            xtitle=deltasym+'[Fe/H]',ytitle=deltasym+'log(Age)',/nodata
       for i=0,cgalms-1 do $
            oplot,[feonly_ms[i].feh-mgfe_ms[i].feh],[alog10(feonly_ms[i].age/mgfe_ms[i].age)],$
            psym=cgsymcat(46),color=fsc_color('crimson'),symsize=mgfe_mssymsize[i]
       for i=0,cgalcl-1 do $
            oplot,[feonly_cl[i].feh-mgfe_cl[i].feh],[alog10(feonly_cl[i].age/mgfe_cl[i].age)],$
            psym=cgsymcat(46),color=fsc_color('goldenrod'),symsize=mgfe_clsymsize[i]
   device,/close
   psname = 'compare_oldnew_deltafeh_mgfe.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plot,mgfe_ms.alphafe,feonly_ms.feh-mgfe_ms.feh,psym=cgsymcat(46),$
            xtitle='[Mg/Fe]',ytitle=deltasym+'[Fe/H]',/nodata
       for i=0,cgalms-1 do $
       oplot,[mgfe_ms[i].alphafe],[feonly_ms[i].feh-mgfe_ms[i].feh],psym=cgsymcat(46),$
           color=fsc_color('crimson'),symsize=mgfe_mssymsize[i]
       for i=0,cgalcl-1 do $
       oplot,[mgfe_cl[i].alphafe],[feonly_cl[i].feh-mgfe_cl[i].feh],psym=cgsymcat(46),$
           color=fsc_color('goldenrod'),symsize=mgfe_clsymsize[i]
   device,/close
end

pro compare_old_new_measurements
   compare_feonly_mgfe
end
