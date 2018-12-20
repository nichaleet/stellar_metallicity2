pro compare_old_new_feonly
;1) compare old Fe only to new Fe only
   oldfe = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_cl0024_ana_afterrev_final.fits',1)
   newfe = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit02.fits.gz',1)
 
   ;match name
   match2,oldfe.objname,newfe.objname,subo,subn
   good = where(subo ne -1,cgal)
   oldfe = oldfe(good)
   newfe = newfe(subo(good))
   set_plot,'ps'
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()

   psname = 'plots/compare_old_new_mass_feh.eps'
   ;calculate deltafeh using montecarlo, assuming two sided gaussian 
   ntime = 1000
   deltafeh = fltarr(cgal)
   deltafeherr = fltarr(cgal)

   for i=0,cgal-1 do begin
      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_old = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_old(locpos) = gauss(locpos)*(oldfe[i].fehupper-oldfe[i].feh)+oldfe[i].feh
      rndm_feh_old(locneg) = gauss(locneg)*(oldfe[i].feh-oldfe[i].fehlower)+oldfe[i].feh 

      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_new = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_new(locpos) = gauss(locpos)*(newfe[i].fehupper-newfe[i].feh)+newfe[i].feh
      rndm_feh_new(locneg) = gauss(locneg)*(newfe[i].feh-newfe[i].fehlower)+newfe[i].feh
   
      deltaarr = rndm_feh_new-rndm_feh_old
      deltafeh[i] =  median(deltaarr)
      deltafeherr[i] = stdev(deltaarr)
   endfor 
   deltasym = '!9' + String("104B) + '!X'
   deltasym = '!9' + String("104B) + '!X'

   device, filename = psname,xsize = 15,ysize = 20, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,1,1]
      xrange=[9.5,12]
      yrange=[-1.,0.3]
      plot,oldfe.logmstar,oldfe.feh,/nodata,xrange=xrange,xstyle=5,$
           yrange=yrange,ystyle=5,position=[0.15,0.5,0.95,0.95]
      ;draw axis
;      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

;      cgerrplot,sci(goodfit).logmstar,sci(goodfit).fehlower,sci(goodfit).fehupper,color='tomato',thick=0.5
      oldsymsize = 1.5/(max(oldfe.snfit)-min(oldfe.snfit))*oldfe.snfit+0.7
      newsymsize = 0.8/(max(newfe.snfit)-min(newfe.snfit))*newfe.snfit+0.7
      for i=0,cgal-1 do begin
        oplot,[oldfe[i].logmstar],[oldfe[i].feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=oldsymsize[i]
        oplot,[newfe[i].logmstar],[newfe[i].feh],psym=cgsymcat(16),color=fsc_Color('limegreen'),symsize=newsymsize[i]
      endfor

      yrange=[-0.5,0.5]
      ploterror,newfe.logmstar,deltafeh,deltafeherr,psym=1,xrange=xrange,xstyle=5,$
           yrange=yrange,ystyle=5,position=[0.15,0.1,0.95,0.5],/noerase
      for i=0,cgal-1 do begin
        oplot,[newfe[i].logmstar],[deltafeh[i]],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=oldsymsize[i]
      endfor
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle=deltasym+'[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

   device,/close
end

pro compare_feonly_mgfe
;compare Fe only to MgFe 
   mgfe = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit03.fits.gz',1)
   feonly = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit02.fits.gz',1)
   if total(mgfe.ra-feonly.ra) ne 0 then stop,'matching problem. Stopped.'

   good = where(mgfe.good eq 1 and mgfe.oiiew gt -5 and mgfe.snfit gt 6.83 and $
                mgfe.logmstar gt 6.,cgal)
   mgfe = mgfe(good)
   feonly = feonly(good)
   cgal = n_Elements(mgfe)
   print,'number of galaxies:',cgal 

   ;calculate deltafeh using montecarlo, assuming two sided gaussian 
   ntime = 1000
   deltafeh = fltarr(cgal)
   deltafeherr = fltarr(cgal)

   for i=0,cgal-1 do begin
      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_mgfe = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_mgfe(locpos) = gauss(locpos)*(mgfe[i].fehupper-mgfe[i].feh)+mgfe[i].feh
      rndm_feh_mgfe(locneg) = gauss(locneg)*(mgfe[i].feh-mgfe[i].fehlower)+mgfe[i].feh 

      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_feonly = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_feonly(locpos) = gauss(locpos)*(feonly[i].fehupper-feonly[i].feh)+feonly[i].feh
      rndm_feh_feonly(locneg) = gauss(locneg)*(feonly[i].feh-feonly[i].fehlower)+feonly[i].feh
   
      deltaarr = rndm_feh_feonly-rndm_feh_mgfe
      deltafeh[i] =  median(deltaarr)
      deltafeherr[i] = stdev(deltaarr)
   endfor 

   set_plot,'ps'
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()

   psname = 'plots/compare_feonly_mgfe_mass_feh.eps'
   deltasym = '!9' + String("104B) + '!X'
   deltasym = '!9' + String("104B) + '!X'

   device, filename = psname,xsize = 15,ysize = 20, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,1,1]
      xrange=[9.5,12]
      yrange=[-1.,0.3]
      plot,mgfe.logmstar,mgfe.feh,/nodata,xrange=xrange,xstyle=5,$
           yrange=yrange,ystyle=5,position=[0.15,0.5,0.95,0.95]
      ;draw axis
;      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

;      cgerrplot,sci(goodfit).logmstar,sci(goodfit).fehlower,sci(goodfit).fehupper,color='tomato',thick=0.5
      mgfesymsize = 1.5/(max(mgfe.snfit)-min(mgfe.snfit))*mgfe.snfit+0.7
      feonlysymsize = 0.8/(max(feonly.snfit)-min(feonly.snfit))*feonly.snfit+0.7
      for i=0,cgal-1 do begin
        oplot,[mgfe[i].logmstar],[mgfe[i].feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=mgfesymsize[i]
        oplot,[feonly[i].logmstar],[feonly[i].feh],psym=cgsymcat(16),color=fsc_Color('limegreen'),symsize=feonlysymsize[i]
      endfor
      al_legend,['mg fe','fe only'],/fill,psym=[46,16],colors=['tomato','limegreen'],box=0,/bottom,/right

      yrange=[-0.5,0.5]
      ploterror,feonly.logmstar,deltafeh,deltafeherr,psym=1,xrange=xrange,xstyle=5,$
           yrange=yrange,ystyle=5,position=[0.15,0.1,0.95,0.5],/noerase
      for i=0,cgal-1 do begin
        oplot,[feonly[i].logmstar],[deltafeh[i]],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=mgfesymsize[i]
      endfor
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle=deltasym+'[Fe/H] (Feonly-MgFe)'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   device,/close

   psname = 'plots/compare_feonly_mgfe_delta_hist.eps'
   ;since it's a histogram, we don't need uncertainty. 
   ;we'll just simple difference (not montecarlo)
   !p.multi = [0,2,2]
   device, filename = psname,xsize = 25,ysize = 20, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plothist,feonly.feh-mgfe.feh,ytitle=deltasym+'[Fe/H] (Feonly-MgFe)',/halfbin
       plothist,alog10(feonly.age/mgfe.age),ytitle=deltasym+'log(Age) (Feonly-MgFe)',/halfbin
       plot,feonly.feh-mgfe.feh,alog10(feonly.age/mgfe.age),psym=cgsymcat(46),$
            xtitle=deltasym+'[Fe/H]',ytitle=deltasym+'log(Age)'
       plot,mgfe.alphafe,feonly.feh-mgfe.feh,psym=cgsymcat(46),$
            xtitle='[Mg/Fe]',ytitle=deltasym+'[Fe/H]'
   device,/close
end

pro compare_fiducial_testmodel,nametest,testfile
;compare testmodel to MgFe fiducial model
   namefiducial = 'MgNFe'
   fiducial = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit03.fits.gz',1)
   testmodel = mrdfits(testfile,1)

   if total(fiducial.ra-testmodel.ra) ne 0 then stop,'matching problem. Stopped.'

   good = where(fiducial.good eq 1 and fiducial.oiiew gt -5 and fiducial.snfit gt 6.83 and $
                fiducial.logmstar gt 6.,cgal)
   fiducial = fiducial(good)
   testmodel = testmodel(good)
   cgal = n_Elements(fiducial)
   print,'number of galaxies:',cgal 

   ;calculate deltafeh using montecarlo, assuming two sided gaussian 
   ntime = 1000
   deltafeh = fltarr(cgal)
   deltafeherr = fltarr(cgal)

   for i=0,cgal-1 do begin
      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_fiducial = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_fiducial(locpos) = gauss(locpos)*(fiducial[i].fehupper-fiducial[i].feh)+fiducial[i].feh
      rndm_feh_fiducial(locneg) = gauss(locneg)*(fiducial[i].feh-fiducial[i].fehlower)+fiducial[i].feh 

      gauss = randomn(seed,ntime) ;lower and upper
      rndm_feh_testmodel = fltarr(ntime)
      locpos = where(gauss gt 0.,complement=locneg)
      rndm_feh_testmodel(locpos) = gauss(locpos)*(testmodel[i].fehupper-testmodel[i].feh)+testmodel[i].feh
      rndm_feh_testmodel(locneg) = gauss(locneg)*(testmodel[i].feh-testmodel[i].fehlower)+testmodel[i].feh
   
      deltaarr = rndm_feh_testmodel-rndm_feh_fiducial
      deltafeh[i] =  median(deltaarr)
      deltafeherr[i] = stdev(deltaarr)
   endfor 

   set_plot,'ps'
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()

   psname = 'plots/compare_'+namefiducial+'_'+nametest+'_mass_feh.eps'
   deltasym = '!9' + String("104B) + '!X'
   deltasym = '!9' + String("104B) + '!X'

   device, filename = psname,xsize = 15,ysize = 20, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,1,1]
      xrange=[9.5,12]
      yrange=[-1.,0.3]
      plot,fiducial.logmstar,fiducial.feh,/nodata,xrange=xrange,xstyle=5,$
           yrange=yrange,ystyle=5,position=[0.15,0.5,0.95,0.95]
      ;draw axis
;      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

;      cgerrplot,sci(goodfit).logmstar,sci(goodfit).fehlower,sci(goodfit).fehupper,color='tomato',thick=0.5
      fiducialsymsize = 1.5/(max(fiducial.snfit)-min(fiducial.snfit))*fiducial.snfit+0.7
      testmodelsymsize = 0.8/(max(testmodel.snfit)-min(testmodel.snfit))*testmodel.snfit+0.7
      for i=0,cgal-1 do begin
        oplot,[fiducial[i].logmstar],[fiducial[i].feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=fiducialsymsize[i]
        oplot,[testmodel[i].logmstar],[testmodel[i].feh],psym=cgsymcat(16),color=fsc_Color('limegreen'),symsize=testmodelsymsize[i]
      endfor
      al_legend,[namefiducial,nametest],/fill,psym=[46,16],colors=['tomato','limegreen'],box=0,/bottom,/right

      yrange=[-0.5,0.5]
      ploterror,testmodel.logmstar,deltafeh,deltafeherr,psym=1,xrange=xrange,xstyle=5,$
           yrange=yrange,ystyle=5,position=[0.15,0.1,0.95,0.5],/noerase
      for i=0,cgal-1 do begin
        oplot,[testmodel[i].logmstar],[deltafeh[i]],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=fiducialsymsize[i]
      endfor
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle=deltasym+'[Fe/H] ('+nametest+'-'+namefiducial+')'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   device,/close

   psname = 'plots/compare_'+nametest+'_'+namefiducial+'_delta_hist.eps'
   ;since it's a histogram, we don't need uncertainty. 
   ;we'll just simple difference (not montecarlo)
   !p.multi = [0,3,2]
   device, filename = psname,xsize = 35,ysize = 20, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plothist,testmodel.feh-fiducial.feh,ytitle=deltasym+'[Fe/H] ('+nametest+'-'+namefiducial+')',/halfbin
       plothist,alog10(testmodel.age/fiducial.age),ytitle=deltasym+'log(Age) ('+nametest+'-'+namefiducial+')',/halfbin
       plothist,testmodel.alphafe-fiducial.alphafe,ytitle=deltasym+'[Mg/Fe] ('+nametest+'-'+namefiducial+')',/halfbin
       plot,testmodel.feh-fiducial.feh,alog10(testmodel.age/fiducial.age),psym=cgsymcat(46),$
            xtitle=deltasym+'[Fe/H]',ytitle=deltasym+'log(Age)'
       plot,testmodel.feh-fiducial.feh,testmodel.alphafe-fiducial.alphafe,psym=cgsymcat(46),$
            xtitle=deltasym+'[Fe/H]',ytitle=deltasym+'[Mg/Fe]'
       plot,testmodel.alphafe-fiducial.alphafe,alog10(testmodel.age/fiducial.age),psym=cgsymcat(46),$
            xtitle=deltasym+'[Mg/Fe]',ytitle=deltasym+'log(Age)'
   device,/close
stop
end

pro compare_old_new_measurements
   compare_feonly_mgfe
   compare_old_new_feonly
   compare_fiducial_testmodel,'Mgonly','/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit01.fits.gz'
end
