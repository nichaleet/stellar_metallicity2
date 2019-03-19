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

      mgfe_mssymsize = 1.5/(36-6)*(mgfe_ms.snfit-0.)+0.7
      feonly_mssymsize = 0.8/(36-6)*(feonly_ms.snfit-0.)+0.7
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

      mgfe_clsymsize = 1.5/(36-6.)*(mgfe_cl.snfit-0.)+0.7
      feonly_clsymsize = 0.8/(36-6.)*(feonly_cl.snfit-0.)+0.7
      for i=0,cgalcl-1 do begin
        oplot,[mgfe_cl[i].logmstar],[mgfe_cl[i].feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=mgfe_clsymsize[i]
        oplot,[feonly_cl[i].logmstar],[feonly_cl[i].feh],psym=cgsymcat(16),color=fsc_Color('slateblue'),symsize=feonly_clsymsize[i]
      endfor
   device,/close

   !p.multi = [0,1,1]
   psname = 'compare_oldnew_deltafeh_deltaage.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[-0.4,0.4]
      yrange=[-0.3,0.2]
      plothist,[feonly_ms.feh-mgfe_ms.feh,feonly_cl.feh-mgfe_cl.feh],$
            /fill,xrange=xrange,xstyle=5,ystyle=5,/halfbin,$
            color=fsc_color('gray'),fcolor=fsc_color('gray'),yrange=[0,100]
      plothist,[alog10(feonly_ms.age/mgfe_ms.age),alog10(feonly_cl.age/mgfe_cl.age)],$
            /fill,yrange=yrange,ystyle=5,xstyle=5,/rotate,xrange=[0,120],/noerase,$
            fcolor=fsc_color('gray'),color=fsc_color('gray'),/halfbin
      plot,feonly_ms.feh-mgfe_ms.feh,alog10(feonly_ms.age/mgfe_ms.age),psym=cgsymcat(46),$
            xtitle=deltasym+'[Fe/H]',ytitle=deltasym+'log(Age)',/nodata,xrange=xrange,$
            xstyle=1,yrange=yrange,ystyle=1,/noerase
      for i=0,cgalms-1 do $
            oplot,[feonly_ms[i].feh-mgfe_ms[i].feh],[alog10(feonly_ms[i].age/mgfe_ms[i].age)],$
            psym=cgsymcat(46),color=fsc_color('crimson'),symsize=mgfe_mssymsize[i]
      for i=0,cgalcl-1 do $
            oplot,[feonly_cl[i].feh-mgfe_cl[i].feh],[alog10(feonly_cl[i].age/mgfe_cl[i].age)],$
            psym=cgsymcat(46),color=fsc_color('goldenrod'),symsize=mgfe_clsymsize[i]
      oplot,[-0.1,0.1,0.1,-0.1,-0.1],[-0.1,-0.1,0.1,0.1,-0.1],psym=0,linestyle=1       
   device,/close
   psname = 'compare_oldnew_deltafeh_mgfe.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       xrange=[-0.4,0.8]
       yrange=[-0.4,0.4]
       plot,mgfe_ms.alphafe,feonly_ms.feh-mgfe_ms.feh,psym=cgsymcat(46),$
            xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,/nodata

       x = [mgfe_ms.alphafe,mgfe_cl.alphafe]
       y = [feonly_ms.feh-mgfe_ms.feh,feonly_cl.feh-mgfe_cl.feh]
       yerr = [sqrt(0.25*(feonly_ms.fehupper-feonly_ms.fehlower)^2+0.25*(mgfe_ms.fehupper-mgfe_ms.fehlower)^2),sqrt(0.25*(feonly_cl.fehupper-feonly_cl.fehlower)^2+0.25*(mgfe_cl.fehupper-mgfe_cl.fehlower)^2)]
       linparall = linfit(x,y,measure_errors=yerr,sigma=sigmaall)

       x = [mgfe_cl.alphafe]
       y = [feonly_cl.feh-mgfe_cl.feh]
       yerr = [sqrt(0.25*(feonly_cl.fehupper-feonly_cl.fehlower)^2+0.25*(mgfe_cl.fehupper-mgfe_cl.fehlower)^2)]
       linparcl = linfit(x,y,measure_errors=yerr,sigma=sigmacl)

       xshade=findgen(51)/50.*(max(xrange)-min(xrange))+min(xrange)
       yshade=fltarr(51,4)
       signm = [+1.,+1.,-1.,-1.]
       signb = [+1.,-1.,+1.,-1.]
       for i=0,3 do begin
          b = linparcl[0]+signb[i]*sigmacl[0]
          m = linparcl[1]+signm[i]*sigmacl[1]
          yshade[*,i]=xshade*m+b
       endfor
       yshade_maxcl = max(yshade,dimension=2)
       yshade_mincl = min(yshade,dimension=2)

       xpoly = [xshade,reverse(xshade)]
       ypoly = [yshade_maxcl,reverse(yshade_mincl)]      
       polyfill,xpoly,ypoly,color=fsc_color('gray')

       for i=0,3 do begin
          b = linparall[0]+signb[i]*sigmaall[0]
          m = linparall[1]+signm[i]*sigmaall[1]
          yshade[*,i]=xshade*m+b
       endfor
       yshade_maxall = max(yshade,dimension=2)
       yshade_minall = min(yshade,dimension=2)

       xpoly = [xshade,reverse(xshade)]
       ypoly = [yshade_maxall,reverse(yshade_minall)]
;       polyfill,xpoly,ypoly,color=fsc_color('rose')
       for i=0,cgalms-1 do $
       oplot,[mgfe_ms[i].alphafe],[feonly_ms[i].feh-mgfe_ms[i].feh],psym=cgsymcat(46),$
           color=fsc_color('crimson'),symsize=mgfe_mssymsize[i]
       for i=0,cgalcl-1 do $
       oplot,[mgfe_cl[i].alphafe],[feonly_cl[i].feh-mgfe_cl[i].feh],psym=cgsymcat(46),$
           color=fsc_color('goldenrod'),symsize=mgfe_clsymsize[i]

;       oplot,!x.crange,linparall[0]+linparall[1]*!x.crange      
       oplot,!x.crange,linparcl[0]+linparcl[1]*!x.crange,linestyle=1 
       axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='[Mg/Fe]'
       axis,yaxis=0,yrange=yrange,ystyle=1,ytitle=deltasym+'[Fe/H]'
       axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
       axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   device,/close
   print,'all par:',linparall,sigmaall
   print,'Cl par:',linparcl,sigmacl

   psname = 'compare_oldnew_deltafeh_mgfe_bysn.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       xrange=[-0.4,0.8]
       yrange=[-0.4,0.4]
       plot,mgfe_ms.alphafe,feonly_ms.feh-mgfe_ms.feh,psym=cgsymcat(46),$
            xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,/nodata

       x = [mgfe_ms.alphafe,mgfe_cl.alphafe]
       y = [feonly_ms.feh-mgfe_ms.feh,feonly_cl.feh-mgfe_cl.feh]
       yerr = [sqrt(0.25*(feonly_ms.fehupper-feonly_ms.fehlower)^2+0.25*(mgfe_ms.fehupper-mgfe_ms.fehlower)^2),sqrt(0.25*(feonly_cl.fehupper-feonly_cl.fehlower)^2+0.25*(mgfe_cl.fehupper-mgfe_cl.fehlower)^2)]
       sn = [mgfe_ms.snfit,mgfe_cl.snfit]
       mediansn = median(sn)
       highsn = where(sn gt mediansn,complement=lowsn)
       linpar_highsn = linfit(x(highsn),y(highsn),measure_errors=yerr(highsn),sigma=sigma_highsn)
       linpar_lowsn = linfit(x(lowsn),y(lowsn),measure_errors=yerr(lowsn),sigma=sigma_lowsn)
       linpar_allsn = linfit(x,y,measure_errors=yerr,sigma=sigma_allsn)
       xshade=findgen(51)/50.*(max(xrange)-min(xrange))+min(xrange)
       yshade=fltarr(51,4)
       signm = [+1.,+1.,-1.,-1.]
       signb = [+1.,-1.,+1.,-1.]

       for i=0,3 do begin
          b = linpar_allsn[0]+signb[i]*sigma_allsn[0]
          m = linpar_allsn[1]+signm[i]*sigma_allsn[1]
          yshade[*,i]=xshade*m+b
       endfor
       yshade_max_allsn = max(yshade,dimension=2)
       yshade_min_allsn = min(yshade,dimension=2)

       xpoly = [xshade,reverse(xshade)]
       ypoly = [yshade_max_allsn,reverse(yshade_min_allsn)]
       polyfill,xpoly,ypoly,color=fsc_color('lightgray')

       for i=0,3 do begin
          b = linpar_highsn[0]+signb[i]*sigma_highsn[0]
          m = linpar_highsn[1]+signm[i]*sigma_highsn[1]
          yshade[*,i]=xshade*m+b
       endfor
       yshade_max_highsn = max(yshade,dimension=2)
       yshade_min_highsn = min(yshade,dimension=2)

       xpoly = [xshade,reverse(xshade)]
       ypoly = [yshade_max_highsn,reverse(yshade_min_highsn)]      
       polyfill,xpoly,ypoly,color=fsc_color('gray')

       for i=0,cgalms-1 do $
       oplot,[mgfe_ms[i].alphafe],[feonly_ms[i].feh-mgfe_ms[i].feh],psym=cgsymcat(46),$
           color=fsc_color('crimson'),symsize=mgfe_mssymsize[i]
       for i=0,cgalcl-1 do $
       oplot,[mgfe_cl[i].alphafe],[feonly_cl[i].feh-mgfe_cl[i].feh],psym=cgsymcat(46),$
           color=fsc_color('goldenrod'),symsize=mgfe_clsymsize[i]

;       oplot,!x.crange,linparall[0]+linparall[1]*!x.crange      
       oplot,!x.crange,linpar_highsn[0]+linpar_highsn[1]*!x.crange,linestyle=1 
       axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='[Mg/Fe]'
       axis,yaxis=0,yrange=yrange,ystyle=1,ytitle=deltasym+'[Fe/H]'
       axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
       axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   device,/close
stop
end

pro compare_old_new_measurements
   compare_feonly_mgfe
end
