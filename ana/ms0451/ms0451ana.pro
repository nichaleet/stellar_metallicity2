pro ms0451ana
;GETTING SCI DATA
   sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline/sps_fit.fits.gz',1)
   ;sample selection
   good = where(sci.oiiew gt -5.,cgood)
   sci = sci(good)
   ageform = (galage(sci.zfit,1000.)/1.e9-sci.age)>0. ;age of universe when it was formed

   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   sunsym = sunsymbol()

   goodfit = where(sci.goodfit eq 1, ngoodfit, complement=badfit,ncomplement=nbadfit)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;S/N histograms
   sn = sci.snfit
   wnofit = where(sci.snfit eq 0,cnofit)
   if cnofit gt 0 then sn(wnofit) = sci(wnofit).sn
   ;This sn is per pixel. With 600 grating, the resolution is 0.65 A per pixel.
   ;To convert to per angstrom divide by sqrt(0.65) = x1.24
   sn = 1.24*sn
   psname='mass_sn.eps'
   device, filename = psname,xsize = 15,ysize = 15, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[9.5,12]
      binsize=0.25
      plot,sci.logmstar,sn,xrange=xrange,/nodata,xstyle=5,ytitle='S/N per Angstrom',$
          position=[0.1,0.45,0.95,0.9],/normal
      axis,xaxis=1,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      cgplot,sci(badfit).logmstar,sn(badfit),psym=16,color=fsc_color('cadetblue'),/overplot
      cgplot,sci(goodfit).logmstar,sn(goodfit),psym=46,color=fsc_color('tomato'),/overplot

      plothist,sci.logmstar,bin=binsize,xrange=xrange,xstyle=5,ystyle=5,/fill,fcolor=fsc_color('gray'), $
           /noerase,position=[0.1,0.1,0.95,0.45],/normal
      plothist,sci(goodfit).logmstar,bin=binsize,/overplot,color=fsc_color('org8'),/fill,fcolor=fsc_Color('tomato')
      plothist,sci(badfit).logmstar,bin=binsize,/overplot,color=fsc_color('seagreen'),/fline,fcolor=fsc_Color('cadetblue'),forientation=45
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,xaxis=1,xrange=xrange,xtickformat='(A1)',xstyle=1
      axis,yaxis=0,ytitle='# quiescent galaxies'
      axis,yaxis=1,ytickformat='(A1)'
   device,/close

   psname='ms0451_feh_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   
   device,/close

end
