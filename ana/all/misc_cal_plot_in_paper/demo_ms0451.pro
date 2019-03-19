pro demo_ms0451
   sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/arxiv/sps_fit03.fits.gz',1)
   sci.oiiew = sci.oiiew*(-1.)
   whiem = where(sci.oiiew gt 22.5)
   sci(whiem).oiiew = 23.
   ;galex fuv is 1540 A*1.55 = 2387
   ;galex nuv is 2315 A
   ;V band is 5510 A *1.55 8540
   ;use F814
   ;alternatively, use KCORRECT

   quiescent = where(sci.oiiew lt 5. and sci.nuv gt 50. and sci.logmstar gt 6.,cquiescent)
   good = where(sci.oiiew lt 5. and sci.nuv gt 50. and sci.snfit gt 6.47 and sci.logmstar gt 6.,cgood)
   set_plot,'ps'
   !p.multi = [0,2,1]
   !p.font = 0
   !p.charsize =1
   ;oii ew histogram
   psname='ms0451_quiescent_histogram.eps'
   device, filename = psname,xsize = 18,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.charsize=1.2

      plothist,sci.oiiew,bin=2.5,xrange=[-7.5,25],xtitle='OII EW',ytitle='Number of galaxies',/fill,fcolor=fsc_color('gray'),xtickformat='(A1)',yrange=[0,90],position=[0.1,0.17,0.48,0.95]
      plothist,sci(good).oiiew,bin=2.5,/overplot,color=fsc_color('grn7'),/fill,fcolor=fsc_Color('grn3')
      plothist,sci(quiescent).oiiew,bin=2.5,/overplot,color=fsc_color('grn7'),/fline,fcolor=fsc_Color('grn7'),forientation=45
      axis,yaxis=1,yrange=[0,90],ystyle=1,ytickformat='(a1)'
      axis,xaxis=0,xticks=7,xtickv=[-5,0,5,10,15,20,25],xtickn=['-5','0','5','10','15','20','>22.5']
      axis,xaxis=0,xrange=[-7.5,25],xstyle=1,xtickformat='(A1)'
      xyouts,2,80,'MS0451, z~0.54',charsize=1.2
      sunsym = sunsymbol()
      plothist,sci.logmstar,xrange=[9,12],xtitle='log[M!L*!N/M!N'+sunsym+'!N]',/fill,fcolor=fsc_color('gray'),yrange=[0,90],bin=0.3,position=[0.57,0.17,0.95,0.95],xtickformat='(A1)'
      plothist,sci(good).logmstar,/overplot,color=fsc_color('grn7'),/fill,fcolor=fsc_Color('grn3'),bin=0.3
      plothist,sci(quiescent).logmstar,/overplot,color=fsc_color('grn7'),/fline,fcolor=fsc_Color('grn7'),forientation=45,bin=0.3
      axis,xaxis=0,xrange=[9,12],xstyle=1
   device,/close

      
   stop
end
