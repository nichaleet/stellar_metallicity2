pro demo_cl0024
   sciold = mrdfits('/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz',1)
   sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit99.fits.gz',1);300 rows
   sci.oiiew = sci.oiiew*(-1.)
   wHiEm = where(sci.oiiew gt 22.5)
   sci(whiem).oiiew = 23.

   quiescent = where(sci.oiiew lt 5. and (sci.fuv_v_rest gt 3. or sci.nuv_mag eq -99.) and sci.logmstar gt 6.,cquiescent)
;   good = where(sci.oiiew lt 5. and (sci.fuv_v_rest gt 3. or sci.nuv_mag eq -99.) and sciold.snfit gt 6.83 and sci.logmstar gt 6.,cgood) ;equivalent to sn=10 per ang in rest frame

   good = where(sci.oiiew lt 5. and (sci.fuv_v_rest gt 3. or sci.nuv_mag eq -99.) and sciold.snfit gt 5.5 and sci.logmstar gt 6.,cgood) ;equivalent to sn=8 per ang in rest frame
   set_plot,'ps'
   !p.multi = [0,2,1]
   !p.font = 0
   !p.charsize =1
   ;oii ew histogram
   psname='cl0024_quiescent_histogram.eps'
   device, filename = psname,xsize = 18,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.charsize=1.2

      plothist,sci.oiiew,bin=2.5,xrange=[-7.5,25],xtitle='OII EW',ytitle='Number of galaxies',/fill,fcolor=fsc_color('gray'),xtickformat='(A1)',yrange=[0,90],position=[0.1,0.17,0.48,0.95]
      plothist,sci(good).oiiew,bin=2.5,/overplot,color=fsc_color('grn7'),/fill,fcolor=fsc_Color('grn3')
      plothist,sci(quiescent).oiiew,bin=2.5,/overplot,color=fsc_color('grn7'),/fline,fcolor=fsc_Color('grn7'),forientation=45
      axis,yaxis=1,yrange=[0,90],ystyle=1,ytickformat='(a1)'
      axis,xaxis=0,xticks=7,xtickv=[-5,0,5,10,15,20,25],xtickn=['-5','0','5','10','15','20','>22.5']
      axis,xaxis=0,xrange=[-7.5,25],xstyle=1,xtickformat='(A1)'
      !p.charsize=1.2
      xyouts,0,80,'Cl0024, z~0.39',charsize=1.2

      sunsym = sunsymbol()
      plothist,sci.logmstar,xrange=[9,12],xtitle='log[M!L*!N/M'+sunsym+'!N]',/fill,fcolor=fsc_color('gray'),yrange=[0,90],bin=0.3,position=[0.57,0.17,0.95,0.95],xtickformat='(A1)'
      plothist,sci(good).logmstar,/overplot,color=fsc_color('grn7'),/fill,fcolor=fsc_Color('grn3'),bin=0.3
      plothist,sci(quiescent).logmstar,/overplot,color=fsc_color('grn7'),/fline,fcolor=fsc_Color('grn7'),forientation=45,bin=0.3
      axis,xaxis=0,xrange=[9,12],xstyle=1
      al_legend,['All spectra','Quiescent galaxies','Quiescent galaxies with'],polycolor=['gray','grn7','grn3'],psym=[10,10,10],box=0,symsize=2,line_orientation=[-200,45,-200],position=[9.,89],polyspace=0.1,charsize=1
      xyouts,10,65,'S/N>10',charsize=1
   device,/close

      
   stop
end
