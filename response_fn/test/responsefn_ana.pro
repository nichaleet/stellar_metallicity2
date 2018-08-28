pro responsefn_ana
   ;make flat original spectrum
   redshift = 0.55
   dwl = 0.65 ;dispersion per pixel in angstrom
   iwave = 3800.
   nwl = 4096
   lambdaobs = findgen(nwl)*dwl+iwave*(1.+redshift)
   spslambda = lambdaobs/(1.+redshift)
   spsspec = fltarr(nwl)+1.

   abund = [0.25]
   zmet = -0.2
   age = 3 ;gyr
   ;add response fn
;   specalpha = add_response(spslambda,spsspec,[zmet,age],['Mg','O','Si','Ca','Ti'],rebin(abund,5))
   specmg = add_response(spslambda,spsspec,[zmet,age],['Mg'],abund)
   speco = add_response(spslambda,spsspec,[zmet,age],['O'],abund)
   specsi = add_response(spslambda,spsspec,[zmet,age],['Si'],abund)
   specca = add_response(spslambda,spsspec,[zmet,age],['Ca'],abund)
   specTi = add_response(spslambda,spsspec,[zmet,age],['Ti'],abund)
   specN = add_response(spslambda,spsspec,[zmet,age],['N'],abund)
   specC = add_response(spslambda,spsspec,[zmet,age],['C'],abund)

   set_plot,'ps'
   !p.font=0
   !p.multi = [0,1,1]
   epsname = 'responsefn_compare.eps'
   device, filename = epsname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,spslambda,spsspec,xtitle='wavelength (Angstrom)',ytitle='response fn',$
           yrange=[0.8,1.05],ystyle=1,xrange=minmax(spslambda),xstyle=1;,$
           ;title='[Z/H]:'+sigfig(zmet,2)+' Age:'+sigfig(age,2)+' Gyr'
      oplot,spslambda,specmg,color=fsc_color('teal')
      oplot,spslambda,speco,color=fsc_color('navy')
      oplot,spslambda,specsi,color=fsc_color('darkorchid')
      oplot,spslambda,specca,color=fsc_color('salmon')
      oplot,spslambda,specTi,color=fsc_color('olivedrab')
      oplot,spslambda,specN,color=fsc_color('maroon')
      al_legend,['Mg','O','Si','Ca','Ti','N'],psym=cgsymcat(15),$
           colors=['teal','navy','darkorchid','salmon','olivedrab','maroon'],/right,/bottom
   device,/close

   abund = findgen(5)*(0.8-0.)/4.+0.
   specNarr = fltarr(nwl,5)
   for i=0,4 do specNarr[*,i] = add_response(spslambda,spsspec,[zmet,age],['N'],[abund[i]])
   
   abund = findgen(5)*(0.2+0.2)/4.-0.2
   specFearr = fltarr(nwl,5)
   for i=0,4 do specFearr[*,i] = add_response(spslambda,spsspec,[zmet,age],['Fe'],[abund[i]])

   color = ['red','salmon','darkgreen','navy','purple']
   epsname = 'responsefn_N.eps'
   device,filename = epsname,xsize = 15,ysize = 12, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plot,spslambda,specNarr[*,2],/nodata,xtitle='lambda',ytitle='N response fn',$
            yrange=[0.8,1.1],ystyle=1,xrange=minmax(spslambda),xstyle=1,$
            title='[Z/H]:'+sigfig(zmet,2)+' Age:'+sigfig(age,2)+' Gyr'
       for i=0,4 do oplot,spslambda,specNarr[*,i],color=fsc_color(color[i])
   device,/close

   epsname = 'responsefn_Fe.eps'
   device,filename = epsname,xsize = 15,ysize = 12, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       plot,spslambda,specFearr[*,2],/nodata,xtitle='lambda',ytitle='Fe response fn',$
            yrange=[0.8,1.1],ystyle=1,xrange=minmax(spslambda),xstyle=1,$
            title='[Z/H]:'+sigfig(zmet,2)+' Age:'+sigfig(age,2)+' Gyr'
       for i=0,4 do oplot,spslambda,specFearr[*,i],color=fsc_color(color[i])
   device,/close

end
