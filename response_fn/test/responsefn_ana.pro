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
   specalpha = add_response(spslambda,spsspec,[zmet,age],['Mg','O','Si','Ca','Ti'],rebin(abund,5))
   specmg = add_response(spslambda,spsspec,[zmet,age],['Mg'],abund)
   speco = add_response(spslambda,spsspec,[zmet,age],['O'],abund)
   specsi = add_response(spslambda,spsspec,[zmet,age],['Si'],abund)
   specca = add_response(spslambda,spsspec,[zmet,age],['Ca'],abund)
   specTi = add_response(spslambda,spsspec,[zmet,age],['Ti'],abund)
   set_plot,'ps'
   !p.font=0
   !p.multi = [0,1,1]
   epsname = 'responsefn_compare.eps'
   device, filename = epsname,xsize = 15,ysize = 12, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,spslambda,spsspec,xtitle='lambda',ytitle='response fn',$
           yrange=[0.8,1.1],ystyle=1,xrange=minmax(spslambda),xstyle=1,$
           title='[Z/H]:'+sigfig(zmet,2)+' Age:'+sigfig(age,2)+' Gyr'
      oplot,spslambda,specalpha,color=fsc_color('maroon')
      oplot,spslambda,specmg,color=fsc_color('teal')
      oplot,spslambda,speco,color=fsc_color('navy')
      oplot,spslambda,specsi,color=fsc_color('darkorchid')
      oplot,spslambda,specca,color=fsc_color('salmon')
      oplot,spslambda,specca,color=fsc_color('olivedrab')
      al_legend,['all','Mg','O','Si','Ca','Ti'],psym=cgsymcat(15),$
           colors=fsc_color(['maroon','teal','navy','darkorchid','salmon','olivedrab']),/right,/bottom
   device,/close

end
