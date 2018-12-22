pro nta_leetho18
;calculate nta as a fn of mass based on Lu+2015 and functional form of MZR from NL18
   solarfeh = 0.02
   R = 0.46 ;return fraction
   yield = 0.07 ;kinda high but follow Lu+2015

   ;solar abundances from asplund2009 (number of atoms)
   solarfeh = 10^(7.54-12.)*55.
   solarmgh = 10^(7.64-12.)*24.
   yieldfeh = 3.*solarfeh
   yieldmgh = 3.*solarmgh

   nmass = 50.
   mass = 9.55+findgen(nmass)/nmass*1.9
   logfeh = -0.05+0.16*(mass-10.)  
   scatter = 0.1
   feh=solarfeh*10.^(logfeh)
   fehhi = solarfeh*10.^(logfeh+scatter)
   fehlo = solarfeh*10.^(logfeh-scatter) ;assume intrinsic scatter=0.1

   nta = (1.-R)*(yieldfeh/feh-1.)
   ntahi = (1.-R)*(yieldfeh/fehlo-1.)
   ntalo = (1.-R)*(yieldfeh/fehhi-1.)

   logmgh = 0.09+0.10*(mass-10.) ;mgh
   scattermgh = 0.2 ;mgh
   mgh = solarmgh*10.^(logmgh)
   mghhi = solarmgh*10.^(logmgh+scattermgh)
   mghlo = solarmgh*10.^(logmgh-scattermgh) ;assume intrinsic scatter=0.1
   ntamgh = (1.-R)*(yieldmgh/mgh-1.)
   ntamghhi = (1.-R)*(yieldmgh/mghlo-1.)
   ntamghlo = (1.-R)*(yieldmgh/mghhi-1.)
 
   logmghlowz = 0.15+0.09+0.10*(mass-10.) ;mghlowz
   mghlowz = solarmgh*10.^(logmghlowz)
   ntamghlowz = (1.-R)*(yieldmgh/mghlowz-1.)

   ;functional form
   lognta = alog10(nta)
   linpar = linfit(mass,lognta,yfit=fitlognta)
   print, 'dependence of nta on M', linpar(1)
   logntahi = alog10(ntahi)
   logntalo = alog10(ntalo)

   logntamgh = alog10(ntamgh)
   linparmgh = linfit(mass,logntamgh,yfit=fitlogntamgh)
   print, 'dependence of ntamgh on M', linparmgh(1)
   logntamghhi = alog10(ntamghhi)
   logntamghlo = alog10(ntamghlo)

   logntamghlowz = alog10(ntamghlowz)
   linparmghlowz = linfit(mass,logntamghlowz,yfit=fitlogntamghlowz)
   print, 'dependence of ntamgh localz on M', linparmghlowz(1)

   ;below is a eyeball fit to Gallazzi2005
   mass2 = 9.05+findgen(nmass)/nmass*2.4
   logfeh2 = -0.05+0.16*(mass2-10.)
   lowmass = where(mass2 lt 10.5)
   logfeh2(lowmass) = 0.03+0.8/2.5*(mass2(lowmass)-10.5+0.03)
   feh2=solarfeh*10.^(logfeh2)
   fehhi2 = solarfeh*10.^(logfeh2+0.1)
   fehlo2 = solarfeh*10.^(logfeh2-0.1) ;assume intrinsic scatter=0.1

   nta2 = (1.-R)*(yield/feh2-1.)
   ntahi2 = (1.-R)*(yield/fehlo2-1.)
   ntalo2 = (1.-R)*(yield/fehhi2-1.)
   ;functional form
   lognta2 = alog10(nta2)
   linpar2 = linfit(mass2,lognta2,yfit=fitlognta2)
   print, 'dependence of nta2 on M', linpar2(1)
   logntahi2 = alog10(ntahi2)
   logntalo2 = alog10(ntalo2)
   
   set_plot,'ps'
   !p.font = 0
   ntaletter = "150B
   proptoletter = "265B
   sunsym = sunsymbol()

   psname='nta_mass_mg.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,mass,nta,xtitle='Log(M/M'+sunsym+')',ytitle='!9'+string(ntaletter)+'!x',/ylog,/nodata,xrange=[9.5,11.5],yrange=[0.1,10]
      x=[mass,reverse(mass)]
      y=[ntahi,reverse(ntalo)]>0.1
      ymg=[ntamghhi,reverse(ntamghlo)]>0.1
      polyfill,x,y,color=fsc_color('rose')
      polyfill,x,ymg,color=fsc_color('powderblue'),/line_fill,orientation=45
      polyfill,x,ymg,color=fsc_color('powderblue'),/line_fill,orientation=135
      oplot,mass,nta
      oplot,mass,ntamgh
      oplot,mass,10.^fitlognta,color=fsc_color('red'),linestyle=5
      oplot,mass,10.^fitlogntamgh,color=fsc_color('navy'),linestyle=5
      ;oplot,mass,ntamghlowz
      ;oplot,mass,10.^fitlogntamghlowz,color=fsc_color('purple'),linestyle=5


      xyouts,11,6,'!9'+string(ntaletter)+' !9'+string(proptoletter)+' !xM!D*!N!E'+sigfig(linpar(1),2)
   device,/close

   psname='nta_mass_hayward.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,mass,nta,xtitle='Log(M/M'+sunsym+')',ytitle='!9'+string(ntaletter)+'!x',/ylog,/nodata,xrange=[6,12],xstyle=1,yrange=[0.01,500],ystyle=1
      x=[mass,reverse(mass)]
      y=[ntahi,reverse(ntalo)]>0.01
      ymg=[ntamghhi,reverse(ntamghlo)]>0.1
      polyfill,x,y,color=fsc_color('rose')
      polyfill,x,ymg,color=fsc_color('powderblue'),/line_fill,orientation=45
      polyfill,x,ymg,color=fsc_color('powderblue'),/line_fill,orientation=135

      oplot,mass,nta
      oplot,mass,10.^fitlognta,color=fsc_color('red'),linestyle=5
      oplot,mass,ntamgh
      oplot,mass,10.^fitlogntamgh,color=fsc_color('navy'),linestyle=5

   device,/close


   psname='nta2_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,mass2,nta2,xtitle='Log(M/M'+sunsym+')',ytitle='!9'+string(ntaletter)+'!x',/ylog,/nodata,xrange=[9.,11.5]
      x=[mass2,reverse(mass2)]
      y=[ntahi2,reverse(ntalo2)]
      polyfill,x,y,color=fsc_color('cyan')
      oplot,mass2,nta2
      oplot,mass2,10.^fitlognta2,color=fsc_color('navy'),linestyle=5
      xyouts,11,6,'!9'+string(ntaletter)+' !9'+string(proptoletter)+' !xM!D*!N!E'+sigfig(linpar2(1),2)
   device,/close

   stop
end
