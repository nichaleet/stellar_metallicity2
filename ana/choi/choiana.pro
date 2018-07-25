pro choiana
  ;READ DATA
  ;read science

  file = '/scr2/nichal/workspace4/sps_fit/data/choi/ages/sps_fit01.fits.gz'
  science = mrdfits(file,1,/silent)

  ;choi data ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  massbinval = [10.5,10.8,11.0,11.3,10.8,11.1,11.3,10.9,11.0,11.3]
  fehchoi    = [-0.11,-0.05,-0.02,-0.03,-0.07,-0.04,-0.05,-0.15,-0.02,-0.05]
  fehchoierr = [0.03,0.01,0.01,0.02,0.02,0.01,0.02,0.07,0.03,0.03]
  agechoi    = [3.27,3.47,4.55,5.61,2.99,3.28,4.00,2.67,2.49,3.06]
  agechoierr = [0.21,0.12,0.13,0.15,0.15,0.13,0.18,0.20,0.12,0.11]
  CFe  = [0.19,0.16,0.18,0.19,0.18,0.18,0.14,0.16,0.24,0.26]
  NFe  = [0.26,0.25,0.24,0.26,0.26,0.34,0.25,0.18,0.58,0.73]
  NFeerr = [0.10,0.05,0.03,0.05,0.08,0.05,0.06,0.30,0.09,0.09]
  MgFe = [0.18,0.22,0.21,0.29,0.23,0.23,0.30,0.05,0.09,0.19]
  MgFeerr = [0.04,0.02,0.01,0.03,0.04,0.02,0.03,0.13,0.05,0.04]
  CaFe = [0.04,0.03,0.04,0.04,0.05,-0.03,0.08,0.06,-0.03,0.00]
  choi = {logmstar:massbinval,feh:fehchoi,feherr:fehchoierr,age:agechoi,ageerr:agechoierr,CFe:CFe,NFe:NFe,NFeerr:NFeerr,MgFe:MgFe,MgFeerr:MgFeerr}
stop

  agediff = science.age-science.agechoi
  fehdiff = science.feh-science.fehchoi

  agediffupper = ageupper-science.agechoi
  fehdiffupper = fehupper-science.fehchoi
  agedifflower = agelower-science.agechoi
  fehdifflower = fehlower-science.fehchoi

  set_plot,'ps'
  !p.charsize=1.2
  basename= file_basename(file)
  typename= strmid(basename,0,strpos(basename,'_sps_fit.fits.gz'))

  psname = typename+'_choi_comparison2.eps'
  device,filename=psname,xsize=21,ysize=7,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font=0
  !p.charsize=0
  !p.multi=[0,2,1]
  logagenl = alog10(science.age)+9.
  logagenl_upper = alog10(ageupper)+9.
  logagenl_lower = alog10(agelower)+9.

  logagejc = alog10(science.agechoi)+9.
  logagejc_lower = alog10(science.agechoi-science.agechoierr)+9.
  logagejc_upper = alog10(science.agechoi+science.agechoierr)+9.

  plot,logagejc,logagenl,psym=1,xtitle='log(Age!DChoi14!N)',ytitle='log(Age!Dthis work!N)',$
     xrange=[9.3,9.8],yrange=[9.4,9.8],/nodata
  cgerrplot,logagejc,logagenl_lower,logagenl_upper
  cgerrplot,logagenl,logagejc_lower,logagejc_upper,/horizontal
  vsym,4,/fill,rot=45
  oplot,logagejc,logagenl,psym=8,color=fsc_color('red')
  oplot,!x.crange,!x.crange,linestyle=2	
  
  plot,science.fehchoi,science.feh,xtitle='[Fe/H]!DChoi14!N',ytitle='[Fe/H]!Dthis work!N',psym=1,$
      xrange=[-0.25,0.05],yrange=[-0.3,0.2],/nodata
  cgerrplot,science.fehchoi,fehlower,fehupper
  cgerrplot,science.feh,science.fehchoi-science.fehchoierr,science.fehchoi+science.fehchoierr,/horizontal
  oplot,science.fehchoi,science.feh,psym=8,color=fsc_color('red')
  oplot,!x.crange,!x.crange,linestyle=2  
device,/close
  

end
