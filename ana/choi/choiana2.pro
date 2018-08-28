pro choiana2
  sunsym = sunsymbol()
  Delta = '!9'+string("104B)+'!x'
  alpha = '!9'+string("141B)+'!x'  ;READ DATA
  ;read science
  dir = '/scr2/nichal/workspace4/sps_fit/data/choi/ages/'
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
  CaFeerr = [0.06,0.03,0.02,0.04,0.05,0.03,0.03,0.15,0.06,0.05]
  choi = {logmstar:massbinval,feh:fehchoi,feherr:fehchoierr,age:agechoi,ageerr:agechoierr,CFe:CFe,NFe:NFe,NFeerr:NFeerr,MgFe:MgFe,MgFeerr:MgFeerr,caFe:CaFe,CaFeerr:CaFeerr}

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;MILES Library abundance (from Conroy17)
  miles = {feh:[-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],$
           ofe:[0.6,0.5,0.5,0.4,0.3,0.2,0.2,0.1,0.0,0.0],$
           mgfe:[0.4,0.4,0.4,0.4,0.34,0.22,0.14,0.11,0.05,0.04],$
           cafe:[0.32,0.30,0.28,0.26,0.26,0.17,0.12,0.06,0.00,0.00]}
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  set_plot,'ps'
  sci1 = mrdfits(dir+'sps_fit01.fits.gz',1,/silent) ;alpha = {Mg,Si,Ca,Ti,O} fixed
  sci3 = mrdfits(dir+'sps_fit03.fits.gz',1,/silent) ;alpha = '[Mg/Fe]
  sci9 = mrdfits(dir+'sps_fit09.fits.gz',1,/silent) ;6para, alpha = Mg, O
  sci10 = mrdfits(dir+'sps_fit10.fits.gz',1,/silent);6para, alpha = Mg, N
  sci11 = mrdfits(dir+'sps_fit11.fits.gz',1,/silent);13 para alpha = Mg,N,Fe,O,C,N,Si,Ca,Ti

  ;test if the str match
  if total(sci1.fehchoi-choi.feh) ne 0. then stop, 'structure might not matched'

  ;fix the library abundance
;  sci1.alphafe = sci1.alphafe+interpol(miles.mgfe,miles.feh,sci1.feh)
;  sci3.alphafe = sci3.alphafe+interpol(miles.mgfe,miles.feh,sci3.feh)
;  sci9.alphafe = sci9.alphafe+interpol(miles.mgfe,miles.feh,sci9.feh)
;  sci10.alphafe = sci10.alphafe+interpol(miles.mgfe,miles.feh,sci10.feh)
;  sci11.alphafe = sci11.alphafe+interpol(miles.mgfe,miles.feh,sci11.feh)
;  sci11.pfit[4] = sci11.pfit[4]+interpol(miles.mgfe,miles.feh,sci11.feh)
;  sci11.pfit[11] = sci11.pfit[11]+interpol(miles.cafe,miles.feh,sci11.feh)

  sci1.age = alog10(sci1.age)+9.
  sci3.age = alog10(sci3.age)+9.
  sci9.age = alog10(sci9.age)+9.
  sci10.age = alog10(sci10.age)+9.
  sci11.age = alog10(sci11.age)+9.
 
  sci1.ageupper = alog10(sci1.ageupper)+9.
  sci3.ageupper = alog10(sci3.ageupper)+9.
  sci9.ageupper = alog10(sci9.ageupper)+9.
  sci10.ageupper = alog10(sci10.ageupper)+9.
  sci11.ageupper = alog10(sci11.ageupper)+9.

  sci1.agelower = alog10(sci1.agelower)+9.
  sci3.agelower = alog10(sci3.agelower)+9.
  sci9.agelower = alog10(sci9.agelower)+9.
  sci10.agelower = alog10(sci10.agelower)+9.
  sci11.agelower = alog10(sci11.agelower)+9.

  psname = 'all_choicompare_age.eps'
  device,filename=psname,xsize=12,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
    !p.font = 0
    !p.charsize = 1.2
    !P.position = [0.15,0.15,0.95,0.9] 
    logagejc = alog10(sci1.agechoi)+9.
    logagejc_lower = alog10(sci1.agechoi-sci1.agechoierr)+9.
    logagejc_upper = alog10(sci1.agechoi+sci1.agechoierr)+9.
  
    plot,logagejc,sci1.age,psym=1,xtitle='log(Age!DChoi14!N)',ytitle='log(Age!Dthis work!N)',$
       xrange=[9.3,9.8],yrange=[9.4,9.8],/nodata
    ;cgerrplot,logagejc,logagenl_lower,logagenl_upper
    cgerrplot,sci1.age,logagejc_lower,logagejc_upper,/horizontal,color='firebrick'
    cgerrplot,sci3.age,logagejc_lower,logagejc_upper,/horizontal,color='chocolate'
    cgerrplot,sci9.age,logagejc_lower,logagejc_upper,/horizontal,color='blu5'
    cgerrplot,sci10.age,logagejc_lower,logagejc_upper,/horizontal,color='lightseagreen'
    cgerrplot,sci11.age,logagejc_lower,logagejc_upper,/horizontal,color='darkorchid'

    cgplot,logagejc,sci1.age,/overplot,psym=18,color=fsc_color('firebrick')
    cgplot,logagejc,sci3.age,/overplot,psym=15,color=fsc_color('chocolate')
    cgplot,logagejc,sci9.age,/overplot,psym=16,color=fsc_color('blu5')
    cgplot,logagejc,sci10.age,/overplot,psym=14,color=fsc_color('lightseagreen')
    cgplot,logagejc,sci11.age,/overplot,psym=46,color=fsc_color('darkorchid')
    
    oplot,!x.crange,!x.crange,linestyle=2	
  device,/close

  psname = 'all_choicompare_feh.eps'
  device,filename=psname,xsize=12,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  
    plot,choi.feh,sci1.feh,psym=1,xtitle='log([Fe/H]!DChoi14!N)',ytitle='log([Fe/H]!Dthis work!N)',$
       xrange=[-0.2,0.1],yrange=[-0.2,0.2],/nodata
    cgplot,choi.feh,sci1.feh,/overplot,psym=18,color=fsc_color('firebrick')
    cgplot,choi.feh,sci3.feh,/overplot,psym=15,color=fsc_color('chocolate')
    cgplot,choi.feh,sci9.feh,/overplot,psym=16,color=fsc_color('blu5')
    cgplot,choi.feh,sci10.feh,/overplot,psym=14,color=fsc_color('lightseagreen')
    cgplot,choi.feh,sci11.feh,/overplot,psym=46,color=fsc_color('darkorchid')
    
    oplot,!x.crange,!x.crange,linestyle=2	
  device,/close

  psname = 'all_choicompare_mgfe.eps'
  device,filename=psname,xsize=12,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

    plot,choi.mgfe,sci1.alphafe,psym=1,xtitle='log([Mg/Fe]!DChoi14!N)',ytitle='log([Mg/Fe]!Dthis work!N)',$
       xrange=[0.0,0.4],yrange=[0.0,0.4],/nodata
    cgplot,choi.mgfe,sci1.alphafe,/overplot,psym=18,color=fsc_color('firebrick')
    cgplot,choi.mgfe,sci3.alphafe,/overplot,psym=15,color=fsc_color('chocolate')
    cgplot,choi.mgfe,sci9.alphafe,/overplot,psym=16,color=fsc_color('blu5')
    cgplot,choi.mgfe,sci10.alphafe,/overplot,psym=14,color=fsc_color('lightseagreen')
    cgplot,choi.mgfe,sci11.pfit[4],/overplot,psym=46,color=fsc_color('darkorchid')

    oplot,!x.crange,!x.crange,linestyle=2
  device,/close

end
