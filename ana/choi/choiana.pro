pro choiana
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
  choi = {logmstar:massbinval,feh:fehchoi,feherr:fehchoierr,age:agechoi,ageerr:agechoierr,CFe:CFe,NFe:NFe,NFeerr:NFeerr,MgFe:MgFe,MgFeerr:MgFeerr,caFe:CaFe,CaFeerr:CaFeerr,alpha:fltarr(10),alphaerr:fltarr(10)}

  ;solar abundance (Asplund2009)
  sol ={mg:7.60,mgerr:0.04,ca:6.34,caerr:0.04,o:8.69,oerr:0.05,si:7.51,sierr:0.03}
  ;calculate alpha abundance from choi data
  ;weight of solar abundance
  dummy = 6.34
  wsol = {mg:10.^(sol.mg-dummy),mgerr:abs(alog(10.)*sol.mgerr*10.^(sol.mg-dummy)),$
          ca:10.^(sol.ca-dummy),caerr:abs(alog(10.)*sol.caerr*10.^(sol.ca-dummy)),$
           o:10.^(sol.o-dummy),oerr:abs(alog(10.)*sol.oerr*10.^(sol.o-dummy)),$
          si:10.^(sol.si-dummy),sierr:abs(alog(10.)*sol.sierr*10.^(sol.si-dummy))}

  tempalpha = (wsol.mg*10.^choi.mgfe+wsol.ca*10.^choi.cafe)/(wsol.mg+wsol.ca)
  totalw = (wsol.mg+wsol.ca)
  temperr = sqrt(choi.mgfeerr^2*(wsol.mg*(10.^choi.mgfe)*alog(10.)/totalw)^2+$
                 wsol.mgerr^2*(wsol.ca*((10.^choi.mgfe)-(10.^choi.cafe))/totalw^2)^2+$
                 choi.cafeerr^2*(wsol.ca*(10.^choi.cafe)*alog(10.)/totalw)^2+$
                 wsol.caerr^2*(wsol.ca*((10.^choi.cafe)-(10.^choi.mgfe))/totalw^2)^2) 
  choi.alpha = alog10(tempalpha)
  choi.alphaerr = abs(temperr/tempalpha/alog(10.))  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  set_plot,'ps'
 
  for i=1,3 do begin
     case i of 
       1:begin
         name = '5para_5elements'
         sci = mrdfits(dir+'sps_fit01.fits.gz',1,/silent)
         doalpha = 1
         end
       2:begin
         name = '4para'
         sci = mrdfits(dir+'sps_fit02.fits.gz',1,/silent)
         doalpha = 0
         end
       3:begin
         name = '5para_mgonly'
         sci = mrdfits(dir+'sps_fit03.fits.gz',1,/silent)
         doalpha = 1
         end
       else:stop,'something is wrong'
     endcase

     ;test if the str match
     if total(sci.fehchoi-choi.feh) ne 0. then stop, 'structure might not matched'

     psname = name+'_choicompare_agemetal.eps'
     device,filename=psname,xsize=20,ysize=7,$
            xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       !p.font=0
       !p.charsize=1.2
       !p.multi=[0,2,1]
       logagenl = alog10(sci.age)+9.
       logagenl_upper = alog10(sci.ageupper)+9.
       logagenl_lower = alog10(sci.agelower)+9.
   
       logagejc = alog10(sci.agechoi)+9.
       logagejc_lower = alog10(sci.agechoi-sci.agechoierr)+9.
       logagejc_upper = alog10(sci.agechoi+sci.agechoierr)+9.
   
       plot,logagejc,logagenl,psym=1,xtitle='log(Age!DChoi14!N)',ytitle='log(Age!Dthis work!N)',$
          xrange=[9.3,9.8],yrange=[9.4,9.8],/nodata
       cgerrplot,logagejc,logagenl_lower,logagenl_upper
       cgerrplot,logagenl,logagejc_lower,logagejc_upper,/horizontal
       vsym,4,/fill,rot=45
       oplot,logagejc,logagenl,psym=8,color=fsc_color('red')
       oplot,!x.crange,!x.crange,linestyle=2	
       
       plot,sci.fehchoi,sci.feh,xtitle='[Fe/H]!DChoi14!N',ytitle='[Fe/H]!Dthis work!N',psym=1,$
           xrange=[-0.25,0.05],yrange=[-0.3,0.2],/nodata
       cgerrplot,sci.fehchoi,sci.fehlower,sci.fehupper
       cgerrplot,sci.feh,sci.fehchoi-sci.fehchoierr,sci.fehchoi+sci.fehchoierr,/horizontal
       oplot,sci.fehchoi,sci.feh,psym=8,color=fsc_color('red')
       oplot,!x.crange,!x.crange,linestyle=2  
     device,/close
    
     if doalpha eq 1 then begin 
     psname = name+'_choicompare_alpha.eps'
     !p.multi=[0,3,1]
     device,filename=psname,xsize=30,ysize=7,$
            xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
        plot,choi.alpha,sci.alphafe,xtitle='[{Mg,Ca}/Fe]!DChoi14!N',$
             ytitle='[{Mg,Si,Ca,Ti,O}/Fe]!Dthis work!N',/nodata
        cgerrplot,choi.alpha,sci.alphafelower,sci.alphafeupper
        cgerrplot,sci.alphafe,choi.alpha-choi.alphaerr,choi.alpha+choi.alphaerr,/horizontal
        oplot,choi.alpha,sci.alphafe,psym=8,color=fsc_color('red')
        oplot,!x.crange,!x.crange,linestyle=2
   
        plot,choi.mgfe,sci.alphafe,xtitle='[Mg/Fe]!DChoi14!N',$
             ytitle='[{Mg,Si,Ca,Ti,O}/Fe]!Dthis work!N',/nodata
        cgerrplot,choi.mgfe,sci.alphafelower,sci.alphafeupper
        cgerrplot,sci.alphafe,choi.mgfe-choi.mgfeerr,choi.mgfe+choi.mgfeerr,/horizontal
        oplot,choi.mgfe,sci.alphafe,psym=8,color=fsc_color('red')
        oplot,!x.crange,!x.crange,linestyle=2
   
        plot,choi.cafe,sci.alphafe,xtitle='[Ca/Fe]!DChoi14!N',$
             ytitle='[{Mg,Si,Ca,Ti,O}/Fe]!Dthis work!N',/nodata
        cgerrplot,choi.cafe,sci.alphafelower,sci.alphafeupper
        cgerrplot,sci.alphafe,choi.cafe-choi.cafeerr,choi.cafe+choi.cafeerr,/horizontal
        oplot,choi.cafe,sci.alphafe,psym=8,color=fsc_color('red')
        oplot,!x.crange,!x.crange,linestyle=2
   
     device,/close
     endif
  endfor
end
