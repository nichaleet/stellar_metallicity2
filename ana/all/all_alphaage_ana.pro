function get_refvalue,dummy,leetho18=leetho18,gallazzi05=gallazzi05,choi14=choi14,$
              conroy14=conroy14,derossi17=derossi17,sybilska17=sybilska17,$
              naiman18=naiman18,saracco19=saracco19,vincenzo18=vincenzo18,$
              zahid12=zahid12,nelson19=nelson19
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   ;REFERENCE VALUES FOR PLOTTING
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if keyword_set(leetho18) then begin
      restore,'/scr2/nichal/workspace4/ana/ms0451/leetho18a_avevals.sav'
      outstr = {hifeh_cl0024:hifeh,lofeh_cl0024:lofeh,bndry_mass_cl0024:bndry_mass,$
                HIFEH_SDSS:HIFEH_SDSS,LOFEH_SDSS:LOFEH_SDSS,BNDRY_MASS_SDSS:BNDRY_MASS_SDSS}
      return,outstr
   endif
 
   ;Get the average values measured in Gallazzi05
   if keyword_set(gallazzi05) then begin
      mass = [9.00,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5];the first mass is actually 8.91
      feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
      ;below are the stars in the figure 8 of Gallazzi05
      feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
      fehlo = g05_feh-g05_feherr
      fehhi = g05_feh+g05_feherr
      outstr = {mass:mass,feh:feh,feherr:feherr,fehlo:fehlo,fehhi:fehhi}
      return,outstr
   endif

   ;Choi's Data
   if keyword_set(choi14) then begin
      Choi14z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],$
                   Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01],$
                   mgfe:[0.04,0.13,0.16,0.21,0.23],mgfeerr:[0.05,0.03,0.01,0.01,0.02],$
                   age:[2.8,3.57,4.25,5.5,5.84],ageerr:[0.13,0.17,0.10,0.6,0.14]}
      Choi14z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],$
                   Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02],$
                   mgfe:[0.17,0.19,0.20,0.22,0.23],mgfeerr:[0.06,0.02,0.01,0.02,0.04],$
                   age:[3.05,3.16,4.24,4.72,6.24],ageerr:[0.24,0.12,0.11,0.14,0.23]}
      Choi14z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3,-99],$
                   Feh:[-0.11,-0.05,-0.02,-0.03,-99],Feherr:[0.03,0.01,0.01,0.02,-99],$
                   mgfe:[0.18,0.22,0.21,0.29,-99],mgfeerr:[0.04,0.02,0.01,0.03,-99],$
                   age:[3.27,3.47,4.55,5.61,-99],ageerr:[0.21,0.12,0.13,0.15,-99]}
      Choi14z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3,-99,-99],$
                   Feh:[-0.07,-0.04,-0.05,-99,-99],Feherr:[0.02,0.01,0.02,-99,-99],$
                   mgfe:[0.23,0.23,0.30,-99,-99],mgfeerr:[0.04,0.02,0.03,-99,-99],$
                   age:[2.99,3.28,4.00,-99,-99],ageerr:[0.15,0.13,0.18,-99,-99]}
      Choi14z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3,-99,-99],$
                   Feh:[-0.15,-0.02,-0.05,-99,-99],Feherr:[0.07,0.03,0.03,-99,-99],$
                   mgfe:[0.05,0.09,0.19,-99,-99],mgfeerr:[0.13,0.05,0.04,-99,-99],$
                   age:[2.67,2.49,3.06,-99,-99],ageerr:[0.2,0.12,0.11,-99,-99]}
      return, [Choi14z01,Choi14z02,Choi14z03,Choi14z04,Choi14z06]
   endif
   
   if keyword_set(conroy14) then begin
      conroy = {mass:[9.63,9.75,9.80,10.08,10.55,10.70,11.07],Feh:[-0.07,-0.05,-0.04,-0.03,-0.01,-0.02,0.],$
                MgFe:[.05,0.08,0.09,0.12,0.15,0.20,0.22]}
      return, conroy
   endif
   if keyword_set(saracco19) then begin ;7 quiescent galaxy cluster members at z=1.22
      saracco = {mass:[10.94,10.54,11.5,10.08,10.8,10.75,11.2],$
                 zh:[0.09,0.16,0.21,-0.27,-0.18,-0.13,0.05],$
                 zh_upper:[0.13,0.06,0.01,0.49,0.44,0.39,0.17],$
                 zh_lower:[-0.42,-0.16,-0.21,-0.12,-0.21,-0.26,-0.05],$,
                 have_em:[0,0,0,0,1,1,0],$
                 age:[2.7,3.4,2.3,1.3,1.7,1.2,1.7]} 
      saracco.zh_upper = saracco.zh+saracco.zh_upper
      saracco.zh_lower = saracco.zh+saracco.zh_lower
      return,saracco
   endif
   if keyword_set(naiman18) then begin
       naiman = {mass_mgfe:[10.17,10.636,10.71,10.90,10.98,11.04],$
                 mgfe:[0.31,0.31,0.32,0.32,0.33,0.33],$
                 mass_feh:[10.30,10.37,10.44,10.83,10.90,10.97,11.03],$	
                 feh:[0.02,0.05,0.04,0.04,0.02,0.04,0.04],$
                 mass_mgh:fltarr(7),$
                 mgh:fltarr(7),$
                 mass_mgfe_contour:[10.1,10.1,10.3,10.3,10.6,10.6,10.9,10.9,$
                                 11.0,11.0,11.1,11.1,10.9,10.9,10.1],$
                 mgfe_contour:[0.30,0.31,0.31,0.34,0.34,0.36,0.36,0.38,$
                           0.38,0.36,0.36,0.32,0.32,0.30,0.30],$
                 mass_feh_contour:[10.27,10.27,10.40,10.40,10.53,10.53,10.67,10.67,10.87,10.87,10.93,10.93,11.07,11.07,11.00,11.00,11.07,11.07,11.01,11.01,10.94,10.94,10.81,10.81,10.73,10.73,10.54,10.54,10.47,10.47,10.40,10.40,10.34,10.34,10.27],$
                 feh_contour:[0.0,0.06,0.06,0.08,0.08,0.11,0.04,0.10,0.10,0.08,0.08,0.09,0.09,0.08,0.08,0.06,0.06,0.02,0.02,-0.02,-0.02,0.00,0.00,-0.02,-0.02,0.0,0.0,0.0,0.02,0.02,0.0,0.0,0.02,0.02,0.00,0.00]}
                 naiman.mgh = interpol(naiman.mgfe,naiman.mass_mgfe,naiman.mass_feh)+naiman.feh
                 naiman.mass_mgh = naiman.mass_feh
      return, naiman
   endif
   if keyword_set(ma16) then begin
      ;Xiangcheng's FIRE data
      readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0pt8.txt',mass08,feh08
      readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0.txt',mass0,feh0
      feh08 = feh08-0.2
      feh0 = feh0-0.2
      outstr = {mass08:mass08,feh08:feh08,mass0:mass0,feh0:feh0}
      return, outstr
   endif

   if keyword_set(derossi17) then begin
   ;De Rossi 2017 (EAGLE)
   ;read off from Figure 5
      DeRossi_z0={z:0.,mass:[9.15,9.48,9.81,10.15,10.5,10.72],feh:[-0.17,-0.08,0.05,0.12,0.17,0.31],feherr:[0.075,0.08,0.07,0.07,0.07,0.03]}
      DeRossi_z1={z:1.,mass:[9.18,9.47,9.85,10.15,10.55,-99],feh:[-0.39,-0.28,-0.09,0.08,0.17,-99],feherr:[0.06,0.06,0.08,0.07,0.085,0]}
      return, [DeRossi_z0,DeRossi_z1]
   endif
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Sybilska et al 2017 (hELENa, IFU from Sauron data)
   if keyword_set(sybilska17) then begin
      aa=read_csv('/scr2/nichal/workspace2/catalogs/sybilska.csv',n_table_header=1,header=header)
      Syb_z0 = {mass:aa.field1,feh:aa.field2}
      return, syb_z0
   endif
   if keyword_set(sybilska17_2) then begin

   endif
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Lu 2014
   if keyword_set(lu14) then begin
      Lu_z0 = {model:'luz0',mass:[8.21,8.61,9.0,9.41,9.8,10.2,10.61,11],$
              feh:[-0.53,-0.43,-0.35,-0.26,-0.16,-0.098,-0.075,-0.09],$
              feherr:[0.13,0.12,0.11,0.12,0.11,0.10,0.095,0.095]}
      Lu_z1 = {model:'luz',mass:[8.21,8.61,9.0,9.4,9.8,10.2,10.6,11],$
              feh:[-0.56,-0.44,-0.35,-0.25,-0.18,-0.13,-0.11,-0.11],$
              feherr:[0.14,0.13,0.11,0.11,0.11,0.10,0.10,0.09]}
      Lu_somerville_z0 = {model:'lusomervillez0',mass:[-99,8.25,8.65,9.05,9.45,9.85,10.25,10.64],$
                         feh:[-99,-1.00,-0.82,-0.64,-0.46,-0.33,-0.15,-0.03],$
                         feherr:[0,0.06,0.1,0.1,0.1,0.1,0.10,0.095]}
      Lu_lu_z0 = {model:'luluz0',mass:[8.73,8.95,9.18,9.47,9.76,10.09,10.46,10.88],$
                 feh:[-.99,-0.84,-0.70,-0.52,-0.31,-0.06,0.18],$
                 feherr:[0.1,0.1,0.1,0.1,0.1,0.1,0.10,0.095]}
      return,[lu_z0,Lu_z1,Lu_somerville_z0,Lu_lu_z0]
   endif
   ;Vincenzo2018
   if keyword_set(vincenzo18) then begin
      vincenzo = {feh:[-1.0,-0.75,-0.49,-0.40,-0.32,-0.26,-0.20,-0.12,-0.06,0.01,$
                       0.07,0.13,0.19,0.25,0.32,0.40,0.64],$
                  mgfe:[0.35,0.33,0.30,0.28,0.26,0.25,0.23,0.21,0.19,$
                        0.17,0.15,0.13,0.11,0.09,0.07,0.047,-0.04]}
      return, vincenzo
   endif
   ;zahid et al. 2012
   if keyword_set(zahid12) then begin
      readcol,'zahid12_eta.txt',mass_low,eta_low,mass_high,eta_high
      zahid = {mass_low:mass_low,eta_low:eta_low,mass_high:mass_high,eta_high:eta_high} 
      return, zahid
   endif

   if keyword_set(nelson19) then begin
      readcol,'nelson19_eta.txt',massv0,etav0,massv150,etav150,massv250,etav250
      return, {massv0:massv0,etav0:10.^(etav0),massv150:massv150,etav150:10.^(etav150),$
               massv250:massv250,etav250:10.^(etav250)}      
   endif
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451,badsn=badsn
   if ~keyword_set(badsn) then begin  
      wonsdss = where(sciall.zfit lt 0.15 and sciall.goodsn eq 1,cwonsdss)
      woncl0024 = where(sciall.zfit gt 0.35 and sciall.zfit lt 0.42 and sciall.goodsn eq 1,cwoncl0024)
      wonms0451 = where(sciall.zfit gt 0.48 and sciall.goodsn eq 1,cwonms0451)
   endif else begin
      wonsdss = where(sciall.zfit lt 0.15 and sciall.goodsn eq 0,cwonsdss)
      woncl0024 = where(sciall.zfit gt 0.35 and sciall.zfit lt 0.42 and sciall.goodsn eq 0,cwoncl0024)
      wonms0451 = where(sciall.zfit gt 0.48 and sciall.goodsn eq 0,cwonms0451)
   endelse

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro prepdata,sdssfile=sdssfile,cl0024file=cl0024file,ms0451file=ms0451file,$
             sdssprob=sdssprob,cl0024prob=cl0024prob,ms0451prob=ms0451prob,$
             outfile=outfile,plot=plot

   scitemp = {objname:'',ra:99.,dec:99.,snfit:99.,zspec:99.,zfit:99.,z:99.,age:99.,$
              ageerr:99.,ageupper:99.,agelower:99.,feh:99.,feherr:99.,$
              fehupper:99.,fehlower:99.,alphafe:99.,alphafeerr:99.,$
              alphafeupper:99.,alphafelower:99.,logmstar:99.,oiiew:99.,$
              oiiewerr:99.,vdisp:99.,vdisperr:99.,pfit:fltarr(13),pfiterr:fltarr(13),$
              ah:99.,ahupper:99.,ahlower:99.,goodsn:1b,ageform:99.} 

   ;GETTING SCI DATA
   sciall=[]
   goodspecall=[]
   proball = []
 
   ;SDSS
   scinow = mrdfits(sdssfile,1)
   probnow = mrdfits(sdssprob,1)
   ;sample selection
   good = where(scinow.haew gt -1. and scinow.logmstar gt 0. and scinow.good eq 1,cgood)
   scinow = scinow(good)
   scisdss = scinow
   probnow = probnow(good)
   goodfit = where(scinow.vdisp ge 0., cgoodfit,complement=badfit,ncomplement=cbadfit)
   goodspec = bytarr(cgood)
   goodspec(goodfit) = 1

   ;fix when the mass is wrong
   badmass = where(scinow.logmstar-scinow.logm50_ag06 lt -0.5)
   scinow(badmass).logmstar = scinow(badmass).logm50_ag06

   str = replicate(scitemp,cgood)
   struct_assign,scinow,str
   sciall = [sciall,str]
   proball = [proball,probnow]
   goodspecall = [goodspecall,goodspec]

   ;cl0024
   scinow = mrdfits(cl0024file,1)
   probnow = mrdfits(cl0024prob,1)

   ;sample selection
   good = where(scinow.good eq 1 and scinow.oiiew gt -5 and scinow.snfit gt 2. and scinow.logmstar gt 6.,cgood)
   scinow = scinow(good)
   scicl0024 = scinow
   probnow = probnow(good)
   goodfit = where(scinow.snfit gt 6.83, cgoodfit, complement=badfit,ncomplement=cbadfit)
   ;goodfit = where(scinow.snfit gt 6., cgoodfit, complement=badfit,ncomplement=cbadfit)
   ;change to gt 6 to reduce the s/n criteria to ~9 per rest Angstrom so that the sample size is similar to Leethochawalit2018
   goodspec = bytarr(cgood)
   goodspec(goodfit) = 1
   
   str = replicate(scitemp,cgood)
   struct_assign,scinow,str
   sciall = [sciall,str]
   proball = [proball,probnow]
   goodspecall = [goodspecall,goodspec]

   ;ms0451
   scinow = mrdfits(ms0451file,1)
   probnow = mrdfits(ms0451prob,1)
   ;sample selection
   good = where(scinow.good eq 1 and scinow.oiiew gt -5. and scinow.snfit gt 2 and scinow.logmstar gt 6.,cgood)
   scinow = scinow(good)
   scims0451=scinow
   probnow = probnow(good)
   goodfit = where(scinow.snfit gt 6.5,cgoodfit,complement=badfit,ncomplement=cbadfit)
   goodspec = bytarr(cgood)
   goodspec(goodfit) = 1

   str = replicate(scitemp,cgood) 
   struct_assign,scinow,str  
   sciall = [sciall,str]
   proball = [proball,probnow]
   goodspecall= [goodspecall,goodspec]

   sciall.goodsn  = goodspecall
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;MILES Library abundance (from Conroy17)
   miles = {feh:[-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],$
            ofe:[0.6,0.5,0.5,0.4,0.3,0.2,0.2,0.1,0.0,0.0],$
            mgfe:[0.4,0.4,0.4,0.4,0.34,0.22,0.14,0.11,0.05,0.04],$
            cafe:[0.32,0.30,0.28,0.26,0.26,0.17,0.12,0.06,0.00,0.00]}

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   ;fix the library abundance
    tofix = interpol(miles.mgfe,miles.feh,sciall.feh)
    sciall.alphafe = sciall.alphafe+tofix
    sciall.alphafeupper = sciall.alphafeupper+tofix
    sciall.alphafelower = sciall.alphafelower+tofix
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   sciall.ageform = (galage(sciall.zfit,1000.)/1.e9-sciall.age)>0. ;age of universe when it was formed
   sciall.ah = sciall.feh+sciall.alphafe
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   outlier = where(sciall.ah lt -0.3 and sciall.logmstar gt 10.1 and sciall.logmstar lt 10.4 $
                   and sciall.goodsn eq 1 and sciall.zfit gt 0.3 and sciall.zfit lt 0.5,coutlier)
   print,'excluding', coutlier, ' outliers from cl0024'
   if coutlier gt 0 then sciall(outlier).goodsn=0
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ngals = n_elements(sciall)
   dfeharr   = fltarr(51,ngals)
   dagearr   = fltarr(51,ngals)
   dalphaarr = fltarr(41,ngals)
   for nn=0,ngals-1 do begin
     dfeharr[*,nn] = proball[nn].feharr-proball[nn].feh50
     dagearr[*,nn] = proball[nn].agearr-proball[nn].age50
     dalphaarr[*,nn]=proball[nn].alphaarr-proball[nn].alpha50      
   endfor
   ;calculate probability of [Mg/H] from probability cubei
   probfehage = fltarr(51,51,ngals)
   probfehalpha = fltarr(51,41,ngals)
   cumprobfehalpha = fltarr(51,41,ngals)
   ndah = 51
   daharr = fltarr(ndah,ngals)
   probdah = fltarr(ndah,ngals)
   for nn=0,ngals-1 do begin 
      dfeh = dfeharr[1,nn]-dfeharr[0,nn]
      dalpha = dalphaarr[1,nn]-dalphaarr[0,nn]
      griddfeh = fltarr(51,41)
      griddalpha = fltarr(51,41)
      for ii=0,50 do for jj=0,40 do begin
         griddfeh[ii,jj] = dfeharr[ii,nn]
         griddalpha[ii,jj] = dalphaarr[jj,nn] 
         probfehalpha[ii,jj,nn] = int_tabulated(dagearr[*,nn],proball[nn].cubeprob[ii,*,jj])
         cumprobfehalpha[ii,jj,nn] = total(probfehalpha[*,*,nn])*dfeh*dalpha
      endfor 
      ;print, int_tabulated_2d(gridfeh,gridalpha,probfehalpha[*,*,nn]),cumprobfehalpha[50,40,nn] 
      probfehalpha[*,*,nn]=probfehalpha[*,*,nn]/int_tabulated_2d(griddfeh,griddalpha,probfehalpha[*,*,nn])
      cumprobfehalpha[*,*,nn] = cumprobfehalpha[*,*,nn]/cumprobfehalpha[50,40,nn]
      griddah = griddfeh+griddalpha
      order = sort(cumprobfehalpha[*,*,nn])
      curcumprob = (cumprobfehalpha[*,*,nn])[order]
      randompop = randomu(seed,5000)
      loc = value_locate(curcumprob,randompop)
      loc = order(loc)
     
      mindah = min(griddah(where(probfehalpha[*,*,nn] gt 0.)))
      maxdah = max(griddah(where(probfehalpha[*,*,nn] gt 0.)))
      curprobdah = histogram(griddah(loc),locations=curdaharr,min=mindah,max=maxdah,nbins=ndah)
      curprobdah = curprobdah/int_tabulated(curdaharr,float(curprobdah))
      probdah[*,nn] = curprobdah
      daharr[*,nn] = curdaharr

      curcumprobdah = fltarr(ndah)
      for jj=1,ndah-1 do curcumprobdah(jj) = int_tabulated(curdaharr[0:jj],curprobdah[0:jj])
      sciall[nn].ahupper = sciall[nn].ah+interpol(curdaharr,curcumprobdah,0.84)
      sciall[nn].ahlower = sciall[nn].ah+interpol(curdaharr,curcumprobdah,0.16)

      ;plot
      if keyword_set(plot) then begin
         set_plot,'x'
         !p.multi=[0,1,1]
          
         if nn eq 0 then cgdisplay, 1000, 600
         image = probfehalpha[*,*,nn]
         minValue = Floor(Min(image))
         maxValue = Ceil(Max(image))
         nLevels = 10
         xtitle = 'delta [Fe/H]'
         ytitle = 'delta [Mg/Fe]'
         position =   [0.1, 0.1, 0.5, 0.5]
         cgLoadCT, 33, CLIP=[30,255]
         cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
             /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
             XRange=minmax(dfeharr[*,nn]), YRange=minmax(dalphaarr[*,nn]), /Keep_Aspect
         oplot,!x.crange,0.-!x.crange,color=fsc_color('red'),linestyle=2
         oplot,!x.crange,0.5-!x.crange,color=fsc_color('red'),linestyle=2
         oplot,!x.crange,-0.5-!x.crange,color=fsc_color('red'),linestyle=2
         oplot,!x.crange,-1.-!x.crange,color=fsc_color('red'),linestyle=2
        
         cgplot,dfeharr[*,nn],proball[nn].probfeh,/noerase,position=[0.15,0.5,0.45,0.85],$
             xrange=minmax(dfeharr[*,nn]),color=fsc_color('black'),xstyle=5,$
             ytitle='prob([Fe/H])'
         axis,xaxis=1,xrange=minmax(dfeharr[*,nn]),color=fsc_color('black'),xstyle=1
         cgplot,proball[nn].probalpha,dalphaarr[*,nn],/noerase,position=[0.45,0.1,0.7,0.5],$
             yrange=minmax(dalphaarr[*,nn]),ystyle=1,color=fsc_color('black'),$
             xtickformat='(A1)',ytickformat='(A1)'
         axis,xaxis=1,xtitle='prob([Mg/Fe])',color=fsc_color('black')
         axis,yaxis=1,ystyle=1,yrange=minmax(dalphaarr[*,nn]),color=fsc_Color('black')
;         plothist,griddah(loc),xtitle='delta [Mg/H]',position=[0.65,0.6,0.95,0.95],/noerase
         cgplot,curdaharr,curprobdah,xtitle='delta [Mg/H]',ytitle='prob[Mg/H]',$
               position=[0.65,0.6,0.95,0.95],/noerase
         xyouts,0.2,0.95,sciall[nn].objname+' '+string(sciall[nn].snfit,format='(f5.2)'),$
                color=fsc_color('black'),/normal
         c=''
         ;if nn mod 30 eq 0 then read,c else wait,1
         wait,1
      endif

      ;sneakily tag along calculation of prob2d of age and feh
      dfeh = dfeharr[1,nn]-dfeharr[0,nn]
      dage = dagearr[1,nn]-dagearr[0,nn]
      griddfeh = fltarr(51,51) 
      griddage = fltarr(51,51)
      for ii=0,50 do for jj=0,50 do begin
         griddfeh[ii,jj] = dfeharr[ii,nn]
         griddage[ii,jj] = dagearr[jj,nn]
         probfehage[ii,jj,nn] = int_tabulated(dalphaarr[*,nn],proball[nn].cubeprob[ii,jj,*])
      endfor
;      print, int_tabulated_2d(griddfeh,griddage,probfehage[*,*,nn])
      probfehage[*,*,nn]=probfehage[*,*,nn]/int_tabulated_2d(griddfeh,griddage,probfehage[*,*,nn])
   endfor ;nn
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;save values
   proball_new = replicate({objname:'',feharr:fltarr(51),probfeh:fltarr(51),agearr:fltarr(51),$
                            probage:fltarr(51),alphaarr:fltarr(41),probalpha:fltarr(41),$
                            cubeprob:dblarr(51,51,41),age50:0.,feh50:0.,alpha50:0.,$
                            dfeharr:fltarr(51),dagearr:fltarr(51),dalphaarr:fltarr(41),$
                            daharr:fltarr(ndah),probdah:fltarr(ndah),probfehage:fltarr(51,51)},$
                            ngals)
   struct_Assign, proball,proball_new
   proball_new.daharr = daharr
   proball_new.probdah = probdah
   proball_new.dfeharr = dfeharr
   proball_new.dagearr = dagearr
   proball_new.dalphaarr = dalphaarr
   proball_new.probfehage = probfehage
   proball = proball_new

   mwrfits, sciall, outfile,/create,/silent
   mwrfits, proball,outfile,/silent
   mwrfits, scisdss,outfile,/silent
   mwrfits,scicl0024,outfile,/silent
   mwrfits,scims0451,outfile,/silent
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro linearfit_mzr,sciall,proball,prefix=prefix
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                 fileancova,deviationparam_file
   won = where(sciall.goodsn eq 1,cwon)
   
   allmzrpar = linfitmc(sciall[won].logmstar-10.,sciall[won].feh,proball[won].dfeharr,proball[won].probfeh,$
                        allmzrerr,chisq=chiallmzr);,/showplot)
   allmarpar = linfitmc(sciall[won].logmstar-10.,sciall[won].ah,proball[won].daharr,proball[won].probdah,$
                        allmarerr,chisq=chiallmar);,/showplot)

   separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451

   sdssmzrpar = linfitmc(sciall[wonsdss].logmstar-10.,sciall[wonsdss].feh,proball[wonsdss].dfeharr,$
                        proball[wonsdss].probfeh,sdssmzrerr,chisq=chisdssmzr);,/showplot)
   sdssmarpar = linfitmc(sciall[wonsdss].logmstar-10.,sciall[wonsdss].ah,proball[wonsdss].daharr,$
                        proball[wonsdss].probdah,sdssmarerr,chisq=chisdssmar);,/showplot)
   cl0024mzrpar = linfitmc(sciall[woncl0024].logmstar-10.,sciall[woncl0024].feh,proball[woncl0024].dfeharr,$
                        proball[woncl0024].probfeh,cl0024mzrerr,chisq=chicl0024mzr);,/showplot)
   cl0024marpar = linfitmc(sciall[woncl0024].logmstar-10.,sciall[woncl0024].ah,proball[woncl0024].daharr,$
                       proball[woncl0024].probdah,cl0024marerr,chisq=chicl0024mar);,/showplot)
   ms0451mzrpar = linfitmc(sciall[wonms0451].logmstar-10.,sciall[wonms0451].feh,proball[wonms0451].dfeharr,$
                       proball[wonms0451].probfeh,ms0451mzrerr,chisq=chims0451mzr);,/showplot)
   ms0451marpar = linfitmc(sciall[wonms0451].logmstar-10.,sciall[wonms0451].ah,proball[wonms0451].daharr,$
                       proball[wonms0451].probdah,ms0451marerr,chisq=chims0451mar);,/showplot)
   chiallmzr = chiallmzr/cwon
   chiallmar = chiallmar/cwon
   chisdssmzr = chisdssmzr/cwonsdss
   chisdssmar = chisdssmar/cwonsdss
   chicl0024mzr = chicl0024mzr/cwoncl0024
   chicl0024mar = chicl0024mar/cwoncl0024
   chims0451mzr = chims0451mzr/cwonms0451
   chims0451mar = chims0451mar/cwonms0451
   save,allmzrpar,allmarpar,sdssmzrpar,sdssmarpar,cl0024mzrpar,cl0024marpar,ms0451mzrpar,ms0451marpar,$
        allmzrerr,allmarerr,sdssmzrerr,sdssmarerr,cl0024mzrerr,cl0024marerr,ms0451mzrerr,ms0451marerr,$
        chiallmzr,chiallmar,chisdssmzr,chisdssmar,chicl0024mzr,chicl0024mar,chims0451mzr,chims0451mar,$
        filename=linfitparam_file
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro test_evolution_regression,sciall,proball,prefix=prefix
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                 fileancova,deviationparam_file
   ;read the linfit param
   restore, linfitparam_file

   ;two-sample t-test for slope and intercept
   tslope_sdcl_feh = (sdssmzrpar[1]-cl0024mzrpar[1])/sqrt(sdssmzrerr[1]^2+cl0024mzrerr[1]^2)
   tslope_sdms_feh = (sdssmzrpar[1]-ms0451mzrpar[1])/sqrt(sdssmzrerr[1]^2+ms0451mzrerr[1]^2)
   tslope_clms_feh = (cl0024mzrpar[1]-ms0451mzrpar[1])/sqrt(cl0024mzrerr[1]^2+ms0451mzrerr[1]^2)
   tslope_sdcl_mgh = (sdssmarpar[1]-cl0024marpar[1])/sqrt(sdssmarerr[1]^2+cl0024marerr[1]^2)
   tslope_sdms_mgh = (sdssmarpar[1]-ms0451marpar[1])/sqrt(sdssmarerr[1]^2+ms0451marerr[1]^2)
   tslope_clms_mgh = (cl0024marpar[1]-ms0451marpar[1])/sqrt(cl0024marerr[1]^2+ms0451marerr[1]^2)

   ;use t values from here https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
   ;https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   ;if cannot reject the null hypothesis then proceed to find common slope
   ;http://www.biostathandbook.com/ancova.html  
   ;http://www.stat.ucla.edu/~cochran/stat10/winter/lectures/lect21.html
   ;find the weighted slope
   fehcommonslope = (ms0451mzrpar[1]/ms0451mzrerr[1]^2+cl0024mzrpar[1]/cl0024mzrerr[1]^2+$
                      sdssmzrpar[1]/sdssmzrerr[1]^2)/(1./ms0451mzrerr[1]^2+1./cl0024mzrerr[1]^2+$
                      1./sdssmzrerr[1]^2)
   fehcommonslope_err = 1./sqrt(1./ms0451mzrerr[1]^2+1./cl0024mzrerr[1]^2+1./sdssmzrerr[1]^2)
   fehcommonslope_stdev = sqrt(((ms0451mzrpar[1]-fehcommonslope)^2/ms0451mzrerr[1]^2+$
                                (cl0024mzrpar[1]-fehcommonslope)^2/cl0024mzrerr[1]^2+$
                                (sdssmzrpar[1]-fehcommonslope)^2/sdssmzrerr[1]^2)/$
                               (1./ms0451mzrerr[1]^2+1./cl0024mzrerr[1]^2+1./sdssmzrerr[1]^2)/3.)

   mghcommonslope = (ms0451marpar[1]/ms0451marerr[1]^2+cl0024marpar[1]/cl0024marerr[1]^2+$
                      sdssmarpar[1]/sdssmarerr[1]^2)/(1./ms0451marerr[1]^2+1./cl0024marerr[1]^2+$
                      1./sdssmarerr[1]^2)
   mghcommonslope_err = 1./sqrt(1./ms0451marerr[1]^2+1./cl0024marerr[1]^2+1./sdssmarerr[1]^2)
   mghcommonslope_stdev = sqrt(((ms0451marpar[1]-mghcommonslope)^2/ms0451marerr[1]^2+$
                                (cl0024marpar[1]-mghcommonslope)^2/cl0024marerr[1]^2+$
                                (sdssmarpar[1]-mghcommonslope)^2/sdssmarerr[1]^2)/$
                               (1./ms0451marerr[1]^2+1./cl0024marerr[1]^2+1./sdssmarerr[1]^2)/3.)

   ; Find difference in the intercepts when fitting with a fixed slope
   separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451

   sdssmzrpar_fixslope = [sdssmzrpar[0],fehcommonslope]
   sdssmzryfit_fixslope = linfit_fixedslope_mc(sciall(wonsdss).logmstar-10.,sciall(wonsdss).feh,$
                               proball(wonsdss).dfeharr,proball(wonsdss).probfeh,sdssmzrpar_fixslope,$
                               sdssmzrerr_fixslope,chisq=sdssmzrchi_fixslope)
   sdssmarpar_fixslope = [sdssmarpar[0],mghcommonslope]
   sdssmaryfit_fixslope = linfit_fixedslope_mc(sciall(wonsdss).logmstar-10.,sciall(wonsdss).ah,$
                               proball(wonsdss).daharr,proball(wonsdss).probdah,sdssmarpar_fixslope,$
                               sdssmarerr_fixslope,chisq=sdssmarchi_fixslope)

   cl0024mzrpar_fixslope = [cl0024mzrpar[0],fehcommonslope]
   cl0024mzryfit_fixslope = linfit_fixedslope_mc(sciall(woncl0024).logmstar-10.,sciall(woncl0024).feh,$
                               proball(woncl0024).dfeharr,proball(woncl0024).probfeh,cl0024mzrpar_fixslope,$
                               cl0024mzrerr_fixslope,chisq=cl0024mzrchi_fixslope)
   cl0024marpar_fixslope = [cl0024marpar[0],mghcommonslope]
   cl0024maryfit_fixslope = linfit_fixedslope_mc(sciall(woncl0024).logmstar-10.,sciall(woncl0024).ah,$
                               proball(woncl0024).daharr,proball(woncl0024).probdah,cl0024marpar_fixslope,$
                               cl0024marerr_fixslope,chisq=cl0024marchi_fixslope)

   ms0451mzrpar_fixslope = [ms0451mzrpar[0],fehcommonslope]
   ms0451mzryfit_fixslope = linfit_fixedslope_mc(sciall(wonms0451).logmstar-10.,sciall(wonms0451).feh,$
                               proball(wonms0451).dfeharr,proball(wonms0451).probfeh,ms0451mzrpar_fixslope,$
                               ms0451mzrerr_fixslope,chisq=ms0451mzrchi_fixslope)
   ms0451marpar_fixslope = [ms0451marpar[0],mghcommonslope]
   ms0451maryfit_fixslope = linfit_fixedslope_mc(sciall(wonms0451).logmstar-10.,sciall(wonms0451).ah,$
                               proball(wonms0451).daharr,proball(wonms0451).probdah,ms0451marpar_fixslope,$
                               ms0451marerr_fixslope,chisq=ms0451marchi_fixslope)
   sdssmzrchi_fixslope = sdssmzrchi_fixslope/cwonsdss
   sdssmarchi_fixslope = sdssmarchi_fixslope/cwonsdss
   cl0024mzrchi_fixslope = cl0024mzrchi_fixslope/cwoncl0024
   cl0024marchi_fixslope = cl0024marchi_fixslope/cwoncl0024
   ms0451mzrchi_fixslope = ms0451mzrchi_fixslope/cwonms0451
   ms0451marchi_fixslope = ms0451marchi_fixslope/cwonms0451
   ;z test comparing two means
   tconst_sdcl_feh = ((sdssmzrpar_fixslope[0]-cl0024mzrpar_fixslope[0])-0.)/$
                  sqrt(sdssmzrerr_fixslope[0]^2+cl0024mzrerr_fixslope[0]^2)
   tconst_sdcl_mgh = ((sdssmarpar_fixslope[0]-cl0024marpar_fixslope[0])-0.)/$
                  sqrt(sdssmarerr_fixslope[0]^2+cl0024marerr_fixslope[0]^2)
   tconst_sdms_feh = ((sdssmzrpar_fixslope[0]-ms0451mzrpar_fixslope[0])-0.)/$
                  sqrt(sdssmzrerr_fixslope[0]^2+ms0451mzrerr_fixslope[0]^2)
   tconst_sdms_mgh = ((sdssmarpar_fixslope[0]-ms0451marpar_fixslope[0])-0.)/$
                  sqrt(sdssmarerr_fixslope[0]^2+ms0451marerr_fixslope[0]^2)
   tconst_clms_feh = ((cl0024mzrpar_fixslope[0]-ms0451mzrpar_fixslope[0])-0.)/$
                  sqrt(cl0024mzrerr_fixslope[0]^2+ms0451mzrerr_fixslope[0]^2)
   tconst_clms_mgh = ((cl0024marpar_fixslope[0]-ms0451marpar_fixslope[0])-0.)/$
                  sqrt(cl0024marerr_fixslope[0]^2+ms0451marerr_fixslope[0]^2)

   save, tconst_sdcl_feh,tconst_sdcl_mgh,tconst_sdms_feh,$
         tconst_sdms_mgh,tconst_clms_feh,tconst_clms_mgh,$
         tslope_sdcl_feh,tslope_sdcl_mgh,tslope_sdms_feh,$
         tslope_sdms_mgh,tslope_clms_feh,tslope_clms_mgh,$
         fehcommonslope,fehcommonslope_err,fehcommonslope_stdev,$
         mghcommonslope,mghcommonslope_err,mghcommonslope_stdev,$
         sdssmzrpar_fixslope,cl0024mzrpar_fixslope,ms0451mzrpar_fixslope,$
         sdssmzrerr_fixslope,cl0024mzrerr_fixslope,ms0451mzrerr_fixslope,$
         sdssmzrchi_fixslope,cl0024mzrchi_fixslope,ms0451mzrchi_fixslope,$
         sdssmarpar_fixslope,cl0024marpar_fixslope,ms0451marpar_fixslope,$
         sdssmarerr_fixslope,cl0024marerr_fixslope,ms0451marerr_fixslope,$
         sdssmarchi_fixslope,cl0024marchi_fixslope,ms0451marchi_fixslope,$
         filename=linfitparam_fixslope_file
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ancova,sciall,proball,fileancova=fileancova
   won = where(sciall.goodsn eq 1,cwon)
   ;fit linear regression with monte carlo technique 
   ;setting up all the neccesary parameters
   nobjs = cwon
   feh  = sciall[won].feh
   dfeh = proball[won].dfeharr
   probdfeh = proball[won].probfeh

   mgh = sciall[won].ah
   dmgh = proball[won].daharr
   probdmgh = proball[won].probdah

   mass = sciall[won].logmstar-10.

   redshift = sciall[won].zfit
   cl = fltarr(cwon)
   ms = fltarr(cwon)
   woncl = where(redshift gt 0.35 and redshift lt 0.42,cwoncl)
   wonms = where(redshift gt 0.48 and redshift lt 0.65,cwonms)
   cl(woncl) = 1
   ms(wonms) = 1

   m_cl = cl*mass
   m_ms = ms*mass

   ;first with the simple model (no interaction terms)
   x = [transpose(mass),transpose(cl),transpose(ms)]
   df1_mo1 = 3
   df2_mo1 = nobjs-df1_mo1-1
   mo1_pname = ['mass    ','cl0024  ','ms0451  '] 
   coeff_feh_mo1 = regress_mc(x,feh,dfeh,probdfeh,const=const_feh_mo1,sigmacoeff=sigma_feh_mo1,$
                      sigmaconst=consterr_feh_mo1,rsq=rsq_feh_mo1,stderror=stderror_feh_mo1,$
                      adj_rsq=adjrsq_feh_mo1)
   coeff_mgh_mo1 = regress_mc(x,mgh,dmgh,probdmgh,const=const_mgh_mo1,sigmacoeff=sigma_mgh_mo1,$
                      sigmaconst=consterr_mgh_mo1,rsq=rsq_mgh_mo1,stderror=stderror_mgh_mo1,$
                      adj_rsq=adjrsq_mgh_mo1)
   ;second with interaction terms
   x = [transpose(mass),transpose(cl),transpose(ms),transpose(m_cl),transpose(m_ms)]
   df1_mo2 = 5
   df2_mo2 = nobjs-df1_mo2-1
   mo2_pname = ['mass      ','cl0024    ','ms0451    ','mass*cl0024','mass*ms0451'] 
   coeff_feh_mo2 = regress_mc(x,feh,dfeh,probdfeh,const=const_feh_mo2,sigmacoeff=sigma_feh_mo2,$
                      sigmaconst=consterr_feh_mo2,rsq=rsq_feh_mo2,stderror=stderror_feh_mo2,$
                      adj_rsq=adjrsq_feh_mo2)
   coeff_mgh_mo2 = regress_mc(x,mgh,dmgh,probdmgh,const=const_mgh_mo2,sigmacoeff=sigma_mgh_mo2,$
                      sigmaconst=consterr_mgh_mo2,rsq=rsq_mgh_mo2,stderror=stderror_mgh_mo2,$
                      adj_rsq=adjrsq_mgh_mo2)
   ;for mgh, only: also do it with a super simple. no evolution terms
   x = [transpose(mass)]
   df1_mo3 = 1
   df2_mo3 = nobjs-df1_mo3-1
   mo3_pname = ['mass   ']
   coeff_mgh_mo3 = regress_mc(x,mgh,dmgh,probdmgh,const=const_mgh_mo3,sigmacoeff=sigma_mgh_mo3,$
                      sigmaconst=consterr_mgh_mo3,rsq=rsq_mgh_mo3,stderror=stderror_mgh_mo3,$
                      adj_rsq=adjrsq_mgh_mo3)
   ;save results
   save, coeff_feh_mo1,const_feh_mo1,sigma_feh_mo1,rsq_feh_mo1,stderror_feh_mo1,adjrsq_feh_mo1,$
         consterr_feh_mo1,$
         coeff_feh_mo2,const_feh_mo2,sigma_feh_mo2,rsq_feh_mo2,stderror_feh_mo2,adjrsq_feh_mo2,$
         consterr_feh_mo2,$
         coeff_mgh_mo1,const_mgh_mo1,sigma_mgh_mo1,rsq_mgh_mo1,stderror_mgh_mo1,adjrsq_mgh_mo1,$
         consterr_mgh_mo1,$
         coeff_mgh_mo2,const_mgh_mo2,sigma_mgh_mo2,rsq_mgh_mo2,stderror_mgh_mo2,adjrsq_mgh_mo2,$
         consterr_mgh_mo2,$
         coeff_mgh_mo3,const_mgh_mo3,sigma_mgh_mo3,rsq_mgh_mo3,stderror_mgh_mo3,adjrsq_mgh_mo3,$
         consterr_mgh_mo3,$
         df1_mo1,df2_mo1,df1_mo2,df2_mo2,df1_mo3,df2_mo3,mo1_pname,mo2_pname,mo3_pname,filename=fileancova
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro anova_massbin,sciall,proball,prefix=prefix,fileanova=fileanova
   minmass = 10.
   maxmass = 11.5
   mbin = 0.3 ;mass bin in logmass
   nbins = floor((maxmass-minmass)/mbin)
   outstr_feh = replicate({mbinmin:0d,mbinmax:0d,f_anova:0d,dof_ssb:0,dof_ssw:0,$
                           SSW:0d,SSB:0d,TSS:0d,metal_mean:[0.,0.,0.],metal_sigmas:[0.,0.,0.]},nbins)
   outstr_mgh = outstr_feh

   separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451
   sdss = sciall(wonsdss)
   ms0451 = sciall(wonms0451)
   cl0024 = sciall(woncl0024)

   probsd = proball(wonsdss)
   probms = proball(wonms0451)
   probcl = proball(woncl0024)
   
   for i=0,nbins-1 do begin
      iminmass = minmass+i*mbin
      imaxmass = minmass+(i+1.)*mbin
      goodsdss = where(sdss.logmstar ge iminmass and sdss.logmstar lt imaxmass, cgoodsdss)
      sdssnow = sdss(goodsdss)
      probsdnow = probsd(goodsdss)

      goodcl0024 = where(cl0024.logmstar ge iminmass and cl0024.logmstar lt imaxmass, cgoodcl0024)
      cl0024now = cl0024(goodcl0024)
      probclnow = probcl(goodcl0024)

      goodms0451 = where(ms0451.logmstar ge iminmass and ms0451.logmstar lt imaxmass, cgoodms0451)
      ms0451now = ms0451(goodms0451)
      probmsnow = probms(goodms0451)

      allnow = [sdssnow,cl0024now,ms0451now]
      proballnow = [probsdnow, probclnow, probmsnow]

      ;FEH ANOVA
      meanmc,sdssnow.feh,probsdnow.dfeharr,probsdnow.probfeh,sdss_meanfeh,fehsigmam,$
             fehsigmad,fehsigmas_sd
      meanmc,cl0024now.feh,probclnow.dfeharr,probclnow.probfeh,cl0024_meanfeh,fehsigmam,$
             fehsigmad,fehsigmas_cl
      meanmc,ms0451now.feh,probmsnow.dfeharr,probmsnow.probfeh,ms0451_meanfeh,fehsigmam,$
             fehsigmad,fehsigmas_ms
      meanmc,allnow.feh,proballnow.dfeharr,proballnow.probfeh,all_meanfeh,fehsigmam,$
             fehsigmad,fehsigmas_all

      SSW = total((sdssnow.feh-sdss_meanfeh)^2)+total((cl0024now.feh-cl0024_meanfeh)^2)+$
            total((ms0451now.feh-ms0451_meanfeh)^2)  ;sum of squares within group 
      TSS = total((allnow.feh-all_meanfeh)^2) ;total sum of squares
      SSB = (sdss_meanfeh-all_meanfeh)^2+(cl0024_meanfeh-all_meanfeh)^2+$
            (ms0451_meanfeh-all_meanfeh)^2 ;sum of squares between groups
      dof_ssb = 2 ;= ngroups-1
      dof_ssw = cgoodsdss+cgoodcl0024+cgoodms0451-3 ;total obs - ngroups
 
      f_anova = (ssb/dof_ssb)/(ssw/dof_ssw)
      outstr_feh[i].mbinmin = iminmass
      outstr_feh[i].mbinmax = imaxmass
      outstr_feh[i].f_anova = f_anova
      outstr_feh[i].dof_ssb = dof_ssb
      outstr_feh[i].dof_ssw = dof_ssw
      outstr_feh[i].SSW = SSW
      outstr_feh[i].SSB = SSB
      outstr_feh[i].TSS = TSS
      outstr_feh[i].metal_mean = [sdss_meanfeh,cl0024_meanfeh,ms0451_meanfeh]
      outstr_feh[i].metal_sigmas = [fehsigmas_sd,fehsigmas_cl,fehsigmas_ms]
      ;MgH ANOVA
      meanmc,sdssnow.ah,probsdnow.daharr,probsdnow.probdah,sdss_meanah,ahsigmam,$
             ahsigmad,ahsigmas_sd
      meanmc,cl0024now.ah,probclnow.daharr,probclnow.probdah,cl0024_meanah,ahsigmam,$
             ahsigmad,ahsigmas_cl
      meanmc,ms0451now.ah,probmsnow.daharr,probmsnow.probdah,ms0451_meanah,ahsigmam,$
             ahsigmad,ahsigmas_ms
      meanmc,allnow.ah,proballnow.daharr,proballnow.probdah,all_meanah,ahsigmam,$
             ahsigmad,ahsigmas

      SSW = total((sdssnow.ah-sdss_meanah)^2)+total((cl0024now.ah-cl0024_meanah)^2)+$
            total((ms0451now.ah-ms0451_meanah)^2)  ;sum of squares within group 
      TSS = total((allnow.ah-all_meanah)^2) ;total sum of squares
      SSB = (sdss_meanah-all_meanah)^2+(cl0024_meanah-all_meanah)^2+$
            (ms0451_meanah-all_meanah)^2 ;sum of squares between groups

      f_anova = (ssb/dof_ssb)/(ssw/dof_ssw)

      outstr_mgh[i].mbinmin = iminmass
      outstr_mgh[i].mbinmax = imaxmass
      outstr_mgh[i].f_anova = f_anova
      outstr_mgh[i].dof_ssb = dof_ssb
      outstr_mgh[i].dof_ssw = dof_ssw
      outstr_mgh[i].SSW = SSW
      outstr_mgh[i].SSB = SSB
      outstr_mgh[i].TSS = TSS
      outstr_mgh[i].metal_mean = [sdss_meanah,cl0024_meanah,ms0451_meanah]
      outstr_mgh[i].metal_sigmas = [ahsigmas_sd,ahsigmas_cl,ahsigmas_ms]
   endfor
   mwrfits,outstr_feh,fileanova,/create,/silent
   mwrfits,outstr_mgh,fileanova,/silent
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro average_mzr,sciall,proball,prefix=prefix
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                 fileancova,deviationparam_file
   ;set values for each sample set (clusters)
   zmin = [0.,0.35,0.48]
   zmax = [0.15,0.42,0.6]
   ngalperbin = [10,8,8]
 
   for jj=0,2 do begin
      nperbin = ngalperbin[jj]
      won = where(sciall.zfit lt zmax[jj] and sciall.zfit gt zmin[jj] and sciall.goodsn eq 1)
      scinow = sciall(won)
      probnow = proball(won)

      order = sort(scinow.logmstar)
      scinow = scinow(order)
      probnow = probnow(order)

      nbins_now = floor(n_elements(scinow)/float(nperbin))
      locbin = findgen(nbins_now)*nperbin
      massbins_now = scinow(locbin).logmstar
      massbins_now = [massbins_now,11.7]

      ave_mass_now    = fltarr(nbins_now)
      ave_mass_dev_now= fltarr(nbins_now)
      ave_feh_now     = fltarr(nbins_now)
      ave_feh_dev_now = fltarr(nbins_now)
      ave_mgh_now     = fltarr(nbins_now)
      ave_mgh_dev_now = fltarr(nbins_now)

      for i=0,nbins_now-1 do begin
         msel = where(scinow.logmstar ge massbins_now[i] and scinow.logmstar lt massbins_now[i+1],cmsel)
         if cmsel lt nperbin then stop,'STOPPED: something is wrong about mass average'
         meanmc,scinow(msel).feh,probnow(msel).dfeharr,probnow(msel).probfeh,fehmean,fehsigmam,$
             fehsigmad,fehsigmas
         meanmc,scinow(msel).ah,probnow(msel).daharr,probnow(msel).probdah,mghmean,mghsigmam,$
             mghsigmad,mghsigmas
         ave_feh_now(i) = fehmean
         ave_feh_dev_now(i) = fehsigmas
         ave_mgh_now(i) = mghmean
         ave_mgh_dev_now(i) =  mghsigmas
         ave_mass_now(i) = mean(scinow(msel).logmstar)
         ave_mass_dev_now(i) = stdev(scinow(msel).logmstar)
      endfor
      bndry_mass_now = massbins_now
      hifeh_now = interpol(ave_feh_now+ave_feh_dev_now,ave_mass_now,bndry_mass_now)
      lofeh_now = interpol(ave_feh_now-ave_feh_dev_now,ave_mass_now,bndry_mass_now)
      himgh_now = interpol(ave_mgh_now+ave_mgh_dev_now,ave_mass_now,bndry_mass_now)
      lomgh_now = interpol(ave_mgh_now-ave_mgh_dev_now,ave_mass_now,bndry_mass_now)
      outstr = {galperbin:nperbin,massbins:massbins_now,ave_mass:ave_mass_now,$
                ave_mass_dev:ave_mass_dev_now,ave_feh:ave_feh_now,$
                ave_feh_dev:ave_feh_dev_now,hifeh:hifeh_now,lofeh:lofeh_now,$
                ave_mgh:ave_mgh_now,ave_mgh_dev:ave_mgh_dev_now,himgh:himgh_now,lomgh:lomgh_now}
      if jj eq 0 then sdss_avestr = outstr
      if jj eq 1 then cl0024_avestr = outstr
      if jj eq 2 then ms0451_avestr = outstr
   endfor
   save, sdss_avestr,cl0024_avestr,ms0451_avestr,filename=aveparam_file
end

pro deviation_fromz0line,sciall,proball,fehcoeff,mghcoeff,fileout=fileout
;calculate deviation as a function of formation age
   ;calculate FEH IDEAL
   fehideal = fltarr(n_Elements(sciall))
   separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451 
   fehideal(wonsdss) = fehcoeff[0]+(sciall(wonsdss).logmstar-10.)*fehcoeff[1]
   fehideal(woncl0024) = fehcoeff[0]+(sciall(woncl0024).logmstar-10.)*(fehcoeff[1]+fehcoeff[2])
   fehideal(wonms0451) = fehcoeff[0]+(sciall(wonms0451).logmstar-10.)*(fehcoeff[1]+fehcoeff[3])

   good = where(sciall.ageform gt 0.5 and sciall.goodsn eq 1,cgood,complement=bad,ncomplement=cbad)
   print,'number of outliers (feh deviation) is', cbad
   x=sciall(good).ageform
   y=sciall(good).feh-fehideal(good)
   dx = 0.5*(sciall(good).ageupper-sciall(good).agelower)
   dy = 0.5*(sciall(good).fehupper-sciall(good).fehlower)
   dxarr = -1.*proball(good).dagearr ;since x axis is ageuni-agegal
   dyarr = proball(good).dfeharr
   probdxdyarr = transpose(proball(good).probfehage,[1,0,2]) ;age feh nobjs
  
   feh_dev_par = linfitprobgrid_general(x,y,dxarr,dyarr,probdxdyarr,feh_dev_par_all,$
                   outfitfile='all_prob_deltafeh_age.fits',outepsfile='all_plot_deltafeh_age.eps',$
                   fitrange=[0,11.5],/plot,initial_guess=[-0.4,0.06],stepsize=[0.0005,0.0005],$
                   nstep=3000)

   mc_feh_dev_par = linfitprobgrid_mc(x,y,dxarr,dyarr,probdxdyarr,mc_feh_dev_par_err)

   simple_feh_dev_par = linfit(x,y,sigma=feh_dev_par_err)

   ;MAR
   ahideal = mghcoeff[0]+(sciall.logmstar-10.)*mghcoeff[1]
   x=sciall.ageform
   y=sciall.ah-ahideal

   good = where(y gt -0.5 and x gt 0.5,cgood,complement=bad,ncomplement=cbad)
   print,'number of outliers (mgh deviation) is', cbad
   simple_mgh_dev_par = linfit(x(good),y(good),sigma=simple_mgh_dev_par_err)
   ;mc fit
   y1 = sciall(good).feh
   y2 = sciall(good).alphafe
   dxarr = -1.*proball(good).dagearr ;since x axis is ageuni-agegal
   dy1arr = proball(good).dfeharr
   dy2arr = proball(good).dalphaarr
   probcube = transpose(proball(good).cubeprob,[1,0,2,3])
   mc_mgh_dev_par = linfitprobgrid_mc_mgh(x(good),y1,y2,dxarr,dy1arr,dy2arr,probcube,mc_mgh_dev_par_err)

   save, feh_dev_par,feh_dev_par_all,mc_feh_dev_par,mc_feh_dev_par_err,simple_feh_dev_par,$
         feh_dev_par_err,mc_mgh_dev_par,mc_mgh_dev_par_err,simple_mgh_dev_par,$
         simple_mgh_dev_par_err,fehideal,ahideal,filename=fileout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro plot_deviation_fromz0line,sciall
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                  fileancova,deviationparam_file
common color_pallete, clustercolor, shadecolor, meancolor

   restore,deviationparam_file
   goodfit = where(sciall.goodsn eq 1)
   symsize = 1.5/(max(sciall(goodfit).snfit)-min(sciall(goodfit).snfit))*sciall.snfit+0.7
   psname='FEH_deviation_fromz0line.eps'

   xfitrange = [0,13]
   yfitrange = [-0.4,0.4]
;   reduceimage,'all_prob_deltafeh_age.fits',xfitrange,yfitrange,img,ximg,yimg
   plotposition=[0.15,0.18,0.95,0.87]
;   set_plot,'x'
;   loadct,60
;   cgdisplay,/pixmap
;   cgimage,img_sdss,position=plotposition
;   bgimage = cgsnapshot()
;   wdelete

   ntaletter = "150B
   deltaletter = "104B
   sunsym = sunsymbol()

   set_plot,'ps'
   !p.font=0
   !p.charsize=1.4
   device, filename = psname,xsize = 15,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;
;      loadct,54
;      good = where(finite(img),cgood,complement=bad)
;      img(bad) = min(img(good))
;      loadct,55
;      cgimage,img,position=plotposition
;      cgimage,img,transparent=50,alphabackgroundimage=bgimage,alphafgposition=plotposition
      plot,sciall.ageform,sciall.feh-fehideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle='!9'+string(deltaletter)+'!x[Fe/H]',xrange=xfitrange,xstyle=9,yrange=yfitrange,/nodata,position=plotposition,/noerase

      axis,xaxis=1,xticks=4,xtickv=[1.513,3.223, 5.747,8.422 ,12.161],xtickn=['4','2','1','0.5','0.1'],xtitle='z'

      for i=0,n_elements(sciall)-1 do begin
          if sciall[i].goodsn eq 1 and sciall[i].ageform gt 0.5 then begin
             if sciall[i].zspec lt 0.3 then colornow = 'skyblue'
             if sciall[i].zspec gt 0.3 and sciall[i].zspec lt 0.45 then colornow = 'goldenrod'
             if sciall[i].zspec gt 0.45 then colornow = 'crimson'
             cgplot,[sciall[i].ageform],[sciall[i].feh-fehideal[i]],psym=46,/overplot,symsize=symsize[i],color=colornow
          endif
      endfor

;      oplot,[0,14],feh_dev_par(0)+feh_dev_par(1)*[0,14],thick=1,linestyle=1
      oplot,[0,14],mc_feh_dev_par(0)+mc_feh_dev_par(1)*[0,14],thick=1,linestyle=2;,color=fsc_color('red') 
      ;typical uncertainties
      zmin =[0.,0.3,0.45]
      zmax =[0.3,0.45,0.7]
      yplot = [0.29]
      xplot = [1.7,5,8.]
      
      for i=0,2 do begin
         sel = where(sciall.goodsn eq 1 and sciall.ageform gt 0.5 and $
                     sciall.zspec gt zmin[i] and sciall.zspec le zmax[i])
         upperlim = [median(sciall(sel).fehupper-sciall(sel).feh)]
         lowerlim = [median(sciall(sel).feh-sciall(sel).fehlower)]
         leftlim  = [median(sciall(sel).ageupper-sciall(sel).age)]
         rightlim = [median(sciall(sel).age-sciall(sel).agelower)]
         cgplot,[xplot[i]],yplot,err_xhigh=rightlim,err_xlow=leftlim,$
                err_yhigh=upperlim,err_ylow=lowerlim,psym=46,$
                color=clustercolor[i],/overplot,symsize=1.2,err_thick=2
      endfor
      ;;;;;;;;;;;;;;;;;;;;;

      xyouts,1,-0.38,'older galaxies',charsize=0.8
      xyouts,9.4,-0.38,'younger galaxies',charsize=0.8
      arrow,0.9,-0.37,0.3,-0.37,/data
      arrow,12.1,-0.37,12.7,-0.37,/data
   device,/close

   set_plot,'ps'
   psname = 'mgh_deviation_fromz0line.eps'
   device, filename = psname,xsize = 15,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;
      yfitrange=[-0.5,0.8]
      plot,sciall.ageform,sciall.ah-ahideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle='!9'+string(deltaletter)+'!x[Mg/H]',xrange=xfitrange,xstyle=9,yrange=yfitrange,/nodata,position=plotposition,/noerase,ystyle=1

      axis,xaxis=1,xticks=4,xtickv=[1.513,3.223, 5.747,8.422 ,12.161],xtickn=['4','2','1','0.5','0.1'],xtitle='z'
      for i=0,n_elements(sciall)-1 do begin
          if sciall[i].goodsn eq 1 and sciall[i].ageform gt 0.5 then begin
             if sciall[i].zspec lt 0.3 then colornow = 'skyblue'
             if sciall[i].zspec gt 0.3 and sciall[i].zspec lt 0.45 then colornow = 'goldenrod'
             if sciall[i].zspec gt 0.45 then colornow = 'crimson'
             cgplot,[sciall[i].ageform],[sciall[i].ah-ahideal[i]],psym=46,/overplot,symsize=symsize[i],color=colornow
          endif
      endfor
;      bad = where(sciall.goodsn eq 0)
;      cgplot,sciall(bad).ageform,sciall(ms0451).ah-ahideal(ms0451),psym=16,/overplot,symsize=1,color='gray'

 ;     oplot,[0,14],simple_mgh_dev_par(0)+simple_mgh_dev_par(1)*[0,14],thick=1,linestyle=2
      oplot,[0,14],mc_mgh_dev_par(0)-0.1+mc_mgh_dev_par(1)*[0,14],thick=1,linestyle=2
      ;oplot,[0,14],ah_dev_par_sdss(0)+ah_dev_par_sdss(1)*[0,14],thick=2,linestyle=2,color=fsc_color('blu5')
      ;oplot,[0,14],ah_dev_par_cl(0)+ah_dev_par_cl(1)*[0,14],thick=2,linestyle=2,color=fsc_color('org6')
;typical uncertainties
      zmin =[0.,0.3,0.45]
      zmax =[0.3,0.45,0.7]
      yplot = [0.65]
      xplot = [1.7,5,8.]

      for i=0,2 do begin
         sel = where(sciall.goodsn eq 1 and sciall.ageform gt 0.5 and $
                     sciall.zspec gt zmin[i] and sciall.zspec le zmax[i])
         upperlim = [median(sciall(sel).ahupper-sciall(sel).ah)]
         lowerlim = [median(sciall(sel).ah-sciall(sel).ahlower)]
         leftlim  = [median(sciall(sel).ageupper-sciall(sel).age)]
         rightlim = [median(sciall(sel).age-sciall(sel).agelower)]
         cgplot,[xplot[i]],yplot,err_xhigh=rightlim,err_xlow=leftlim,$
                err_yhigh=upperlim,err_ylow=lowerlim,psym=46,$
                color=clustercolor[i],/overplot,symsize=1.2,err_thick=2
      endfor
      ;;;;;;;;;;;;;;;;;;;;;

      xyouts,1,-0.48,'older galaxies',charsize=0.8
      xyouts,9.3,-0.48,'younger galaxies',charsize=0.8
      arrow,0.9,-0.46,0.3,-0.46,/data
      arrow,12.1,-0.46,12.7,-0.46,/data
   device,/close
end


pro plot_masshist,sciall
common color_pallete, clustercolor, shadecolor, meancolor

   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()
   psname= 'mass_hist_allsample.eps'

   separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451

   device, filename = psname,xsize = 15,ysize = 10, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plothist,sciall(wonsdss).logmstar,xrange=[9.7,12],xtitle='log[M!L*!N/M!N'+sunsym+'!N]',/fill,$
           fcolor=fsc_color(zcol[0]),bin=0.3,xtickformat='(A1)',yrange=[0,40]
      plothist,sciall(wonms0451).logmstar,/overplot,color=fsc_color('grn7'),/fline,$
           fcolor=fsc_Color(zcol[2]),bin=0.3,forientation=135
      plothist,sciall(woncl0024).logmstar,/overplot,color=fsc_color('grn7'),/fline,$
           fcolor=fsc_Color(zcol[1]),forientation=45,bin=0.3
      axis,xaxis=0,xrange=[9.7,12],xstyle=1

   device,/close
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro plot_mzr,sciall,prefix=prefix
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                fileancova,deviationparam_file
common color_pallete, clustercolor, shadecolor, meancolor
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'

   restore, aveparam_file  
;     output: sdss_avestr,cl0024_avestr,ms0451_avestr
;             = {galperbin:nperbin,massbins:massbins_now,ave_mass:ave_mass_now,$
;                ave_mass_dev:ave_mass_dev_now,ave_feh:ave_feh_now,$
;                ave_feh_dev:ave_feh_dev_now,hifeh:hifeh_now,lofeh:lofeh_now,$
;                ave_mgh:ave_mgh_now,ave_mgh_dev:ave_mgh_dev_now,himgh:himgh_now,lomgh:lomgh_now}
   restore,linfitparam_fixslope_file
;         tconst_sdcl_feh,tconst_sdcl_mgh,tconst_sdms_feh,$
;         tconst_sdms_mgh,tconst_clms_feh,tconst_clms_mgh,$
;         tslope_sdcl_feh,tslope_sdcl_mgh,tslope_sdms_feh,$
;         tslope_sdms_mgh,tslope_clms_feh,tslope_clms_mgh,$
;         fehcommonslope,fehcommonslope_err,fehcommonslope_stdev,$
;         mghcommonslope,mghcommonslope_err,mghcommonslope_stdev,$
;         sdssmzrpar_fixslope,cl0024mzrpar_fixslope,ms0451mzrpar_fixslope,$
;         sdssmzrerr_fixslope,cl0024mzrerr_fixslope,ms0451mzrerr_fixslope,$
;         sdssmzrchi_fixslope,cl0024mzrchi_fixslope,ms0451mzrchi_fixslope,$
;         sdssmarpar_fixslope,cl0024marpar_fixslope,ms0451marpar_fixslope,$
;         sdssmarerr_fixslope,cl0024marerr_fixslope,ms0451marerr_fixslope,$
;         sdssmarchi_fixslope,cl0024marchi_fixslope,ms0451marchi_fixslope,$

   restore, linfitparam_file
;        allmzrpar,allmarpar,sdssmzrpar,sdssmarpar,cl0024mzrpar,cl0024marpar,ms0451mzrpar,ms0451marpar,$
;        allmzrerr,allmarerr,sdssmzrerr,sdssmarerr,cl0024mzrerr,cl0024marerr,ms0451mzrerr,ms0451marerr,$
;        chiallmzr,chiallmar,chisdssmzr,chisdssmar,chicl0024mzr,chicl0024mar,chims0451mzr,chims0451mar,$

   anova_feh = mrdfits(fileanova,1)
   anova_mgh = mrdfits(fileanova,2)
;        outstr_feh = replicate({mbinmin:0d,mbinmax:0d,f_anova:0d,dof_ssb:0,dof_ssw:0,$
;                           SSW:0d,SSB:0d,TSS:0d,metal_mean:[0.,0.,0.],metal_sigmas:[0.,0.,0.]},nbins)

   goodfit = where(sciall.goodsn eq 1,cgoodfit,complement=badfit,ncomplement=cbadfit)
   symsizepar = [0.7,1.5/(max(sciall(goodfit).snfit)-min(sciall(goodfit).snfit))]
   symsize = symsizepar[1]*(sciall.snfit-min(sciall(goodfit).snfit))+symsizepar[0]
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   psname= 'MZR_'+prefix+'.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange = [9.,12]
      yrange = [-1.,0.3]
      plot,sciall.logmstar,sciall.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5

;      polyfill,[sdss_avestr.massbins,reverse(sdss_avestr.massbins)],$
;          [sdss_avestr.hifeh,sdss_avestr.lofeh],color=fsc_color(shadecolor[0])
;      polyfill,[cl0024_avestr.massbins,reverse(cl0024_avestr.massbins)],$
;          [cl0024_avestr.hifeh,cl0024_avestr.lofeh],color=fsc_color(shadecolor[1])
;      polyfill,[ms0451_avestr.massbins,reverse(ms0451_avestr.massbins)],$
;          [ms0451_avestr.hifeh,ms0451_avestr.lofeh],color=fsc_color(shadecolor[2])


      for i=0,n_elements(sciall)-1 do begin
        if sciall[i].zspec lt 0.3 then symcol = clustercolor[0]
        if sciall[i].zspec gt 0.3 and sciall[i].zspec lt 0.5 then symcol = clustercolor[1] 
        if sciall[i].zspec gt 0.5 and sciall[i].zspec lt 0.6 then symcol = clustercolor[2]

        if sciall[i].goodsn eq 1 then symnum = 46 else symnum = 16
        if sciall[i].goodsn eq 1 then begin   
            cgerrplot,[sciall[i].logmstar],[sciall[i].fehlower],[sciall[i].fehupper],color=symcol,$
               thick=0.5
        endif
        cgPlot,[sciall[i].logmstar],[sciall[i].feh], Color=symcol, PSym=symnum, SymColor=symcol, $
           SymSize=symsize[i],/overplot
;        oplot,[sciall[i].logmstar],[sciall[i].feh],psym=cgsymcat(46),color=fsc_Color(symcol),$
;              symsize=symsize[i]
      endfor
      ;overplot with linear fit
      oplot,!x.crange,(!x.crange-10.)*sdssmzrpar[1]+sdssmzrpar[0],linestyle=2,$
            thick=2,color=fsc_color(meancolor[0])
      oplot,!x.crange,(!x.crange-10.)*cl0024mzrpar[1]+cl0024mzrpar[0],linestyle=2,$
            thick=2,color=fsc_color(meancolor[1])
      oplot,!x.crange,(!x.crange-10.)*ms0451mzrpar[1]+ms0451mzrpar[0],linestyle=2,$
            thick=2,color=fsc_color(meancolor[2])

      oplot,!x.crange,(!x.crange-10.)*sdssmzrpar_fixslope[1]+sdssmzrpar_fixslope[0],$
            thick=2,color=fsc_color(meancolor[0])
      oplot,!x.crange,(!x.crange-10.)*cl0024mzrpar_fixslope[1]+cl0024mzrpar_fixslope[0],$
            thick=2,color=fsc_color(meancolor[1])
      oplot,!x.crange,(!x.crange-10.)*ms0451mzrpar_fixslope[1]+ms0451mzrpar_fixslope[0],$
            thick=2,color=fsc_color(meancolor[2])

      ;overplot with means used in anova
      for i=0,2 do begin
         cgerrplot,anova_feh.metal_mean[i],anova_feh.mbinmin,anova_feh.mbinmax,/horizontal,color=meancolor[i]
         cgerrplot,(anova_feh.mbinmin+anova_feh.mbinmax)*0.5,anova_feh.metal_mean[i]-anova_feh.metal_sigmas[i],$
                   anova_feh.metal_mean[i]+anova_feh.metal_sigmas[i],color=meancolor[i]
         cgplot,(anova_feh.mbinmin+anova_feh.mbinmax)*0.5,anova_feh.metal_mean[i],psym=14,$
                symcolor=meancolor[i],/overplot,symsize=3
      endfor

      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

      al_legend,['z~0 SDSS','z~0.4 Cl0024','z~0.54 MS0451'],psym=cgsymcat(46),$
                 box=0,colors=['skyblue','goldenrod','crimson'],$
                 /bottom,/right,charsize=1,symsize=2
   device,/close

   psname= 'MMgR_'+prefix+'.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange = [9.,12]
      yrange = [-1.,0.5]
      plot,sciall.logmstar,sciall.ah,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5

;      polyfill,[sdss_avestr.massbins,reverse(sdss_avestr.massbins)],$
;          [sdss_avestr.himgh,sdss_avestr.lomgh],color=fsc_color(shadecolor[0])
;      polyfill,[cl0024_avestr.massbins,reverse(cl0024_avestr.massbins)],$
;          [cl0024_avestr.himgh,cl0024_avestr.lomgh],color=fsc_color(shadecolor[1])
;      polyfill,[ms0451_avestr.massbins,reverse(ms0451_avestr.massbins)],$
;          [ms0451_avestr.himgh,ms0451_avestr.lomgh],color=fsc_color(shadecolor[2])

      for i=0,n_elements(sciall)-1 do begin
        if sciall[i].zspec lt 0.3 then symcol = clustercolor[0]
        if sciall[i].zspec gt 0.3 and sciall[i].zspec lt 0.5 then symcol = clustercolor[1] 
        if sciall[i].zspec gt 0.5 and sciall[i].zspec lt 0.6 then symcol = clustercolor[2]

        if sciall[i].goodsn eq 1 then symnum = 46 else symnum = 16
        if sciall[i].goodsn eq 1 then begin   
            cgerrplot,[sciall[i].logmstar],[sciall[i].ahlower],[sciall[i].ahupper],$
            color=symcol,thick=0.5
        endif
        cgPlot,[sciall[i].logmstar],[sciall[i].ah], Color=symcol, PSym=symnum, SymColor=symcol, $
           SymSize=symsize[i],/overplot
;        oplot,[sciall[i].logmstar],[sciall[i].ah],psym=cgsymcat(46),color=fsc_Color(symcol),$
;              symsize=symsize[i]
      endfor
      ;overplot with linear fit
      oplot,!x.crange,(!x.crange-10.)*sdssmarpar[1]+sdssmarpar[0],linestyle=2,$
            thick=2,color=fsc_color(meancolor[0])
      oplot,!x.crange,(!x.crange-10.)*cl0024marpar[1]+cl0024marpar[0],linestyle=2,$
            thick=2,color=fsc_color(meancolor[1])
      oplot,!x.crange,(!x.crange-10.)*ms0451marpar[1]+ms0451marpar[0],linestyle=2,$
            thick=2,color=fsc_color(meancolor[2])

      oplot,!x.crange,(!x.crange-10.)*sdssmarpar_fixslope[1]+sdssmarpar_fixslope[0],$
            thick=2,color=fsc_color(meancolor[0])
      oplot,!x.crange,(!x.crange-10.)*cl0024marpar_fixslope[1]+cl0024marpar_fixslope[0],$
            thick=2,color=fsc_color(meancolor[1])
      oplot,!x.crange,(!x.crange-10.)*ms0451marpar_fixslope[1]+ms0451marpar_fixslope[0],$
            thick=2,color=fsc_color(meancolor[2])

      ;overplot with means used in anova
      for i=0,2 do begin
         cgerrplot,anova_mgh.metal_mean[i],anova_mgh.mbinmin,anova_mgh.mbinmax,/horizontal,color=meancolor[i]
         cgerrplot,(anova_mgh.mbinmin+anova_mgh.mbinmax)*0.5,anova_mgh.metal_mean[i]-anova_mgh.metal_sigmas[i],$
                   anova_mgh.metal_mean[i]+anova_mgh.metal_sigmas[i],color=meancolor[i]
         cgplot,(anova_mgh.mbinmin+anova_mgh.mbinmax)*0.5,anova_mgh.metal_mean[i],psym=14,$
                symcolor=meancolor[i],/overplot,symsize=3
      endfor

      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

      al_legend,['z~0 SDSS','z~0.4 Cl0024','z~0.55 MS0451'],psym=cgsymcat(46),$
                 box=0,colors=['skyblue','goldenrod','crimson'],$
                 /bottom,/right,charsize=1,symsize=2
   device,/close
end

pro plot_byage,sciall,prefix=prefix
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                  fileancova,deviationparam_file
common color_pallete, clustercolor, shadecolor, meancolor

   plotfixslope = 0
   plotancova = 1

   restore, linfitparam_file
   restore,linfitparam_fixslope_file
   restore, fileancova
   choi = get_refvalue(5,/choi14)
   choi14z01 = choi[where(choi.zlow eq 0.1)]
   choi14z04 = choi[where(choi.zlow eq 0.4)]
   choi14z02 = choi[where(choi.zlow eq 0.2)]
   choi14z06 = choi[where(choi.zlow eq 0.55)]
   
   if plotfixslope then begin
     sdssmzrpar = sdssmzrpar_fixslope
     cl0024mzrpar = cl0024mzrpar_fixslope
     cl0024mzrpar[0] = cl0024mzrpar[0]-0.02
     ms0451mzrpar = ms0451mzrpar_fixslope
     sdssmarpar = sdssmarpar_fixslope
     cl0024marpar = cl0024marpar_fixslope
     ms0451marpar = ms0451marpar_fixslope
   endif
   if plotancova then begin
     sdssmzrpar = [const_feh_mo2,coeff_feh_mo2[0]]
     cl0024mzrpar = [const_feh_mo2+coeff_feh_mo2[1],coeff_feh_mo2[0]+coeff_feh_mo2[3]]
     ms0451mzrpar = [const_feh_mo2+coeff_feh_mo2[2],coeff_feh_mo2[0]+coeff_feh_mo2[4]]
     sdssmarpar = [const_mgh_mo3,coeff_mgh_mo3[0]]
     cl0024marpar = [const_mgh_mo3,coeff_mgh_mo3[0]]
     ms0451marpar = [const_mgh_mo3,coeff_mgh_mo3[0]]
   endif

   naiman = get_refvalue(55,/naiman18)
   saracco = get_refvalue(55,/saracco19)
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
   agerange = [0,13,0,2.,4.,6.,8.,13.]  
   for n=0,1 do begin 
;   for n=0,n_Elements(agerange)-2 do begin
     inrange = where(sciall.ageform ge agerange(n) and sciall.ageform lt agerange(n+1),cinrange)
     if cinrange eq 0 then continue
     sci = sciall(inrange)
     goodness =sciall(inrange).goodsn
     goodfit = where(goodness eq 1,cgoodfit, complement=badfit,ncomplement=cbadfit)
     name = 'univage_'+strtrim(string(agerange(n),format='(I)'),2)+'_'+strtrim(string(agerange(n+1),format='(I)'),2) 
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      psname= name+'_feh_mass.eps'
      device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[9.,12]
         yrange=[-1.,0.3]
         plot,sci.logmstar,sci.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M!L*!N/M'+sunsym+')'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   
;         cgerrplot,sci(goodfit).logmstar,sci(goodfit).fehlower,sci(goodfit).fehupper,color='crimson',thick=0.5
;         if cbadfit gt 0 then oplot,sci(badfit).logmstar,sci(badfit).feh,psym=cgsymcat(16),color=fsc_Color('darkgray'),symsize=0.7
         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'crimson'

           oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).feh],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         oplot,!x.crange,(!x.crange-10.)*sdssmzrpar[1]+sdssmzrpar[0],linestyle=0,thick=2,color=fsc_color('navyblue')
         oplot,!x.crange,(!x.crange-10.)*cl0024mzrpar[1]+cl0024mzrpar[0],linestyle=0,thick=2,color=fsc_color('peru')
         oplot,!x.crange,(!x.crange-10.)*ms0451mzrpar[1]+ms0451mzrpar[0],linestyle=0,thick=2,color=fsc_color('darkred')
         oplot,naiman.mass_feh,naiman.feh,color=fsc_color('darkorchid'),thick=3, linestyle=2 
         oplot,choi14z01.mass,choi14z01.feh,psym=cgsymcat(16),color=fsc_color('ygb5'),symsize=1
         oplot,choi14z04.mass,choi14z04.feh,psym=cgsymcat(16),color=fsc_color('darkgoldenrod'),symsize=1
         ;add Kriek2016
         oploterror,[11.5],[-0.25],[0.1],[0.11],color=fsc_color('firebrick'),psym=16,symsize=1.5
         ;oplot,choi14z04.mass,choi14z04.feh,color=fsc_color('darkgoldenrod'),linestyle=0
;         cgerrplot,saracco.mass,saracco.zh_lower,saracco.zh_upper,color='darkred',psym=16
;         oplot,saracco.mass,saracco.zh,psym=cgsymcat(16),color=fsc_color('darkred'),symsize=1
            
         ;xyouts,11.3,-0.8,strtrim(string(agerange(n),format='(I)'),2)+'-'+strtrim(string(agerange(n+1),format='(I)'),2)+'Gyr'

         ;typical error bars
         locsdss = where(goodness eq 1 and sci.zspec lt 0.3)
         uppersdss= median(sci(locsdss).fehupper-sci(locsdss).feh)
         lowersdss = median(sci(locsdss).feh-sci(locsdss).fehlower)
         loccl0024 = where(goodness eq 1 and sci.zspec gt 0.3 and sci.zspec lt 0.5)
         uppercl0024 = median(sci(loccl0024).fehupper-sci(loccl0024).feh)
         lowercl0024 = median(sci(loccl0024).feh-sci(loccl0024).fehlower)
         locms0451 = where(goodness eq 1 and sci.zspec gt 0.5)
         upperms0451 = median(sci(locms0451).fehupper-sci(locms0451).feh)
         lowerms0451 = median(sci(locms0451).feh-sci(locms0451).fehlower)
         ;labelling
         ylab = -0.85
         xlab = 9.6
         cgplot,[xlab],[ylab],err_xhigh=[0.1],err_xlow=[0.1],err_yhigh=uppersdss,err_ylow=lowersdss,$
                color='skyblue',psym=46,symsize=2,/overplot
         xyouts,xlab+0.15,ylab-0.025,'z~0.05',charsize=1
         xlab = 10.4
         cgplot,[xlab],[ylab],err_xhigh=[0.1],err_xlow=[0.1],err_yhigh=uppercl0024,err_ylow=lowercl0024,$
                color='goldenrod',psym=46,symsize=2,/overplot
         xyouts,xlab+0.15,ylab-0.025,'z~0.39',charsize=1
         xlab = 11.2
         cgplot,[xlab],[ylab],err_xhigh=[0.1],err_xlow=[0.1],err_yhigh=upperms0451,err_ylow=lowerms0451,$
                color='crimson',psym=46,symsize=2,/overplot
         xyouts,xlab+0.15,ylab-0.025,'z~0.54',charsize=1
         ;label
;      Fraction under estimatedymcat(46),box=0,colors=['skyblue','goldenrod','crimson'],$,psym=1
;                   charsize=1,sy!xmsize=2,position=[10.6,-0.6]

      device,/close
      vincenzo = get_refvalue(55,/vincenzo18) 
      psname=name+'_feh_alpha.eps'
      device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[-0.6,0.3]
         yrange=[-0.5,0.9]
         plot,sci.feh,sci.alphafe,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='[Fe/H]'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/Fe]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   
   ;      cgerrplot,sci(goodfit).alphafe,sci(goodfit).alphafelower,sci(goodfit).alphafeupper,$
   ;                color='crimson',thick=0.5
   ;      cgerrplot,sci(goodfit).feh,sci(goodfit).fehlower,sci(goodfit).fehupper,$
   ;                color='crimson',thick=0.5,/horizontal
;         if cbadfit gt 0 then oplot,sci(badfit).feh,sci(badfit).alphafe,psym=cgsymcat(16),color=fsc_Color('darkgray'),symsize=0.7
         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'crimson'

           oplot,[sci(goodfit[i]).feh],[sci(goodfit[i]).alphafe],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         oplot,vincenzo.feh,vincenzo.mgfe,linestyle=0,color=fsc_color('darkorchid'),thick=4
         ;add Kriek2016
         oploterror,[-0.25],[0.59],[0.11],[0.11],color=fsc_color('firebrick'),psym=16,symsize=1.5
         oplot,choi14z01.feh,choi14z01.mgfe,psym=cgsymcat(16),color=fsc_color('ygb5'),symsize=1.5
         oplot,choi14z04.feh,choi14z04.mgfe,psym=cgsymcat(16),color=fsc_color('darkgoldenrod'),symsize=1.5
;         xyouts,-0.9,-0.4,strtrim(string(agerange(n),format='(I)'),2)+'-'+strtrim(string(agerange(n+1),format='(I)'),2)+'Gyr'
      device,/close
   
      psname = name+'_alpha_mass.eps'
      device,filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[9.,12]
         yrange=[-1.,0.9]
         plot,sci.logmstar,sci.ah,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M!L*!N/M'+sunsym+')'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/H]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
  ;       for i=0,cgoodfit-1 do begin
  ;         oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).feh],psym=cgsymcat(46),color=fsc_Color('lightgray'),symsize=goodsymsize[i]
  ;       endfor
   
  ;       if cbadfit gt 0 then oplot,sci(badfit).logmstar,sci(badfit).ah,psym=cgsymcat(46),color=fsc_Color('darkgray'),symsize=0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'crimson'
           oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).ah],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         ;Add Choi's data
   ;      oploterror,choi14z01.mass,choi14z01.feh+choi14z01.mgfe,sqrt(choi14z01.feherr^2+choi14z01.mgfe^2),color=fsc_color('ygb5'),linethick=2,errcolor=fsc_color('ygb5')
         oplot,choi14z01.mass,choi14z01.feh+choi14z01.mgfe,psym=cgsymcat(16),color=fsc_color('ygb5'),symsize=1.
 
   ;      oploterror,choi14z04.mass,choi14z04.feh+choi14z04.mgfe,sqrt(choi14z04.feherr^2+choi14z04.mgfe^2),color=fsc_color('darkgoldenrod'),linethick=2,errcolor=fsc_color('darkgoldenrod')
         oplot,choi14z04.mass,choi14z04.feh+choi14z04.mgfe,psym=cgsymcat(16),color=fsc_color('darkgoldenrod'),symsize=1.
         
         ;add Kriek2016
         oploterror,[11.5],[0.59-0.25],[0.1],[0.16],color=fsc_color('firebrick'),psym=16,symsize=1.5
         oplot,naiman.mass_mgh,naiman.mgh,color=fsc_color('darkorchid'),thick=3,linestyle=2      
         ;add best fit lines
         oplot,!x.crange,(!x.crange-10.)*sdssmarpar[1]+sdssmarpar[0],linestyle=0,thick=2,color=fsc_color('navyblue')
         oplot,!x.crange,(!x.crange-10.)*cl0024marpar[1]+cl0024marpar[0],linestyle=0,thick=2,color=fsc_color('peru')
         oplot,!x.crange,(!x.crange-10.)*ms0451marpar[1]+ms0451marpar[0],linestyle=0,thick=2,color=fsc_color('indianred')
         oplot,!x.crange,(!x.crange-10.)*sdssmarpar[1]+sdssmarpar[0],linestyle=0,thick=2,color=fsc_color('black') 
;         al_legend,['z~0','z~0.39','z~0.55'],psym=cgsymcat(46),box=0,symsize=2,$
;                colors=['skyblue','goldenrod','crimson'],charsize=1.1,position=[11.1,-0.4]
;         al_legend,['Choi+2014'],psym=cgsymcat(16),box=0,symsize=1.5,colors=['darkgray'],charsize=1.1,position=[11.1,-0.8]
;         xyouts,11.5,-0.4,'Observed Redshift',alignment=0.5,charsize=1.1
         ;typical error bars
         uppersdss = median(sci(locsdss).ahupper-sci(locsdss).ah)
         lowersdss = median(sci(locsdss).ah-sci(locsdss).ahlower)
         uppercl0024 = median(sci(loccl0024).ahupper-sci(loccl0024).ah)
         lowercl0024 = median(sci(loccl0024).ah-sci(loccl0024).ahlower)
         upperms0451 = median(sci(locms0451).ahupper-sci(locms0451).ah)
         lowerms0451 = median(sci(locms0451).ah-sci(locms0451).ahlower)
         ;labelling
         ylab = -0.7
         xlab = 9.6
         cgplot,[xlab],[ylab],err_xhigh=[0.1],err_xlow=[0.1],err_yhigh=uppersdss,err_ylow=lowersdss,$
                color='skyblue',psym=46,symsize=2,/overplot
         xyouts,xlab+0.15,ylab-0.025,'z~0.05',charsize=1
         xlab = 10.4
         cgplot,[xlab],[ylab],err_xhigh=[0.1],err_xlow=[0.1],err_yhigh=uppercl0024,err_ylow=lowercl0024,$
                color='goldenrod',psym=46,symsize=2,/overplot
         xyouts,xlab+0.15,ylab-0.025,'z~0.39',charsize=1
         xlab = 11.2
         cgplot,[xlab],[ylab],err_xhigh=[0.1],err_xlow=[0.1],err_yhigh=upperms0451,err_ylow=lowerms0451,$
                color='crimson',psym=46,symsize=2,/overplot
         xyouts,xlab+0.15,ylab-0.025,'z~0.54',charsize=1


      device,/close

      psname=name+'_mass_alphafe.eps'
      device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[9.5,12]
         yrange=[-0.5,0.9]
         plot,sci.logmstar,sci.alphafe,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M!L*!N/M'+sunsym+')'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/Fe]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

   ;      cgerrplot,sci(goodfit).alphafe,sci(goodfit).alphafelower,sci(goodfit).alphafeupper,$
   ;                color='crimson',thick=0.5
   ;      cgerrplot,sci(goodfit).feh,sci(goodfit).fehlower,sci(goodfit).fehupper,$
   ;                color='crimson',thick=0.5,/horizontal
         if cbadfit gt 0 then oplot,sci(badfit).logmstar,sci(badfit).alphafe,psym=cgsymcat(16),color=fsc_Color('darkgray'),symsize=0.7
         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'crimson'

           oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).alphafe],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         xyouts,-0.9,-0.4,strtrim(string(agerange(n),format='(I)'),2)+'-'+strtrim(string(agerange(n+1),format='(I)'),2)+'Gyr'
      device,/close

      psname=name+'_age_alphafe.eps'
      device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[0,13]
         yrange=[-0.5,0.9]
         plot,sci.age,sci.alphafe,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Age(Gyr)'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/Fe]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'crimson'

           oplot,[sci(goodfit[i]).age],[sci(goodfit[i]).alphafe],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         ;overplot with walcher 2016
         oplot,[0,9],[0,9]*0.009-0.01,thick=2
         oplot,[9,13],[9,13]*0.031-0.199,thick=2
         ;overplot with choi14
         oplot,choi14z01.age,choi14z01.mgfe,psym=cgsymcat(16),color=fsc_color('ygb5'),symsize=1.2
;         oplot,choi14z02.age,choi14z02.mgfe,psym=cgsymcat(16),color=fsc_color('darkgreen'),symsize=1.
         oplot,choi14z04.age,choi14z04.mgfe,psym=cgsymcat(16),color=fsc_color('darkgoldenrod'),symsize=1.2
;         oplot,choi14z06.age,choi14z06.mgfe,psym=cgsymcat(16),color=fsc_color('indianred'),symsize=1.
         ;add Kriek2016
         oploterror,[2.71],[0.59],[0.22],[0.11],color=fsc_color('firebrick'),psym=16,symsize=1.5
      device,/close
   endfor
end

pro plot_feh_alpha_byage,sciall
   agerange = [0.5,2.,4.,6.,8.,13.]
   colors = ['slateblue','thistle','gold','coral','firebrick']
   vincenzo = get_refvalue(55,/vincenzo18) 

   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()

   psname='feh_alpha_byage.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   xrange=[-0.6,0.3]
   yrange=[-0.5,0.9]
   firstplot = 0
   
   goodall = where(sciall.ageform gt 0.5 and sciall.goodsn eq 1)
   symsize = 1.5/(max(sciall(goodall).snfit)-min(sciall(goodall).snfit))*sciall.snfit+0.7   

   for n=n_Elements(agerange)-2,0,-1 do begin
      inrange = where(sciall.ageform ge agerange(n) and sciall.ageform lt agerange(n+1)$
                      and sciall.goodsn eq 1,cinrange)
      if cinrange eq 0 then continue
      sci = sciall(inrange)
      symsizenow = symsize(inrange)
      if firstplot eq 0 then begin
         plot,sci.feh,sci.alphafe,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='[Fe/H]'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/Fe]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
         firstplot = 99
      endif

      for i=0,cinrange-1 do begin
        if sci[i].zspec lt 0.3 then symnum = 6
        if sci[i].zspec gt 0.3 and sci[i].zspec lt 0.5 then symnum = 5
        if sci[i].zspec gt 0.5 and sci[i].zspec lt 0.6 then symnum = 4

        oplot,[sci[i].feh],[sci[i].alphafe],psym=cgsymcat(symnum),color=fsc_Color(colors[n]),symsize=symsizenow[i]
      endfor
    endfor
    oplot,vincenzo.feh,vincenzo.mgfe,linestyle=0,color=fsc_color('darkorchid'),thick=4
    device,/close
stop
end

pro calculate_nta_data, sciall,proball,coeff,coeff_err
   ;set up
   R = 0.46 ;return fraction
   ;solar abundances from asplund2009 (number of atoms)
   solarmgh = 10^(7.64-12.)*24. ;times 24 to convert number of atoms to mass. cause Mg is 24 mass number.
   yieldmgh = 3.*solarmgh
;   yieldmgh = 6.*solarmgh
   ;prep data
   won = where(sciall.goodsn eq 1,cwon)
   sci = sciall[won]
   prob = proball[won]
   mghtest = solarmgh*10.^(sci.ah)
   ntamghtest = (1.-R)*(yieldmgh/mghtest-1.)
   ;take out outliers
   outliers = where((ntamghtest lt 0.1 and sci.logmstar lt 10.2) or (ntamghtest lt 0.01),$
                 coutliers,complement=good,ncomplement=nobjs)
   print, 'there are', coutliers,'outliders'

   ;get final data
   sci = sci(good)
   prob = prob(good)
   mgh = sci.ah
   dmgh = prob.daharr
   probdmgh = prob.probdah
   mass = sci.logmstar
   ;calcualte for plotting
   ;real data
   mghreal = solarmgh*10.^(mgh)
   ntamgh = (1.-R)*(yieldmgh/mghreal-1.)
   goodsymsize = 1.5/(max(sci.snfit)-min(sci.snfit))*sci.snfit+0.7

   ;functional form monte carlo
   ;getting cumulative probability
   dummy = size(dmgh,/dimensions)
   ngrid = dummy[0]
   cumprob = fltarr(ngrid,nobjs)
   for i=0,nobjs-1 do begin
      curyarr = dmgh[*,i]
      curprob = probdmgh[*,i]
      for j=1,ngrid-1 do cumprob[j,i] = int_tabulated(curyarr[0:j],curprob[0:j])
   endfor

   nmc = 1000
   linparmc = fltarr(2,nmc)
   ntaarray = fltarr(nobjs,nmc)
   for i = 0,nmc-1 do begin
      logmghrdm = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         logmghrdm(j) = mgh(j)+interpol(dmgh[*,j],cumprob[*,j],randomarr(j))
      endfor

      mghrdm = (solarmgh*10.^(logmghrdm))<(yieldmgh*0.999)
      ntardm = (1.-R)*(yieldmgh/mghrdm-1.)
      ntaarray[*,i] = ntardm
      ntardm = alog10(ntardm)
      linpar_rdm = linfit(mass-10.,ntardm)
      linparmc[*,i] = linpar_rdm
   endfor
   ntalims = fltarr(3,nobjs)
   for i=0,nobjs-1 do begin
      ntalims[*,i] = percentiles(ntaarray[i,*],value=[0.16,0.5,0.85])
   endfor
   ntaerr = stddev(ntaarray,dimension=2)
   linparmgh = mean(linparmc,dimension=2)
   linparmgh_err = stddev(linparmc,dimension=2)
   print, 'CALCULATING NTA FROM DATA DIRECTLY'
   print, 'dependence of ntamgh on M [mean,stdev]:', linparmgh[1],linparmgh_err[1]

   save,sci,ntamgh,ntalims,ntaerr,filename='nta_data.sav'
   ;create high and low shading
   nmass = 50.
   massline = 9.05+findgen(nmass)/nmass*2.45
   fitlogntamgh = (massline-10.)*linparmgh[1]+linparmgh[0]

   scattermgh=0.2
   logmghline = coeff[0]+coeff[1]*(massline-10.)
   mghline = solarmgh*10.^(logmghline)
   mghlinehi = solarmgh*10.^(logmghline+scattermgh)
   mghlinelo = solarmgh*10.^(logmghline-scattermgh) ;assume intrinsic scatter=0.1
   ntamghline = (1.-R)*(yieldmgh/mghline-1.)
   ntamghlinehi = (1.-R)*(yieldmgh/mghlinelo-1.)
   ntamghlinelo = (1.-R)*(yieldmgh/mghlinehi-1.)

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;literature values
   ;spitoni2010
   spitoni_low = 0.11
   spitoni_high = 5.4

   ;upper limits from Lu2015
   readcol,'lu2015_upperlim.dat',logmass_lu,logeta_lu,comment='#'
   ;Chrisholm 2017
   mass_jc = [9.3,9.5,10.1,10.3,10.7]
   nta_jc = [0.91,2,1.2,0.23,0.54]
   ntaerr_jc = [0.33,0.57,0.51,0.08,0.18]
   ;Muratov2015
   mass_muratov = [9.,9.5,10.,10.5,11.,11.5]
   nta_muratov = 3.6*(10.^(-0.35*(mass_muratov-10.)))
   ;bulatto2013
   mass_bulatto = [9.4]
   nta_bulatto = [3.2]
   ntaupper_bulatto = [10.7]
   ntalower_bulatto = [1.1]

   nelson = get_refvalue(55,/nelson19) 
   ;;;;;;;;;;;;;;;;;;
;   set_plot,'ps'
;   !p.font = 0
;   !p.charsize=1.2
;   ntaletter = "150B
;   proptoletter = "265B
;   proptoletter = "265B
;   pmsymbol = "261B
;   sunsym = sunsymbol()
;   psname='nta_mass_data.eps'
;   device, filename = psname,xsize = 15,ysize = 10, $
;        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;      xrange=[9,11.5]
;      yrange=[0.01,10.]
;      plot,mass,ntamgh,/ylog,/nodata,xstyle=5,ystyle=5,xrange=xrange,yrange=yrange
;      polyfill,[xrange,reverse(xrange)],[spitoni_high,spitoni_high,spitoni_low,spitoni_low],color=fsc_color('lightgray')
;      x=[massline,reverse(massline)]
;      y=[ntamghlinehi,reverse(ntamghlinelo)]
;      polyfill,x,y,color=fsc_color('pink')
;
;      for i=0,nobjs-1 do begin
;           if sci[i].zspec lt 0.3 then symcol = 'skyblue'
;           if sci[i].zspec gt 0.3 and sci[i].zspec lt 0.5 then symcol = 'goldenrod'
;           if sci[i].zspec gt 0.5 and sci[i].zspec lt 0.6 then symcol = 'crimson'
;           oplot,[mass[i]],[ntamgh[i]],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
;      endfor
;
;      oplot,massline,10.^fitlogntamgh,color=fsc_color('darkred')
;      oplot,logmass_lu, 10.^logeta_lu,linestyle=2,color=fsc_color('navy')
;      oplot,mass_muratov,nta_muratov,psym=0,linestyle=1,thick=2
;;      zahid=get_refvalue(55,/zahid12)
;;      oplot,zahid.mass_low,zahid.eta_low,linestyle=2,color=fsc_color('darkorchid')
;;      oplot,zahid.mass_high,zahid.eta_high,linestyle=2,color=fsc_color('darkorchid')
;      oploterror,mass_jc,nta_jc,ntaerr_jc,psym=1,color=fsc_color('navy')
;      cgplot,mass_jc,nta_jc,psym=16,color=fsc_color('navy'),/overplot
;      
;;      oplot,nelson.massv0,nelson.etav0,color=fsc_color('blue')
;      oplot,nelson.massv150,nelson.etav150,color=fsc_color('orange'),linestyle=3
;      oplot,nelson.massv250,nelson.etav250,color=fsc_color('darkred'),linestyle=3
;      xyouts,9.1,0.02,'!9'+string(ntaletter)+' !9'+string(proptoletter)+' !xM!D*!N!E'+sigfig(linparmgh(1),2)+'!9'+string(pmsymbol)+'!x'+sigfig(linparmgh_err(1),1),charsize=1.5
;      axis,xaxis=0,xrange=xrange,xtitle='Log(M!L*!N/M'+sunsym+')'
;      axis,xaxis=1,xrange=xrange,xtickformat='(A1)'
;      axis,yaxis=0,yrange=yrange,ytitle='!9'+string(ntaletter)+'!x'
;      axis,yaxis=1,yrange=yrange,ytickformat='(A1)'
;   device,/close
;calculate median uncertainties
separate_cluster,sci,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451  
print, 'uncertainties in nta, mean ,median'
print,'all',alog10(mean(ntaerr/ntamgh)),alog10(median(ntaerr/ntamgh))
print,'sdss',alog10(mean(ntaerr(wonsdss)/ntamgh(wonsdss))),alog10(median(ntaerr(wonsdss)/ntamgh(wonsdss)))
print,'cl0024',alog10(mean(ntaerr(woncl0024)/ntamgh(woncl0024))),alog10(median(ntaerr(woncl0024)/ntamgh(woncl0024)))
print,'ms0451',alog10(mean(ntaerr(wonms0451)/ntamgh(wonms0451))),alog10(median(ntaerr(wonms0451)/ntamgh(wonms0451)))
stop
end

pro calculate_nta_function, coeff, coefferr
;input: coeff = [intercept, slope] of the Mg-mass relation

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
   mass = 9.05+findgen(nmass)/nmass*2.45

   logmgh = coeff[0]+coeff[1]*(mass-10.) ;mgh
   scattermgh = 0.2 ;mgh
   mgh = solarmgh*10.^(logmgh)
   mghhi = solarmgh*10.^(logmgh+scattermgh)
   mghlo = solarmgh*10.^(logmgh-scattermgh) ;assume intrinsic scatter=0.1
   ntamgh = (1.-R)*(yieldmgh/mgh-1.)
   ntamghhi = (1.-R)*(yieldmgh/mghlo-1.)
   ntamghlo = (1.-R)*(yieldmgh/mghhi-1.)

   ;functional form monte carlo
   nmc = 1000
   linparmc = fltarr(2,nmc)
   scatter = 0.2
   for i = 0,nmc-1 do begin
      coeffrdm = randomn(seed,2)*coefferr+coeff
      massrdm = randomu(seed,180)*2.5+9.
      logmghrdm = coeffrdm[0]+coeffrdm[1]*(massrdm-10.)+randomn(seed,nmc)*scatter
      mghrdm = (solarmgh*10.^(logmghrdm))<(yieldmgh*0.999)
      ntardm = alog10((1.-R)*(yieldmgh/mghrdm-1.))
      linpar_rdm = linfit(massrdm-10.,ntardm)
      linparmc[*,i] = linpar_rdm
   endfor
   linparmgh = mean(linparmc,dimension=2)
   linparmgh_err = stddev(linparmc,dimension=2)
   fitlogntamgh = mass*linparmgh[1]+linparmgh[0]
   print, 'CALCULATING NTA FROM LINEAR MZR'
   print, 'dependence of ntamgh on M [mean,stdev]:', linparmgh[1],linparmgh_err[1]
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   spitoni_low = 0.11
   spitoni_high = 5.4
   set_plot,'ps'
   !p.font = 0
   ntaletter = "150B
   proptoletter = "265B
   sunsym = sunsymbol()
   psname='nta_mass_function.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[9,11.5]
      yrange=[0.01,10.]
      ;plot,mass,ntamgh,xtitle='Log(M/M'+sunsym+')',ytitle='!9'+string(ntaletter)+'!x',/ylog,/nodata,xrange=xrange,yrange=[0.01,10.],xstyle=5,ystyle=5
      plot,mass,ntamgh,xtitle='Log(M/M'+sunsym+')',/ylog,/nodata,xstyle=5,ystyle=5,xrange=xrange,yrange=yrange
      polyfill,[xrange,reverse(xrange)],[spitoni_high,spitoni_high,spitoni_low,spitoni_low],color=fsc_color('lightgray')
      x=[mass,reverse(mass)]
      y=[ntamghhi,reverse(ntamghlo)]
      polyfill,x,y,color=fsc_color('salmon')
      oplot,mass,ntamgh
      oplot,mass,10.^fitlogntamgh,color=fsc_color('darkred'),linestyle=5
      xyouts,10.9,4,'!9'+string(ntaletter)+' !9'+string(proptoletter)+' !xM!D*!N!E'+sigfig(linparmgh(1),2),charsize=1.5
      axis,xaxis=0,xrange=xrange,xtitle='Log(M/M'+sunsym+')'
      axis,xaxis=1,xrange=xrange,xtickformat='(A1)'
      axis,yaxis=0,yrange=yrange,ytitle='!9'+string(ntaletter)+'!x'
      axis,yaxis=1,yrange=yrange,ytickformat='(A1)'
   device,/close

   psname='nta_mass_function_hayward.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,mass,ntamgh,xtitle='Log(M/M'+sunsym+')',ytitle='!9'+string(ntaletter)+'!x',/ylog,/nodata,xrange=[6,12],xstyle=1,yrange=[0.01,500],ystyle=1
      x=[mass,reverse(mass)]
      ymg=[ntamghhi,reverse(ntamghlo)]
      polyfill,x,ymg,color=fsc_color('salmon')

      oplot,mass,ntamgh
      oplot,mass,10.^fitlogntamgh,color=fsc_color('darkred'),linestyle=5
   device,/close

end

pro writetable,sciall
    separate_cluster,sciall,wonsdss,woncl0024,wonms0451,cwonsdss,cwoncl0024,cwonms0451 
    cl0024 = sciall(woncl0024)
    ms0451 = sciall(wonms0451)
    wbad1 = where(strmatch(strtrim(cl0024.objname,2),'p30i4c3') eq 1,cbad1)
    ;find matching nta
    restore,'nta_data.sav'
    match,cl0024.objname,sci.objname,subcl,subsci,count=ctmatch
    cl0024_nta = fltarr(4,cwoncl0024)-999.
    cl0024_nta[0,subcl] = ntamgh(subsci)
    cl0024_nta[1:3,subcl] = ntalims[*,subsci]

    match,ms0451.objname,sci.objname,subcl,subsci,count=ctmatch
    ms0451_nta = fltarr(4,cwonms0451)-999.
    ms0451_nta[0,subcl] = ntamgh(subsci)
    ms0451_nta[1:3,subcl] = ntalims[*,subsci]
    if cbad1 eq 1 then begin
       cl0024(wbad1).ra= 6.6246667
       cl0024(wbad1).dec= 17.168556
    endif
    wbad2 = where(strmatch(strtrim(cl0024.objname,2),'p35i143c2') eq 1,cbad2)
    if cbad2 eq 1 then begin
       cl0024(wbad2).ra= 6.6552917
       cl0024(wbad2).dec=  17.166722
    endif
;    for jj=0,1 do begin
;       if jj eq 0 then sci = cl0024
;       if jj eq 1 then sci = ms0451
;       ra = sci.ra
;       ord = sort(ra)
;       ra = sci(ord).ra
;       dec = sci(ord).dec
;       logmstar = sci(ord).logmstar
;       feh = sci(ord).feh
;       fehupper = sci(ord).fehupper
;       fehlower = sci(ord).fehlower
;       age = sci(ord).age
;       ageupper = sci(ord).ageupper
;       agelower = sci(ord).agelower    
;       mgfe = sci(ord).alphafe
;       mgfeupper = sci(ord).alphafeupper
;       mgfelower = sci(ord).alphafelower
;       redshift = sci(ord).zfit
;       if jj eq 0 then begin
;          openw, lun2, 'table2_paper2.tex', /get_lun
;          printf, lun2, '\begin{deluxetable*}{llllCCC}'
;          printf, lun2, '\tablewidth{0pt}'
;          printf, lun2, '\tablecolumns{9}'
;          printf, lun2, '\tablecaption{Catalog of Measured Age and Metalicities\label{tab:catalog}}'
;          printf, lun2, '\tablehead{\colhead{No.} & \colhead{RA} & \colhead{DEC} & \colhead{$log(M_*/M_\odot)$} & \colhead{[Fe/H]} & \colhead{Age (Gyr)} & \colhead{S/N}'
;          printf, lun2, '\startdata'
;       endif
; 
;       for i=0,9 do begin
;           fes = string(feh[i],fehupper[i],fehlower[i], format='("$",D+5.2,"^{",D+5.2,"}_{",D+5.2,"}$")')
;           ages = string(age[i],ageupper[i],agelower[i],format='("$",D5.1,"^{",D5.1,"}_{",D5.1,"}$")')
;           mgfes = string(mgfe[i],mgfeupper[i],mgfelower[i],format='("$",D+5.2,"^{",D+5.2,"}_{",D+5.2,"}$")')
;           radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
;           ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
;           decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
;           mass = string(logmstar[i],format='(D4.1)')
;           if i ge 9 then nos = string(i+1,format='(I02)') else nos= string(i+1,format='(I01)')
;   
;           printf, lun2,nos,ras,decs,mass,ages,fes,mgfes, format='(A3,2(" & ",A11)," & ",A5,4(" & ",A24),"\\")'
;       endfor
;    endfor
;    close,lun2
;    free_lun, lun2

;    for jj=0,1 do begin
;       if jj eq 0 then sci = cl0024
;       if jj eq 1 then sci = ms0451
;       ra = sci.ra
;       ord = sort(ra)
;       ra = sci(ord).ra
;       dec = sci(ord).dec
;       logmstar = sci(ord).logmstar
;       feh = sci(ord).feh
;       fehupper = sci(ord).fehupper
;       fehlower = sci(ord).fehlower
;       age = sci(ord).age
;       ageupper = sci(ord).ageupper
;       agelower = sci(ord).agelower    
;       mgfe = sci(ord).alphafe
;       mgfeupper = sci(ord).alphafeupper
;       mgfelower = sci(ord).alphafelower
;       redshift = sci(ord).zfit
;       if jj eq 0 then begin
;          openw, lun2, 'table_thesis.tex', /get_lun
;          printf, lun2, '\begin{deluxetable*}{llllCCC}'
;          printf, lun2, '\tablewidth{0pt}'
;          printf, lun2, '\tablecolumns{9}'
;          printf, lun2, '\tablecaption{Catalog of Measured Age and Metalicities\label{tab:catalog}}'
;          printf, lun2, '\tablehead{\colhead{No.} & \colhead{RA} & \colhead{DEC} & \colhead{$log(M_*/M_\odot)$} & \colhead{[Fe/H]} & \colhead{Age (Gyr)} & \colhead{S/N}'
;          printf, lun2, '\startdata'
;       endif
; 
;       for i=0,n_elements(feh)-1 do begin
;           fes = string(feh[i],fehupper[i],fehlower[i], format='("$",D+5.2,"^{",D+5.2,"}_{",D+5.2,"}$")')
;           ages = string(age[i],ageupper[i],agelower[i],format='("$",D5.1,"^{",D5.1,"}_{",D5.1,"}$")')
;           mgfes = string(mgfe[i],mgfeupper[i],mgfelower[i],format='("$",D+5.2,"^{",D+5.2,"}_{",D+5.2,"}$")')
;           radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
;           ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
;           decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
;           mass = string(logmstar[i],format='(D4.1)')
;           if i ge 9 then nos = string(i+1,format='(I02)') else nos= string(i+1,format='(I01)')
;   
;           printf, lun2,nos,ras,decs,mass,ages,fes,mgfes, format='(A3,2(" & ",A11)," & ",A5,4(" & ",A24),"\\")'
;       endfor
;    endfor
;    close,lun2
;    free_lun, lun2

;    openw, lun2, 'table1_catalog.txt', /get_lun
;    printf, lun2, 'No.    RA    DEC    logM    FeH50   FeH84   FeH16    AgeGyr50   '+$
;                     'AgeGyr84 AgeGyr16 Mgfe50  MgFe84   Mgfe16'
;    for jj=0,1 do begin
;       if jj eq 0 then sci = cl0024
;       if jj eq 1 then sci = ms0451
;       ra = sci.ra
;       ord = sort(ra)
;       ra = sci(ord).ra
;       dec = sci(ord).dec
;       logmstar = sci(ord).logmstar
;       feh = sci(ord).feh
;       fehupper = sci(ord).fehupper
;       fehlower = sci(ord).fehlower
;       age = sci(ord).age
;       ageupper = sci(ord).ageupper
;       agelower = sci(ord).agelower
;       mgfe = sci(ord).alphafe
;       mgfeupper = sci(ord).alphafeupper
;       mgfelower = sci(ord).alphafelower
;       redshift = sci(ord).zfit
;
;       n = n_elements(ra)
;       for i=0, n-1 do begin
;           radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
;           ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
;           decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
;           printf,lun2,i+1,ras,decs,logmstar[i],feh[i],fehupper[i],fehlower[i],age[i],ageupper[i],agelower[i],$
;                  mgfe[i],mgfeupper[i],mgfelower[i],$
;                  format='(I02,1x,A,1x,A,1x,D4.1,1x,D+5.2,1x,D+5.2,1x,D+5.2,1x,D5.1,1x,D5.1,1x,D5.1,1x,D+5.2,1x,D+5.2,1x,D+5.2)'
;       endfor          
;    endfor
;    close,lun2
;    free_lun, lun2
;    stop

    openw, lun2, 'table1_catalog_revised.txt', /get_lun
    printf, lun2, 'No.    RA    DEC    logM    FeH50   FeH84   FeH16    AgeGyr50   '+$
                     'AgeGyr84 AgeGyr16 Mgfe50  MgFe84   Mgfe16  eta  eta84  eta16'
    for jj=0,1 do begin
       if jj eq 0 then begin
          sci = cl0024
          nta = cl0024_nta
       endif
       if jj eq 1 then begin
          sci = ms0451
          nta = ms0451_nta
       endif
       ra = sci.ra
       ord = sort(ra)
       ra = sci(ord).ra
       dec = sci(ord).dec
       logmstar = sci(ord).logmstar
       feh = sci(ord).feh
       fehupper = sci(ord).fehupper
       fehlower = sci(ord).fehlower
       age = sci(ord).age
       ageupper = sci(ord).ageupper
       agelower = sci(ord).agelower
       mgfe = sci(ord).alphafe
       mgfeupper = sci(ord).alphafeupper
       mgfelower = sci(ord).alphafelower
       redshift = sci(ord).zfit
       ntabest = nta[0,ord]
       ntaupper = nta[3,ord]
       ntalower = nta[1,ord]
       n = n_elements(ra)
       for i=0, n-1 do begin
           radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
           ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
           decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
           printf,lun2,i+1,ras,decs,logmstar[i],feh[i],fehupper[i],fehlower[i],age[i],ageupper[i],agelower[i],$
                  mgfe[i],mgfeupper[i],mgfelower[i],ntabest[i],ntaupper[i],ntalower[i],$
                  format='(I02,1x,A,1x,A,1x,D4.1,1x,D+5.2,1x,D+5.2,1x,D+5.2,1x,D5.1,1x,D5.1,1x,D5.1,1x,D+5.2,1x,D+5.2,1x,D+5.2,1x,D+5.2,1x,D+5.2,1x,D+5.2)'
       endfor          
    endfor
    close,lun2
stop
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;MAIN PROGRAM;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro all_alphaage_ana,domgo=domgo
common filenames, aveparam_file,linfitparam_file,linfitparam_fixslope_file,fileanova,$
                fileancova,deviationparam_file

common color_pallete, clustercolor, shadecolor, meancolor
   redodatamgn = 0
   redodatamgo = 0
   redoaverage = 0
   redolinearfit = 0
   retest_evolution_regression = 0 ;old method of ANCOVA
   redoancova = 0
   redoanova = 0
   redodeviation_fromz0line = 0
   redocalculate_nta = 0

   doplot_masshist = 0
   doplot_mzr = 0
   doplot_byage = 0
   doplot_deviation = 0
   doplot_alphafe = 0
   dowritetable = 1

   clustercolor = ['skyblue','goldenrod','crimson']
   shadecolor = ['honeydew','wt3','red1']
   meancolor=['navyblue','peru','indianred']
;set up
   if ~keyword_set(domgo) then begin
      domgn=1
      domgo=0
      prefix='mgn'
   endif else begin
      domgn=0
      domgo=1
      prefix='mgo'
   endelse

   if redodatamgn or redodatamgo then begin
      redolinearfit = 1
      redoaverage = 1
      redo_evolution_regression = 1
    endif

;file names
   aveparam_file = 'aveparam_'+prefix+'.sav'
   linfitparam_file = 'linfitpar_'+prefix+'.sav'
   linfitparam_fixslope_file = 'linfitpar_fixslope_'+prefix+'.sav'
   fileancova = 'ancova_outstr_'+prefix+'.fits'
   fileanova = 'anova_evolution_outstr_'+prefix+'.fits'
;   deviationparam_file = 'deviationparam_'+prefix+'.sav'
   deviationparam_file = 'deviationparam_mgn_realz0.sav' ;model with interaction terms
;  deviationparam_file = 'deviationparam_mgn_model1.sav' ;model without interaction terms
;   deviationparam_file = 'deviationparam_mgn_rmmass.sav' ;model with interaction terms, remove all mass
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;PREP DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if redodatamgn eq 1 then begin
      outfile = '/scr2/nichal/workspace4/ana/all/all_sci_mgn.fits'
      sdssfile= '/scr2/nichal/workspace4/sps_fit/data/gallazzi_allmass2/sps_fit03.fits.gz'
      sdssprob= '/scr2/nichal/workspace4/sps_fit/data/gallazzi_allmass2/probdist/probdist_03.fits'
      cl0024file = '/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit03.fits.gz'
      cl0024prob = '/scr2/nichal/workspace4/sps_fit/data/cl0024/probdist/probdist_03.fits'
      ms0451file = '/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz'
      ms0451prob = '/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/probdist/probdist_03.fits'
      prepdata,sdssfile=sdssfile,cl0024file=cl0024file,ms0451file=ms0451file,$
             sdssprob=sdssprob,cl0024prob=cl0024prob,ms0451prob=ms0451prob,outfile=outfile
   endif 

   if redodatamgo eq 1 then begin
      outfile = '/scr2/nichal/workspace4/ana/all/all_sci_mgo.fits'
      sdssfile= '/scr2/nichal/workspace4/sps_fit/data/gallazzi_allmass2/sps_fit06.fits.gz'
      sdssprob= '/scr2/nichal/workspace4/sps_fit/data/gallazzi_allmass2/probdist/probdist_06.fits'
      cl0024file = '/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit06.fits.gz'
      cl0024prob = '/scr2/nichal/workspace4/sps_fit/data/cl0024/probdist/probdist_06.fits'
      ms0451file = '/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit06.fits.gz'
      ms0451prob = '/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/probdist/probdist_06.fits'
      prepdata,sdssfile=sdssfile,cl0024file=cl0024file,ms0451file=ms0451file,$
             sdssprob=sdssprob,cl0024prob=cl0024prob,ms0451prob=ms0451prob,outfile=outfile
   endif 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;READ DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if domgo then outfile = '/scr2/nichal/workspace4/ana/all/all_sci_mgo.fits'
   if domgn then outfile = '/scr2/nichal/workspace4/ana/all/all_sci_mgn.fits'

   sciall = mrdfits(outfile,1)
   proball = mrdfits(outfile,2)
   scisdss = mrdfits(outfile,3)
   scicl0024 = mrdfits(outfile,4)
   scims0451 = mrdfits(outfile,5)

   ;redo the goodsn cut for cl0024
   newgood = where(sciall.zfit gt 0.3 and sciall.zfit lt 0.45 and sciall.snfit gt 5.5)
   sciall(newgood).goodsn = 1
;   sciall([185,213,178,195]).goodsn = 0
   sciall([185,213,195]).goodsn = 0
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;average and get the highlight
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if redoaverage eq 1 then average_mzr,sciall,proball,prefix=prefix
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;linear fit to all data
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if redolinearfit eq 1 then linearfit_mzr,sciall,proball,prefix=prefix 
   restore,linfitparam_file 
   print, ':::FITTING WITH LINEAR FUNCTION(MONTE CARLO):::'
   print, 'MZR PARAMETERS'
   print, 'cluster       intercept         slope         chisq'
   print, 'ALL     '+string(allmzrpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(allmzrerr[0],format='(f5.2)'),2)+$
          '        '+string(allmzrpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(allmzrerr[1],format='(f5.2)'),2)+$
          '        '+string(chiallmzr,format='(f5.2)')
   print, 'SDSS    '+string(sdssmzrpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(sdssmzrerr[0],format='(f5.2)'),2)+$
          '        '+string(sdssmzrpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(sdssmzrerr[1],format='(f5.2)'),2)+$
          '        '+string(chisdssmzr,format='(f5.2)')
   print, 'Cl0024  '+string(cl0024mzrpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(cl0024mzrerr[0],format='(f5.2)'),2)+$
          '        '+string(cl0024mzrpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(cl0024mzrerr[1],format='(f5.2)'),2)+$
          '        '+string(chicl0024mzr,format='(f5.2)')
   print, 'MS0451  '+string(ms0451mzrpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(ms0451mzrerr[0],format='(f5.2)'),2)+$
          '        '+string(ms0451mzrpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(ms0451mzrerr[1],format='(f5.2)'),2)+$
          '        '+string(chims0451mzr,format='(f5.2)')
   print, 'MAR PARAMETERS ([Mg/H])'
   print, 'cluster       intercept         slope'
   print, 'ALL     '+string(allmarpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(allmarerr[0],format='(f5.2)'),2)+$
          '        '+string(allmarpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(allmarerr[1],format='(f5.2)'),2)+$
          '        '+string(chiallmar,format='(f5.2)')
   print, 'SDSS    '+string(sdssmarpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(sdssmarerr[0],format='(f5.2)'),2)+$
          '        '+string(sdssmarpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(sdssmarerr[1],format='(f5.2)'),2)+$
          '        '+string(chisdssmar,format='(f5.2)')
   print, 'Cl0024  '+string(cl0024marpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(cl0024marerr[0],format='(f5.2)'),2)+$
          '        '+string(cl0024marpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(cl0024marerr[1],format='(f5.2)'),2)+$
          '        '+string(chicl0024mar,format='(f5.2)')
   print, 'MS0451  '+string(ms0451marpar[0],format='(f5.2)')+'+/-'+$
                    strtrim(string(ms0451marerr[0],format='(f5.2)'),2)+$
          '        '+string(ms0451marpar[1],format='(f5.2)')+'+/-'+$
                    strtrim(string(ms0451marerr[1],format='(f5.2)'),2)+$
          '        '+string(chims0451mar,format='(f5.2)')
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   ;TEST EVOLUTION;;;;;;;;;;;;;;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   ;FIRST, TEST WITH REGRESSION
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if retest_evolution_regression eq 1 then test_evolution_regression,sciall,proball,prefix=prefix 
   restore,linfitparam_fixslope_file
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, ':::TEST EVOLUTION WITH ANCOVA (OLD METHOD):::'
   print, 'FEH:'
   print, 'normal t value for slope of SDSS and Cl0024 =',tslope_sdcl_feh
   print, 'normal t value for slope of SDSS and MS0451 =',tslope_sdms_feh
   print, 'normal t value for slope of Cl0024 and MS0451 =',tslope_clms_feh
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'MgH:'
   print, 'normal t value for slope of SDSS and Cl0024 =',tslope_sdcl_mgh
   print, 'normal t value for slope of SDSS and MS0451 =',tslope_sdms_mgh
   print, 'normal t value for slope of Cl0024 and MS0451 =',tslope_clms_mgh

   ;use t values from here https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
   ;https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'Two sample t-test'
   print, 'Reject the null hypothesis that the two pop have the same slope if t values is larger than'
   print, '(for large degree of freedom but less than 100)'
   print, '2 sigma: |t| > 2'
   print, '3 sigma: |t| > 3.2'
   ;if cannot reject the null hypothesis then proceed to find common slope
   ;http://www.biostathandbook.com/ancova.html  
   ;http://www.stat.ucla.edu/~cochran/stat10/winter/lectures/lect21.html
   ;find the weighted slope
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print,'Common Slopes are: commonslope, uncertainty, standard deviation'
   print, '   For [Fe/H] = ',fehcommonslope,fehcommonslope_err,fehcommonslope_stdev
   print, '   For [Mg/H] = ',mghcommonslope,mghcommonslope_err,mghcommonslope_stdev
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'Linear fit with common fixed slope:'
   print, 'FeH'
   print, 'cluster    slope       intercept      intercepterr   reduced chisq'
   print, 'SDSS ',sdssmzrpar_fixslope[1], sdssmzrpar_fixslope[0], sdssmzrerr_fixslope[0],sdssmzrchi_fixslope
   print, 'Cl0024 ',cl0024mzrpar_fixslope[1], cl0024mzrpar_fixslope[0], cl0024mzrerr_fixslope[0],cl0024mzrchi_fixslope
   print, 'MS0451 ',ms0451mzrpar_fixslope[1], ms0451mzrpar_fixslope[0], ms0451mzrerr_fixslope[0],ms0451mzrchi_fixslope
   print, 'MgH'
   print, 'cluster    slope       intercept      intercepterr    reduced chisq'
   print, 'SDSS ',sdssmarpar_fixslope[1], sdssmarpar_fixslope[0], sdssmarerr_fixslope[0],sdssmarchi_fixslope
   print, 'Cl0024 ',cl0024marpar_fixslope[1], cl0024marpar_fixslope[0], cl0024marerr_fixslope[0],cl0024marchi_fixslope
   print, 'MS0451 ',ms0451marpar_fixslope[1], ms0451marpar_fixslope[0], ms0451marerr_fixslope[0],ms0451marchi_fixslope
   ;t test comparing two intercepts with fixed slope
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'FEH:'
   print, 'normal t value for intercept of SDSS and Cl0024 =',tconst_sdcl_feh
   print, 'normal t value for intercept of SDSS and MS0451 =',tconst_sdms_feh
   print, 'normal t value for intercept of Cl0024 and MS0451 =',tconst_clms_feh
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'MgH:'
   print, 'normal t value for intercept of SDSS and Cl0024 =',tconst_sdcl_mgh
   print, 'normal t value for intercept of SDSS and MS0451 =',tconst_sdms_mgh
   print, 'normal t value for intercept of Cl0024 and MS0451 =',tconst_clms_mgh

   ;ALSO TEST WITH THE FULL METHOD OF ANCOVA
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, ':::TEST EVOLUTION WITH ANCOVA (FULL METHOD):::'
   if redoancova then ancova,sciall,proball,fileancova=fileancova
   restore, fileancova   
;   goodmass = where(sciall.logmstar gt 9.6)
;   sciall = sciall(goodmass)
;   proball = proball(goodmass)
;   if redoancova then ancova,sciall,proball,fileancova='ancova_outstr_limmass9.6.fits'
;   restore, 'ancova_outstr_limmass9.6.fits'
   fchange_mo1_feh = (rsq_feh_mo1/df1_mo1)/((1.-rsq_feh_mo1)/df2_mo1)
   fchange_mo2_feh = ((rsq_feh_mo2-rsq_feh_mo1)/(df1_mo2-df1_mo1))/((1.-rsq_feh_mo2)/df2_mo2)
   dummy = mpftest(fchange_mo1_feh,df1_mo1,df2_mo1)
   dummy = t_pdf(0.95,86)
   print, ';;;;;;;;;;;;;;;[Fe/H];;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'model    R^2     adjusted R^2      std. error of the estimate'
   print, '1    ',rsq_feh_mo1,adjrsq_feh_mo1,stderror_feh_mo1
   print, '2    ',rsq_feh_mo2,adjrsq_feh_mo2,stderror_feh_mo2
   print, ''
   print, '        Change Statistics'
   print, 'model    R^2 change    F change    df1    df2      Sig F change'
   print, '1    ',rsq_feh_mo1,fchange_mo1_feh,df1_mo1,df2_mo1,mpftest(fchange_mo1_feh,df1_mo1,df2_mo1)
   print, '2    ',rsq_feh_mo2-rsq_feh_mo1,fchange_mo2_feh,df1_mo2-df1_mo1,df2_mo2,$
                mpftest(fchange_mo2_feh,df1_mo2-df1_mo1,df2_mo2)
   print, ''
   print, 'Model Parameters'
   print, 'model  parameter     value     std.error      t      sig.'
   print, '1    ','constant',const_feh_mo1,consterr_feh_mo1,const_feh_mo1/consterr_feh_mo1,$
                           (1.-t_pdf(abs(const_feh_mo1/consterr_feh_mo1),df2_mo1)) 
   for i=0,2 do $
   print, '     ',mo1_pname[i],coeff_feh_mo1[i],sigma_feh_mo1[i],coeff_feh_mo1[i]/sigma_feh_mo1[i],$
                (2*(1.-t_pdf(abs(coeff_feh_mo1[i]/sigma_feh_mo1[i]),df2_mo1)))
   print, '2    ','constant   ',const_feh_mo2,consterr_feh_mo2,const_feh_mo2/consterr_feh_mo2,$
                           (1.-t_pdf(abs(const_feh_mo2/consterr_feh_mo2),df2_mo2))
   for i=0,4 do $
   print, '     ',mo2_pname[i],coeff_feh_mo2[i],sigma_feh_mo2[i],coeff_feh_mo2[i]/sigma_feh_mo2[i],$
                (1.-t_pdf(abs(coeff_feh_mo2[i]/sigma_feh_mo2[i]),df2_mo2))
   ;one tail test (the effect is only 1 direction)
   fchange_mo1_mgh = (rsq_mgh_mo1/df1_mo1)/((1.-rsq_mgh_mo1)/df2_mo1)
   fchange_mo2_mgh = ((rsq_mgh_mo2-rsq_mgh_mo1)/(df1_mo2-df1_mo1))/((1.-rsq_mgh_mo2)/df2_mo2)
   fchange_mo3_mgh = ((rsq_mgh_mo1-rsq_mgh_mo3)/(df1_mo1-df1_mo3))/((1.-rsq_mgh_mo1)/df2_mo1)
   print, ';;;;;;;;;;;;;;;[Mg/H];;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, 'model    R^2     adjusted R^2      std. error of the estimate'
   print, '1    ',rsq_mgh_mo1,adjrsq_mgh_mo1,stderror_mgh_mo1
   print, '2    ',rsq_mgh_mo2,adjrsq_mgh_mo2,stderror_mgh_mo2
   print, '3    ',rsq_mgh_mo3,adjrsq_mgh_mo3,stderror_mgh_mo3
   
   print, ''
   print, '        Change Statistics'
   print, 'model         R^2 change    F change    df1    df2      Sig F change'
   print, '1 to none',rsq_mgh_mo1,fchange_mo1_mgh,df1_mo1,df2_mo1,mpftest(fchange_mo1_mgh,df1_mo1,df2_mo1)
   print, '2 to 1   ',rsq_mgh_mo2-rsq_mgh_mo1,fchange_mo2_mgh,df1_mo2-df1_mo1,df2_mo2,$
                mpftest(fchange_mo2_mgh,df1_mo2-df1_mo1,df2_mo2)
   print, '1 to 3   ',rsq_mgh_mo1-rsq_mgh_mo3,fchange_mo3_mgh,df1_mo1-df1_mo3,df2_mo1,$
                mpftest(fchange_mo3_mgh,df1_mo1-df1_mo3,df2_mo1)
   print, ''
   print, 'Model Parameters'
   print, 'model  parameter     value     std.error      t      sig.'
   print, '1    ','constant',const_mgh_mo1,consterr_mgh_mo1,const_mgh_mo1/consterr_mgh_mo1,$
                           (1.-t_pdf(abs(const_mgh_mo1/consterr_mgh_mo1),df2_mo1))
   for i=0,2 do $
   print, '     ',mo1_pname[i],coeff_mgh_mo1[i],sigma_mgh_mo1[i],coeff_mgh_mo1[i]/sigma_mgh_mo1[i],$
                (2*(1.-t_pdf(abs(coeff_mgh_mo1[i]/sigma_mgh_mo1[i]),df2_mo1)))
   print, '2    ','constant  ',const_mgh_mo2,consterr_mgh_mo2,const_mgh_mo2/consterr_mgh_mo2,$
                           (1.-t_pdf(abs(const_mgh_mo2/consterr_mgh_mo2),df2_mo2))
   for i=0,4 do $
   print, '     ',mo2_pname[i],coeff_mgh_mo2[i],sigma_mgh_mo2[i],coeff_mgh_mo2[i]/sigma_mgh_mo2[i],$
                (1.-t_pdf(abs(coeff_mgh_mo2[i]/sigma_mgh_mo2[i]),df2_mo2))
   print, '3    ','constant  ',const_mgh_mo3,consterr_mgh_mo3,const_mgh_mo3/consterr_mgh_mo3,$
                           (1.-t_pdf(abs(const_mgh_mo3/consterr_mgh_mo3),df2_mo3))
   print, '     ',mo3_pname[0],coeff_mgh_mo3[0],sigma_mgh_mo3[0],coeff_mgh_mo3[0]/sigma_mgh_mo3[0],$
                (1.-t_pdf(abs(coeff_mgh_mo3[0]/sigma_mgh_mo3[0]),df2_mo3))
   ;SECOND, TEST WITH ANOVA FOR EACH MASS BIN
   ;group into mass bins then compare the three redshifts from each mass bin with ANOVA
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if redoanova then anova_massbin,sciall,proball,prefix=prefix,fileanova=fileanova
   anova_feh = mrdfits(fileanova,1)
   anova_mgh = mrdfits(fileanova,2)
  
   print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print, ':::TEST EVOLUTION IN EACH MASS BIN WITH ANOVA:::'
   for i=0,n_elements(anova_feh)-1 do begin
      if anova_feh[i].mbinmin ne anova_mgh[i].mbinmin then stop, 'might have to redo anova, STOPPED'
      print, 'mass bin',anova_feh[i].mbinmin,anova_feh[i].mbinmax
      print, 'FeH ANOVA'
      print, 'TSS  =   SSB    +     SSW'
      print, anova_feh[i].TSS,'=?',anova_feh[i].SSB,'+',anova_feh[i].SSW
      print, anova_feh[i].TSS,'=?',anova_feh[i].SSB+anova_feh[i].SSW
      print, 'f('+strtrim(anova_feh[i].dof_ssb,2)+','+strtrim(anova_feh[i].dof_ssw,2)+$
             ') = ',anova_feh[i].f_anova
      print, 'MgH ANOVA'
      print, 'TSS  =   SSB    +     SSW'
      print, anova_mgh[i].TSS,'=?',anova_mgh[i].SSB,'+',anova_mgh[i].SSW
      print, anova_mgh[i].TSS,'=?',anova_mgh[i].SSB+anova_mgh[i].SSW
      print, 'f('+strtrim(anova_mgh[i].dof_ssb,2)+','+strtrim(anova_mgh[i].dof_ssw,2)+$
             ') = ',anova_mgh[i].f_anova
      print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   endfor

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;DEVIATION FROM Z0 LINE
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   fehcoeff = [const_feh_mo2,coeff_feh_mo2[0],coeff_feh_mo2[3],coeff_feh_mo2[4]] ;(rmmass)
;   fehcoeff = [const_feh_mo1,coeff_feh_mo2[0],0.,0.] ;model with interaction terms (realz0)
   fehcoeff = [const_feh_mo1,coeff_feh_mo1[0],0.,0.] ;model with no interaction terms (model1)
   mghcoeff = [const_mgh_mo3,coeff_mgh_mo3[0]]
   if redodeviation_fromz0line then deviation_fromz0line,sciall,proball,fehcoeff,$
                                    mghcoeff,fileout=deviationparam_file
   print, 'RESTORING ',deviationparam_file
   restore,deviationparam_file
   print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   print,'DEVIATION FROM Z0 LINE'
   print,'FeH:'
   print,'intercept, slope'
   print,'probgrid fit:',feh_dev_par_all
   print,'mc fit',mc_feh_dev_par,mc_feh_dev_par_err
   print,'simple fit',simple_feh_dev_par,simple_mgh_dev_par_err

   print,'MgH:'      
   print,'intercept, slope'
   print,'mc fit',mc_mgh_dev_par,mc_mgh_dev_par_err
   print,'simple fit',simple_mgh_dev_par,simple_mgh_dev_par_err

   wonsdss = where(sciall.zfit lt 0.3 and sciall.goodsn eq 1 and sciall.ageform gt 0.5)
   woncl0024 = where(sciall.zfit gt 0.3 and sciall.zfit lt 0.45 and sciall.goodsn eq 1 and sciall.ageform gt 0.5)
   wonms0451 = where(sciall.zfit gt 0.45 and sciall.goodsn eq 1 and sciall.ageform gt 0.5)
   meanage_sdss = mean(sciall(wonsdss).ageform)
   meanage_cl0024 = mean(sciall(woncl0024).ageform)
   meanage_ms0451 = mean(sciall(wonms0451).ageform)
   print, 'sample          mean age form          (mean age-SDSS age)*Feh_slope     (meanage-SDSS)*Mgh slope'
   print, 'SDSS',meanage_sdss
   print, 'Cl0024',meanage_cl0024,(meanage_cl0024-meanage_sdss)*feh_dev_par[1],(meanage_cl0024-meanage_sdss)*simple_mgh_dev_par[1]
   print, 'MS0451',meanage_ms0451,(meanage_ms0451-meanage_sdss)*feh_dev_par[1],(meanage_ms0451-meanage_sdss)*simple_mgh_dev_par[1]

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;CALCULATE MASS LOADING FACTOR
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if redocalculate_nta then begin
       coeff = [const_mgh_mo3,coeff_mgh_mo3[0]]
       coefferr = [consterr_mgh_mo3,sigma_mgh_mo3[0]]
       calculate_nta_function,coeff,coefferr
       calculate_nta_data,sciall,proball,coeff,coefferr
   endif
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if doplot_masshist then plot_masshist,sciall
   if doplot_mzr then plot_mzr,sciall,prefix=prefix
   if doplot_byage then plot_byage,sciall,prefix=prefix
   if doplot_deviation then plot_deviation_fromz0line,sciall
   if doplot_alphafe then plot_feh_alpha_byage,sciall
   if dowritetable then writetable,sciall 
stop
end
