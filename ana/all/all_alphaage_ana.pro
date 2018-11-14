pro all_alphaage_ana

   redodata = 0
   if redodata eq 1 then begin
      scitemp = {objname:'',ra:99.,dec:99.,snfit:99.,zspec:99.,zfit:99.,z:99.,age:99.,$
                 ageerr:99.,ageupper:99.,agelower:99.,feh:99.,feherr:99.,$
                 fehupper:99.,fehlower:99.,alphafe:99.,alphafeerr:99.,$
                 alphafeupper:99.,alphafelower:99.,logmstar:99.,oiiew:99.,$
                 oiiewerr:99.,vdisp:99.,vdisperr:99.,pfit:fltarr(13),pfiterr:fltarr(13),$
                 ah:99.,ahupper:99.,ahlower:99,goodsn:1b,ageform:99.} 

      ;GETTING SCI DATA
      sciall=[]
      goodspecall=[]
      proball = []
 
      ;SDSS
      scinow = mrdfits('/scr2/nichal/workspace4/sps_fit/data/gallazzi_allmass2/sps_fit03.fits.gz',1)
      probnow = mrdfits('/scr2/nichal/workspace4/sps_fit/data/gallazzi_allmass2/probdist/probdist_03.fits',1)
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
      scinow = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/sps_fit03.fits.gz',1)
      probnow = mrdfits('/scr2/nichal/workspace4/sps_fit/data/cl0024/probdist/probdist_03.fits',1)
      ;sample selection
      good = where(scinow.good eq 1 and scinow.oiiew gt -5 and scinow.snfit gt 2. and scinow.logmstar gt 6.,cgood)
      scinow = scinow(good)
      scicl0024 = scinow
      probnow = probnow(good)
      goodfit = where(scinow.snfit gt 6.83, cgoodfit, complement=badfit,ncomplement=cbadfit)
      goodspec = bytarr(cgood)
      goodspec(goodfit) = 1

      str = replicate(scitemp,cgood)
      struct_assign,scinow,str
      sciall = [sciall,str]
      proball = [proball,probnow]
      goodspecall = [goodspecall,goodspec]

      ;ms0451
      scinow = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
      probnow = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/probdist/probdist_03.fits',1)
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
      sciall.ageform = (galage(sciall.zfit,1000.)/1.e9-sciall.age)>0. ;age of universe when it was formed
      sciall.ah = sciall.feh+sciall.alphafe
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;MILES Library abundance (from Conroy17)
      miles = {feh:[-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],$
               ofe:[0.6,0.5,0.5,0.4,0.3,0.2,0.2,0.1,0.0,0.0],$
               mgfe:[0.4,0.4,0.4,0.4,0.34,0.22,0.14,0.11,0.05,0.04],$
               cafe:[0.32,0.30,0.28,0.26,0.26,0.17,0.12,0.06,0.00,0.00]}
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      ;fix the library abundance
       sciall.alphafe = sciall.alphafe+interpol(miles.mgfe,miles.feh,sciall.feh)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
         toplot = 0
         if toplot eq 1 then begin
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
;            plothist,griddah(loc),xtitle='delta [Mg/H]',position=[0.65,0.6,0.95,0.95],/noerase
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
;         print, int_tabulated_2d(griddfeh,griddage,probfehage[*,*,nn])
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

      outfile = 'all_sci.fits'
      mwrfits, sciall, outfile,/create,/silent
      mwrfits, proball,outfile,/silent
      mwrfits, scisdss,outfile,/silent
      mwrfits,scicl0024,outfile,/silent
      mwrfits,scims0451,outfile,/silent
   endif else begin
      outfile = '/scr2/nichal/workspace4/ana/all/all_sci.fits'
      sciall = mrdfits(outfile,1)
      proball = mrdfits(outfile,2)
      scisdss = mrdfits(outfile,3)
      scicl0024 = mrdfits(outfile,4)
      scims0451 = mrdfits(outfile,5)
   endelse
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;linear fit to all data
   redolinearfit = 0;1
   if redolinearfit eq 1 then begin 
      sciall([185,213]).goodsn = 0
      won = where(sciall.goodsn eq 1)
     
      allmzrpar = linfitmc(sciall[won].logmstar-10.,sciall[won].feh,proball[won].dfeharr,proball[won].probfeh,$
                           allmzrerr,chisq=chiallmzr);,/showplot)
      allmarpar = linfitmc(sciall[won].logmstar-10.,sciall[won].ah,proball[won].daharr,proball[won].probdah,$
                           allmarerr,chisq=chiallmar);,/showplot)
      won = where(sciall.zfit lt 0.15 and sciall.goodsn eq 1)
      sdssmzrpar = linfitmc(sciall[won].logmstar-10.,sciall[won].feh,proball[won].dfeharr,proball[won].probfeh,$
                           sdssmzrerr,chisq=chisdssmzr);,/showplot)
      sdssmarpar = linfitmc(sciall[won].logmstar-10.,sciall[won].ah,proball[won].daharr,proball[won].probdah,$
                           sdssmarerr,chisq=chisdssmar);,/showplot)
      won = where(sciall.zfit gt 0.35 and sciall.zfit lt 0.42 and sciall.goodsn eq 1)
      cl0024mzrpar = linfitmc(sciall[won].logmstar-10.,sciall[won].feh,proball[won].dfeharr,proball[won].probfeh,$
                           cl0024mzrerr,chisq=chicl0024mzr);,/showplot)
      cl0024marpar = linfitmc(sciall[won].logmstar-10.,sciall[won].ah,proball[won].daharr,proball[won].probdah,$
                           cl0024marerr,chisq=chicl0024mar);,/showplot)
      won = where(sciall.zfit gt 0.48 and sciall.goodsn eq 1)
      ms0451mzrpar = linfitmc(sciall[won].logmstar-10.,sciall[won].feh,proball[won].dfeharr,proball[won].probfeh,$
                           ms0451mzrerr,chisq=chims0451mzr);,/showplot)
      ms0451marpar = linfitmc(sciall[won].logmstar-10.,sciall[won].ah,proball[won].daharr,proball[won].probdah,$
                           ms0451marerr,chisq=chims0451mar);,/showplot)
      print, 'MZR PARAMETERS'
      print, 'cluster       intercept         slope'
      print, 'ALL     '+string(allmzrpar[0],format='(f5.2)')+'+/-'+strtrim(string(allmzrerr[0],format='(f5.2)'),2)+$
             '        '+string(allmzrpar[1],format='(f5.2)')+'+/-'+strtrim(string(allmzrerr[1],format='(f5.2)'),2)
      print, 'SDSS    '+string(sdssmzrpar[0],format='(f5.2)')+'+/-'+strtrim(string(sdssmzrerr[0],format='(f5.2)'),2)+$
             '        '+string(sdssmzrpar[1],format='(f5.2)')+'+/-'+strtrim(string(sdssmzrerr[1],format='(f5.2)'),2)
      print, 'Cl0024  '+string(cl0024mzrpar[0],format='(f5.2)')+'+/-'+strtrim(string(cl0024mzrerr[0],format='(f5.2)'),2)+$
             '        '+string(cl0024mzrpar[1],format='(f5.2)')+'+/-'+strtrim(string(cl0024mzrerr[1],format='(f5.2)'),2)
      print, 'MS0451  '+string(ms0451mzrpar[0],format='(f5.2)')+'+/-'+strtrim(string(ms0451mzrerr[0],format='(f5.2)'),2)+$
             '        '+string(ms0451mzrpar[1],format='(f5.2)')+'+/-'+strtrim(string(ms0451mzrerr[1],format='(f5.2)'),2)
      print, 'MAR PARAMETERS'
      print, 'cluster       intercept         slope'
      print, 'ALL     '+string(allmarpar[0],format='(f5.2)')+'+/-'+strtrim(string(allmarerr[0],format='(f5.2)'),2)+$
             '        '+string(allmarpar[1],format='(f5.2)')+'+/-'+strtrim(string(allmarerr[1],format='(f5.2)'),2)
      print, 'SDSS    '+string(sdssmarpar[0],format='(f5.2)')+'+/-'+strtrim(string(sdssmarerr[0],format='(f5.2)'),2)+$
             '        '+string(sdssmarpar[1],format='(f5.2)')+'+/-'+strtrim(string(sdssmarerr[1],format='(f5.2)'),2)
      print, 'Cl0024  '+string(cl0024marpar[0],format='(f5.2)')+'+/-'+strtrim(string(cl0024marerr[0],format='(f5.2)'),2)+$
             '        '+string(cl0024marpar[1],format='(f5.2)')+'+/-'+strtrim(string(cl0024marerr[1],format='(f5.2)'),2)
      print, 'MS0451  '+string(ms0451marpar[0],format='(f5.2)')+'+/-'+strtrim(string(ms0451marerr[0],format='(f5.2)'),2)+$
             '        '+string(ms0451marpar[1],format='(f5.2)')+'+/-'+strtrim(string(ms0451marerr[1],format='(f5.2)'),2)
      save,allmzrpar,allmarpar,sdssmzrpar,sdssmarpar,cl0024mzrpar,cl0024marpar,ms0451mzrpar,ms0451marpar,filename='linfitpar.sav'
   endif else restore,'linfitpar.sav'
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   ;REFERENCE VALUES
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   restore,'/scr2/nichal/workspace4/ana/ms0451/leetho18a_avevals.sav'
   hifeh_cl0024 = hifeh
   lofeh_cl0024 = lofeh
   bndry_mass_cl0024 = bndry_mass

   ;Get the average values measured in Gallazzi05
   g05_mass = [9.00,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5];the first mass is actually 8.91
   g05_feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
   ;below are the stars in the figure 8 of Gallazzi05
   g05_feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
   g05_fehlo = g05_feh-g05_feherr
   g05_fehhi = g05_feh+g05_feherr
   toolow = where(g05_fehlo lt -1.,ctoolow)
   if ctoolow gt 0 then g05_fehlo(toolow) = -1

   ;Choi's Data
   Choi14z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
   Choi14z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
   Choi14z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
   Choi14z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
   Choi14z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}

   ;Xiangcheng's FIRE data
   readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0pt8.txt',xma_mass08,xma_feh08
   readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0.txt',xma_mass0,xma_feh0
   xma_feh08 = xma_feh08-0.2
   xma_feh0 = xma_feh0-0.2
   
   ;De Rossi 2017 (EAGLE)
   ;read off from Figure 5
   DeRossi_z0={mass:[9.15,9.48,9.81,10.15,10.5,10.72],feh:[-0.17,-0.08,0.05,0.12,0.17,0.31],feherr:[0.075,0.08,0.07,0.07,0.07,0.03]}
   DeRossi_z1={mass:[9.18,9.47,9.85,10.15,10.55],feh:[-0.39,-0.28,-0.09,0.08,0.17],feherr:[0.06,0.06,0.08,0.07,0.085]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Sybilska et al 2017 (hELENa, IFU from Sauron data)
   aa=read_csv('/scr2/nichal/workspace2/catalogs/sybilska.csv',n_table_header=1,header=header)
   Syb_z0 = {mass:aa.field1,feh:aa.field2}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Lu 2014
   Lu_z0 = {mass:[8.21,8.61,9.0,9.41,9.8,10.2,10.61,11],feh:[-0.53,-0.43,-0.35,-0.26,-0.16,-0.098,-0.075,-0.09],$
            feherr:[0.13,0.12,0.11,0.12,0.11,0.10,0.095,0.095]}
   Lu_z1 = {mass:[8.21,8.61,9,9.4,9.8,10.2,10.6,11],feh:[-0.56,-0.44,-0.35,-0.25,-0.18,-0.13,-0.11,-0.11],$
             feherr:[0.14,0.13,0.11,0.11,0.11,0.10,0.10,0.09]}
   Lu_somerville_z0 = {mass:[8.25,8.65,9.05,9.45,9.85,10.25,10.64],feh:[-1.00,-0.82,-0.64,-0.46,-0.33,-0.15,-0.03],$
                      feherr:[0.06,0.1,0.1,0.1,0.1,0.10,0.095]}
   Lu_lu_z0 = {mass:[8.73,8.95,9.18,9.47,9.76,10.09,10.46,10.88],feh:[-.99,-0.84,-0.70,-0.52,-0.31,-0.06,0.18],$
              feherr:[0.1,0.1,0.1,0.1,0.1,0.1,0.10,0.095]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   agerange = [0,13,0,2.,4.,6.,8.,13.]  
   for n=0,1 do begin 
   ;for n=0,n_Elements(agerange)-2 do begin
     inrange = where(sciall.ageform gt agerange(n) and sciall.ageform lt agerange(n+1),cinrange)
     if cinrange eq 0 then continue
     sci = sciall(inrange)
     goodness =sciall(inrange).goodsn
     goodfit = where(goodness eq 1,cgoodfit, complement=badfit,ncomplement=cbadfit)
     name = 'univage_'+strtrim(string(agerange(n),format='(I)'),2)+'_'+strtrim(string(agerange(n+1),format='(I)'),2) 
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      psname= name+'_alpha_feh_mass.eps'
      device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[9.,12]
         yrange=[-1.,0.3]
         plot,sci.logmstar,sci.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   
;         cgerrplot,sci(goodfit).logmstar,sci(goodfit).fehlower,sci(goodfit).fehupper,color='tomato',thick=0.5
;         if cbadfit gt 0 then oplot,sci(badfit).logmstar,sci(badfit).feh,psym=cgsymcat(16),color=fsc_Color('darkgray'),symsize=0.7
         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'tomato'

           oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).feh],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         oplot,!x.crange,(!x.crange-10.)*sdssmzrpar[1]+sdssmzrpar[0],linestyle=2,thick=2,color=fsc_color('navyblue')
         oplot,!x.crange,(!x.crange-10.)*cl0024mzrpar[1]+cl0024mzrpar[0],linestyle=2,thick=2,color=fsc_color('peru')
         oplot,!x.crange,(!x.crange-10.)*ms0451mzrpar[1]+ms0451mzrpar[0],linestyle=2,thick=2,color=fsc_color('indianred')
         ;xyouts,11.3,-0.8,strtrim(string(agerange(n),format='(I)'),2)+'-'+strtrim(string(agerange(n+1),format='(I)'),2)+'Gyr'
         al_legend,['z~0 SDSS','z~0.4 Cl0024','z~0.55 MS0451'],psym=cgsymcat(46),box=0,colors=['skyblue','goldenrod','tomato'],$
                    /bottom,/right,charsize=1,symsize=2
      device,/close
   
      psname=name+'_alpha_feh_alpha.eps'
      device, filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[-1.,0.3]
         yrange=[-0.5,0.9]
         plot,sci.feh,sci.alphafe,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
         ;draw axis
         axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='[Fe/H]'
         axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/Fe]'
         axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
         axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
   
   ;      cgerrplot,sci(goodfit).alphafe,sci(goodfit).alphafelower,sci(goodfit).alphafeupper,$
   ;                color='tomato',thick=0.5
   ;      cgerrplot,sci(goodfit).feh,sci(goodfit).fehlower,sci(goodfit).fehupper,$
   ;                color='tomato',thick=0.5,/horizontal
         if cbadfit gt 0 then oplot,sci(badfit).feh,sci(badfit).alphafe,psym=cgsymcat(16),color=fsc_Color('darkgray'),symsize=0.7
         goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'tomato'

           oplot,[sci(goodfit[i]).feh],[sci(goodfit[i]).alphafe],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         xyouts,-0.9,-0.4,strtrim(string(agerange(n),format='(I)'),2)+'-'+strtrim(string(agerange(n+1),format='(I)'),2)+'Gyr'
      device,/close
   
      psname = name+'_alpha_alpha_mass.eps'
      device,filename = psname,xsize = 15,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         xrange=[9.5,12]
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
   
 ;        if cbadfit gt 0 then oplot,sci(badfit).logmstar,sci(badfit).ah,psym=cgsymcat(16),color=fsc_Color('darkgray'),symsize=0.7
         for i=0,cgoodfit-1 do begin
           if sci(goodfit[i]).zspec lt 0.3 then symcol = 'skyblue'
           if sci(goodfit[i]).zspec gt 0.3 and sci(goodfit[i]).zspec lt 0.5 then symcol = 'goldenrod'
           if sci(goodfit[i]).zspec gt 0.5 and sci(goodfit[i]).zspec lt 0.6 then symcol = 'tomato'
           oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).ah],psym=cgsymcat(46),color=fsc_Color(symcol),symsize=goodsymsize[i]
         endfor
         oplot,!x.crange,(!x.crange-10.)*sdssmarpar[1]+sdssmarpar[0],linestyle=2,thick=2,color=fsc_color('navyblue')
         oplot,!x.crange,(!x.crange-10.)*cl0024marpar[1]+cl0024marpar[0],linestyle=2,thick=2,color=fsc_color('peru')
         oplot,!x.crange,(!x.crange-10.)*ms0451marpar[1]+ms0451marpar[0],linestyle=2,thick=2,color=fsc_color('indianred')
         al_legend,['z~0','z~0.39','z~0.55'],psym=cgsymcat(46),box=0,symsize=2,$
                colors=['skyblue','goldenrod','tomato'],charsize=1.2,position=[11.1,-0.4]
         xyouts,11.5,-0.4,'Observed Redshift',alignment=0.5,charsize=1.2
      device,/close
   endfor
stop
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;plot deviation as a function of formation age
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='FEH_deviation_fromz0line.eps'
   sdss_fixslope_pars= sdssmzrpar
   fehideal = sdss_fixslope_pars[0]+(sciall.logmstar-10.)*sdss_fixslope_pars[1]
   x=sciall.ageform
   y=sciall.feh-fehideal
   dx = 0.5*(sciall.ageupper-sciall.agelower)
   dy = 0.5*(sciall.fehupper-sciall.fehlower)
   dxarr = -1.*proball.dagearr ;since x axis is ageuni-agegal
   dyarr = proball.dfeharr
   probdxdyarr = transpose(proball.probfehage,[1,0,2]) ;age feh nobjs
   redolinfit = 0
   if redolinfit eq 1 then begin
     feh_dev_par = linfitprobgrid_general(x,y,dxarr,dyarr,probdxdyarr,sigma_a_b,$
                   outfitfile='all_prob_deltafeh_age.fits',outepsfile='all_plot_deltafeh_age.eps',$
                   fitrange=[0,11.5],/plot,initial_guess=[-0.4,0.06],stepsize=[0.0005,0.0005],$
                   badid=[34,248],nstep=3000)

      save, feh_dev_par,sigma_a_b,filename='linfitprobgrid_param_all.sav'
   endif else begin
      restore,'linfitprobgrid_param_all.sav'
      feh_dev_par = sigma_a_b[1,*]
   endelse
   redolinfitmc = 0
   if redolinfitmc eq 1 then begin
      mcpar = linfitprobgrid_mc(x,y,dxarr,dyarr,probdxdyarr,sigmamc)
      save, mcpar,sigmamc,filename='linfitmc_param_all.sav'
   endif else begin
      restore,'linfitmc_param_all.sav'
   endelse

   print,'feh deviation (intercept, slope)'
   print,'probgrid fit',sigma_a_b
   print,'mc fit',mcpar,sigmamc
   ;stupid linfit
   good = where(x gt 0.)
   simplepar = linfit(x(good),y(good),sigma=simpleparerr)
   print,'simple fit',simplepar,simpleparerr      

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

   set_plot,'ps'
   device, filename = psname,xsize = 15,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;
;      loadct,54
;      good = where(finite(img),cgood,complement=bad)
;      img(bad) = min(img(good))
;      loadct,55
;      cgimage,img,position=plotposition
;      cgimage,img,transparent=50,alphabackgroundimage=bgimage,alphafgposition=plotposition
      plot,sciall.ageform,sciall.feh-fehideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle=delta+'[Fe/H]',xrange=xfitrange,xstyle=9,yrange=yfitrange,/nodata,position=plotposition,/noerase

      axis,xaxis=1,xticks=4,xtickv=[1.513,3.223, 5.747,8.422 ,12.161],xtickn=['4','2','1','0.5','0.1'],xtitle='z'

      sdss = where(sciall.zspec lt 0.3)
      cgplot,sciall(sdss).ageform,sciall(sdss).feh-fehideal(sdss),psym=16,/overplot,symsize=0.5,color='skyblue'
      cl0024 = where(sciall.zspec gt 0.3 and sciall.zspec lt 0.45)
      cgplot,sciall(cl0024).ageform,sciall(cl0024).feh-fehideal(cl0024),psym=16,/overplot,symsize=0.5,color='goldenrod'
      ms0451 = where(sciall.zspec gt 0.45)
      cgplot,sciall(ms0451).ageform,sciall(ms0451).feh-fehideal(ms0451),psym=16,/overplot,symsize=0.5,color='tomato'
      bad = where(sciall.goodsn eq 0)
      cgplot,sciall(bad).ageform,sciall(ms0451).feh-fehideal(ms0451),psym=16,/overplot,symsize=0.5,color='gray'

;      oplot,[0,14],feh_dev_par(0)+feh_dev_par(1)*[0,14],thick=1,linestyle=2
      oplot,[0,14],mcpar(0)+mcpar(1)*[0,14],thick=1,linestyle=2;,color=fsc_color('red') 
      ;oplot,[0,14],feh_dev_par_sdss(0)+feh_dev_par_sdss(1)*[0,14],thick=2,linestyle=2,color=fsc_color('blu5')
      ;oplot,[0,14],feh_dev_par_cl(0)+feh_dev_par_cl(1)*[0,14],thick=2,linestyle=2,color=fsc_color('org6')

  ;    xyouts,0.5,0.32,'(c) real observations',charsize=1.2
      xyouts,1,-0.38,'older galaxies',charsize=0.8
      xyouts,9.4,-0.38,'younger galaxies',charsize=0.8
      arrow,0.9,-0.37,0.3,-0.37,/data
      arrow,12.1,-0.37,12.7,-0.37,/data
   device,/close

   ahideal = sdssmarpar[0]+(sciall.logmstar-10.)*sdssmarpar[1]
   x=sciall.ageform
   y=sciall.ah-ahideal
   ;stupid linfit
   good = where(y gt -0.5)
   simplepar = linfit(x(good),y(good),sigma=simpleparerr)
   ;mc fit
   redolinfitmc = 1
   if redolinfitmc eq 1 then begin
      good = where(y gt -0.5)
      y1 = sciall(good).feh
      y2 = sciall(good).alphafe
      dxarr = -1.*proball(good).dagearr ;since x axis is ageuni-agegal
      dy1arr = proball(good).dfeharr
      dy2arr = proball(good).dalphaarr
      probcube = transpose(proball(good).cubeprob,[1,0,2,3])
      mcpar = linfitprobgrid_mc_mgh(x(good),y1,y2,dxarr,dy1arr,dy2arr,probcube,sigmamc)
      save, mcpar,sigmamc,filename='linfitmc_mghparam_all.sav'
   endif else begin
      restore,'linfitmc_mghparam_all.sav'
   endelse

   print,'delta Mg'
   print,'simple fit',simplepar,simpleparerr
   print,'mc fit',mcpar,sigmamc
   set_plot,'ps'
   psname = 'mgh_deviation_fromz0line.eps'
   device, filename = psname,xsize = 15,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;
      yfitrange=[-0.5,0.6]
      plot,sciall.ageform,sciall.ah-ahideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle=delta+'[Mg/H]',xrange=xfitrange,xstyle=9,yrange=yfitrange,/nodata,position=plotposition,/noerase,ystyle=1

      axis,xaxis=1,xticks=4,xtickv=[1.513,3.223, 5.747,8.422 ,12.161],xtickn=['4','2','1','0.5','0.1'],xtitle='z'

      sdss = where(sciall.zspec lt 0.3)
      cgplot,sciall(sdss).ageform,sciall(sdss).ah-ahideal(sdss),psym=16,/overplot,symsize=0.5,color='skyblue'
      cl0024 = where(sciall.zspec gt 0.3 and sciall.zspec lt 0.45)
      cgplot,sciall(cl0024).ageform,sciall(cl0024).ah-ahideal(cl0024),psym=16,/overplot,symsize=0.5,color='goldenrod'
      ms0451 = where(sciall.zspec gt 0.45)
      cgplot,sciall(ms0451).ageform,sciall(ms0451).ah-ahideal(ms0451),psym=16,/overplot,symsize=0.5,color='tomato'
      bad = where(sciall.goodsn eq 0)
      cgplot,sciall(bad).ageform,sciall(ms0451).ah-ahideal(ms0451),psym=16,/overplot,symsize=0.5,color='gray'

      oplot,[0,14],simplepar(0)+simplepar(1)*[0,14],thick=1,linestyle=2
      ;oplot,[0,14],ah_dev_par_sdss(0)+ah_dev_par_sdss(1)*[0,14],thick=2,linestyle=2,color=fsc_color('blu5')
      ;oplot,[0,14],ah_dev_par_cl(0)+ah_dev_par_cl(1)*[0,14],thick=2,linestyle=2,color=fsc_color('org6')

  ;    xyouts,0.5,0.32,'(c) real observations',charsize=1.2
      xyouts,1,-0.48,'older galaxies',charsize=0.8
      xyouts,9.4,-0.48,'younger galaxies',charsize=0.8
      arrow,0.9,-0.37,0.3,-0.37,/data
      arrow,12.1,-0.37,12.7,-0.37,/data
   device,/close


end
