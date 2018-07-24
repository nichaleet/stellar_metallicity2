pro cl1604alphaana
   ;GETTING SCI DATA
   scims = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit01.fits.gz',1)
   scicl = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_cl1604/sps_fit01.fits.gz',1)
   stack = mrdfits('/scr2/nichal/workspace4/sps_fit/data/stacked_cl1604/sps_fit01.fits.gz',1)
   ;sample selection
   good = where(scims.oiiew gt -5.,cgood)
   scims = scims(good)

   good = where(scicl.oiiew gt -5. and scicl.good eq 1,cgood)
   scicl = scicl(good)   

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
   sunsym = sunsymbol()
   alpha = '!9'+string("141B)+'!x'
   Delta = '!9'+string("104B)+'!x'

   ngals = n_elements(scicl)
   flagmosfire = bytarr(ngals) 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   sn = scicl.snfit
   dispersion = fltarr(ngals)
   for i=0,ngals-1 do begin
      lamnow = scicl[i].lambda
      goodlam = where(lamnow ne 0. and scicl[i].contdivivar ne 0.)
      lamnow = lamnow(goodlam)
      lamdiff = -1.*ts_diff(lamnow(goodlam),1)
      dispersion[i] = median(lamdiff)
      if max(lamnow)/(1.+scicl[i].zspec) gt 5300. then flagmosfire[i]=1
   endfor
   wmosfire = where(flagmosfire eq 1,cmosfire,complement=wnomosfire,ncomplement=cnomosfire)
    
   wnofit = where(scicl.snfit eq 0,cnofit)
   if cnofit gt 0 then sn(wnofit) = scicl(wnofit).sn
   ;This sn is per pixel. With 600 grating, the resolution is 0.65 A per pixel.
   ;To convert to per angstrom divide by sqrt(0.65) = x1.24
   sn = sn/sqrt(dispersion)

   goodmsfit = where(scims.goodfit eq 1 and scims.snfit gt 10, cgoodmsfit, complement=badmsfit,ncomplement=cbadmsfit)
   goodsymsize = 1.5/(max(scims(goodmsfit).snfit)-min(scims(goodmsfit).snfit))*scims(goodmsfit).snfit+0.7

   goodclfit = where(scicl.good eq 1, cgoodclfit)
   goodwmosfire = where(scicl.good eq 1 and flagmosfire eq 1,cgoodwmosfire)
   goodwnomosfire = where(scicl.good eq 1 and flagmosfire eq 0, cgoodwnomosfire)
   goodwmosfiresymsize = 1.5/(max(scims(goodmsfit).snfit)-min(scims(goodmsfit).snfit))*scicl(goodwmosfire).snfit+0.7
   goodwnomosfiresymsize = 1.5/(max(scims(goodmsfit).snfit)-min(scims(goodmsfit).snfit))*scicl(goodwnomosfire).snfit+0.7

   psname='cl1604alpha_feh_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[9.5,12]
      yrange=[-1.,0.3]
      plot,scims.logmstar,scims.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      ;colorarr=[10,50,203,230,254]
      colorarr= fsc_color(['royalblue','darkorchid','deeppink','maroon','red8'])
      x=[bndry_mass_cl0024,reverse(bndry_mass_cl0024)]
      y=[hifeh_cl0024,reverse(lofeh_cl0024)]
      polyfill,x,y,color=fsc_color('grn2')

      ylefthi = interpol(hifeh_sdss,bndry_mass_sdss,[xrange[0]])
      yleftlo = interpol(lofeh_sdss,bndry_mass_sdss,[xrange[0]])
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      goodx = where(x gt xrange[0])
      x= x(goodx)
      y= y(goodx)
      x= [xrange[0],x,xrange[0]]
      y= [ylefthi,y,yleftlo]
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135

      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
     
     ;oplot ms stuff
      cgerrplot,scims(goodmsfit).logmstar,scims(goodmsfit).fehlower,scims(goodmsfit).fehupper,color='tomato',thick=0.5
      ;oplot,scims(badmsfit).logmstar,scims(badmsfit).feh,psym=cgsymcat(16),color=fsc_Color('tomato'),symsize=0.7
      for i=0,cgoodmsfit-1 do begin
        oplot,[scims(goodmsfit[i]).logmstar],[scims(goodmsfit[i]).feh],psym=cgsymcat(46),color=fsc_Color('tomato'),$
              symsize=goodsymsize[i]
        oplot,[scims(goodmsfit[i]).logmstar],[scims(goodmsfit[i]).feh],psym=cgsymcat(45),color=fsc_color('brown'),$
              symsize=goodsymsize[i]
     endfor

     ;oplot cl stuff
     for i=0,cgoodwnomosfire-1 do oplot,[scicl(goodwnomosfire[i]).logmstar_sed_bl],[scicl(goodwnomosfire[i]).feh],psym=cgsymcat(46),color=fsc_Color('thistle'),symsize=goodwnomosfiresymsize[i]
     for i=0,cgoodwmosfire-1 do oplot,[scicl(goodwmosfire[i]).logmstar_sed_bl],[scicl(goodwmosfire[i]).feh],psym=cgsymcat(46),color=fsc_Color('darkorchid'),symsize=goodwmosfiresymsize[i]
     oplot,stack.logmstar_sed_bl,stack.feh,psym=cgsymcat(16),color=fsc_color('darkorchid'),symsize=2
   device,/close

  psname='cl1604alpha_feh_alpha.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[-1.,0.3]
      yrange=[-0.5,0.5]
      plot,scicl.feh,scicl.alphafe,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='[Zmet/H]'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='['+alpha+'/Zmet]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

;      cgerrplot,scicl(goodfit).alphafe,scicl(goodfit).alphafelower,scicl(goodfit).alphafeupper,$
;                color='tomato',thick=0.5
;      cgerrplot,scicl(goodfit).feh,scicl(goodfit).fehlower,scicl(goodfit).fehupper,$
;                color='tomato',thick=0.5,/horizontal
      for i=0,cgoodwnomosfire-1 do $
        oplot,[scicl(goodwnomosfire[i]).feh],[scicl(goodwnomosfire[i]).alphafe],psym=cgsymcat(46),color=fsc_Color('thistle'),symsize=goodwnomosfiresymsize[i]
      for i=0,cgoodwmosfire-1 do $
        oplot,[scicl(goodwmosfire[i]).feh],[scicl(goodwmosfire[i]).alphafe],psym=cgsymcat(46),color=fsc_color('darkorchid'),symsize=goodwmosfiresymsize[i]
      oplot,stack.feh,stack.alphafe,psym=cgsymcat(16),color=fsc_color('darkorchid'),symsize=2
   device,/close



   stop
end
