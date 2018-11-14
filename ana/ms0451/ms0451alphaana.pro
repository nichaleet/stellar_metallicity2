pro ms0451alphaana
   ;GETTING SCI DATA
   sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
   stack = mrdfits('/scr2/nichal/workspace4/sps_fit/data/stacked_ms0451/sps_fit03.fits.gz',1)
   ;sample selection
   good = where(sci.snfit gt 2 and sci.logmstar gt 6. and sci.oiiew gt -5.,cgood)
   sci = sci(good)
   ageform = (galage(sci.zfit,1000.)/1.e9-sci.age)>0. ;age of universe when it was formed

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;REFERENCE VALUES
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   restore,'leetho18a_avevals.sav'
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
  ;MILES Library abundance (from Conroy17)
  miles = {feh:[-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],$
           ofe:[0.6,0.5,0.5,0.4,0.3,0.2,0.2,0.1,0.0,0.0],$
           mgfe:[0.4,0.4,0.4,0.4,0.34,0.22,0.14,0.11,0.05,0.04],$
           cafe:[0.32,0.30,0.28,0.26,0.26,0.17,0.12,0.06,0.00,0.00]}
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  ;fix the library abundance
  sci.alphafe = sci.alphafe+interpol(miles.mgfe,miles.feh,sci.feh)
  stack.alphafe = stack.alphafe+interpol(miles.mgfe,miles.feh,stack.feh) 
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.5
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'

   goodfit = where(sci.snfit gt 6.5 and sci.logmstar gt 6.,cgoodfit,complement=badfit,ncomplement=cbadfit)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   psname='ms0451alpha_feh_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[9.5,12]
      yrange=[-1.,0.3]
      plot,sci.logmstar,sci.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
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

      cgerrplot,sci(goodfit).logmstar,sci(goodfit).fehlower,sci(goodfit).fehupper,color='tomato',thick=0.5
      oplot,sci(badfit).logmstar,sci(badfit).feh,psym=cgsymcat(16),color=fsc_Color('cadetblue'),symsize=0.7
      goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
      for i=0,cgoodfit-1 do begin
        oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).feh],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=goodsymsize[i]
        oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).feh],psym=cgsymcat(45),color=fsc_color('maroon'),symsize=goodsymsize[i]
      endfor
      ;draw stacked results
      oplot,stack.logmstar,stack.feh,psym=cgsymcat(16),color=fsc_color('indianred'),symsize=2
   device,/close

   psname='ms0451alpha_feh_alpha.eps'
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
      oplot,sci(badfit).feh,sci(badfit).alphafe,psym=cgsymcat(16),color=fsc_Color('cadetblue'),symsize=0.7
      goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
      for i=0,cgoodfit-1 do begin
        oplot,[sci(goodfit[i]).feh],[sci(goodfit[i]).alphafe],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=goodsymsize[i]
        oplot,[sci(goodfit[i]).feh],[sci(goodfit[i]).alphafe],psym=cgsymcat(45),color=fsc_color('maroon'),symsize=goodsymsize[i]
      endfor
      oplot,stack.feh,stack.alphafe,psym=cgsymcat(16),color=fsc_color('indianred'),symsize=2
   device,/close

   psname = 'ms0451alpha_alpha_mass.eps'
   ah = sci.alphafe+sci.feh
   ahupper = ah+sqrt((sci.alphafeupper-sci.alphafe)^2+(sci.fehupper-sci.feh)^2)
   ahlower = ah-sqrt((sci.alphafe-sci.alphafelower)^2+(sci.feh-sci.fehlower)^2)
   device,filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange=[9.5,12]
      yrange=[-1.,0.9]
      plot,sci.logmstar,ah,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M!L*!N/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Mg/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      goodsymsize = 1.5/(max(sci(goodfit).snfit)-min(sci(goodfit).snfit))*sci(goodfit).snfit+0.7
      for i=0,cgoodfit-1 do begin
        oplot,[sci(goodfit[i]).logmstar],[sci(goodfit[i]).feh],psym=cgsymcat(46),color=fsc_Color('lightgray'),symsize=goodsymsize[i]
      endfor

      oplot,sci(badfit).logmstar,ah(badfit),psym=cgsymcat(16),color=fsc_Color('cadetblue'),symsize=0.7
      for i=0,cgoodfit-1 do begin
        oplot,[sci(goodfit[i]).logmstar],[ah(goodfit[i])],psym=cgsymcat(46),color=fsc_Color('tomato'),symsize=goodsymsize[i]
        oplot,[sci(goodfit[i]).logmstar],[ah(goodfit[i])],psym=cgsymcat(45),color=fsc_color('maroon'),symsize=goodsymsize[i]
      endfor
      oplot,stack.logmstar,stack.feh+stack.alphafe,psym=cgsymcat(16),color=fsc_color('indianred'),symsize=2
   device,/close

   ;below is to print out RA and DECs of the high Fe/H vs lower Fe/H of the galaxies with M gt 11 Msun
   hifegroup = where(sci.logmstar gt 11. and sci.feh gt 0.)
   lofegroup = where(sci.logmstar gt 11. and sci.feh gt -0.4 and sci.feh lt 0.)
   write_ds9_regionfile, sci(hifegroup).ra, sci(hifegroup).dec, filename='ms0451_high_fe_group.reg',color='red'
   write_ds9_regionfile, sci(lofegroup).ra, sci(lofegroup).dec, filename='ms0451_low_fe_group.reg',color='green'

   himggroup = where(sci.logmstar gt 9.5 and ah gt sci.logmstar*1.4/2.5-5.82)
   lomggroup = where(sci.logmstar gt 9.5 and ah lt sci.logmstar*1.4/2.5-5.82)
   write_ds9_regionfile, sci(himggroup).ra, sci(himggroup).dec, filename='ms0451_high_mg_group.reg',color='red',symbol='diamond'
   write_ds9_regionfile, sci(lomggroup).ra, sci(lomggroup).dec, filename='ms0451_low_mg_group.reg',color='green',symbol='diamond'
   stop
end
