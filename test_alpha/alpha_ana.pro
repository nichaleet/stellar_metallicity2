pro alpha_ana
   ;read data
   nwalkers = 5

   for ff=0,3 do begin         
      if ff eq 0 then begin
         str =  mrdfits('/scr2/nichal/workspace4/test_alpha/datasps/mockdata_alpha/sps_fit01.fits.gz',1)
         name='01alphafit_5para'
      endif
      if ff eq 1 then begin
         str = mrdfits('/scr2/nichal/workspace4/test_alpha/datasps/mockdata_alpha/sps_fit02.fits.gz',1)
         str.alphafe = 0
         name='02feonlyfit_4para'
      endif
      if ff eq 2 then begin
         str = mrdfits('/scr2/nichal/workspace4/test_alpha/datasps/mockdata_alpha/sps_fit03.fits.gz',1)
         name = '03fitmgtoalpha'
      endif
      if ff eq 3 then begin
         str = mrdfits('/scr2/nichal/workspace4/test_alpha/datasps/mockdata_alpha/sps_fit04.fits.gz',1)
         str.alphafe = 0
         str.alphafeerr = 0
         name='04feonlyfit_4para_maskemlines'
      endif
      ;get correct/input values
      nstr = n_elements(str)
      izmet = fltarr(nstr) ;i is for initial
      iage = fltarr(nstr)
      isn = fltarr(nstr)
      ialpha = fltarr(nstr)
      for i=0,nstr-1 do begin
          exts = strsplit(str[i].objname,'_',/extract)
          mop = (strmid(exts[0],4,1) eq 'm') ? -1. : 1.
          izmet[i] = mop*(float(strmid(exts[0],5,1))+0.01*float(strmid(exts[0],7,2)))
          iage[i] = float(strmid(exts[1],0,1))
          mop = (strmid(exts[2],5,1) eq 'm') ? -1. : 1
          ialpha[i] = mop*(float(strmid(exts[2],6,1))+0.01*float(strmid(exts[2],8,2)))
          isn[i] = float(strmid(exts(3),2))
      endfor
      ;set colors and symbols in the plots
      zmet_uniq = izmet[uniq(izmet,sort(izmet))] ;should be 3 values
      age_uniq  = iage[uniq(iage,sort(iage))]    ;should be 1 value
      sn_uniq = isn[uniq(isn,sort(isn))]         ;should be 5 values
      alpha_uniq = ialpha[uniq(ialpha,sort(ialpha))] ;should be 5 values
      zmet_symbols = [14,16,15,46,24,22] ;for cgsymcat
      zmet_symbols = zmet_symbols[0:n_elements(zmet_uniq)-1]
      alpha_colors = ['blu6','teal','org4','red6','deeppink','gray','dodgerblue'] ;for fsc_colors
      alpha_colors = alpha_colors[0:n_elements(alpha_uniq)-1]
      ;calculate differences
      dage = alog10(str.age/iage)
      dageu = alog10(str.ageupper/iage)
      dagel = alog10(str.agelower/iage)
      dz  = str.feh-izmet 
      dzu = str.fehupper-izmet
      dzl = str.fehlower-izmet
      dalpha = str.alphafe-ialpha
      dalphau = str.alphafeupper-ialpha
      dalphal = str.alphafelower-ialpha

      ;symbols
      Delta = '!9'+string("104B)+'!x'
      Delta = '!9'+string("104B)+'!x'
    
      ;plotting
      if n_elements(age_uniq) gt 1 then stop,'!!!!!WARNING: Plots might not be meaningful. There are more than 1 age values!!!!!!!'
      set_plot,'ps'
      !p.font=0
      !p.charsize=2
      psname='plots/'+name+'deltaage_sn_fixzmet.eps'
      device, filename = psname,xsize = 30,ysize = 8, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         !p.multi = [0,n_elements(zmet_uniq),1]
         for i=0,n_elements(zmet_uniq)-1 do begin
            sel = where(izmet eq zmet_uniq(i),csel)
            if csel eq 0 then stop,'oops1'
            plot,isn(sel),dage(sel),/nodata,xtitle='SN per A',ytitle= delta+'log(Age/Gyr)',xrange=[min(isn)-5,max(isn)+5],$
                 title='input zmet: '+strtrim(string(zmet_uniq(i),format='(f5.2)'),2)
            oplot,!x.crange,[0,0],color=fsc_color('gray')
            for j=0,n_elements(alpha_uniq)-1 do begin
               sel = where(izmet eq zmet_uniq(i) and ialpha eq alpha_uniq(j), csel)
               if csel eq 0 then stop,'oop2'
               cgerrplot,isn(sel)+(j-2)*0.5,dagel(sel),dageu(sel),color=alpha_colors[j]
               oplot,isn(sel)+(j-2)*0.5,dage(sel),psym=cgsymcat(zmet_symbols[i]),color=fsc_color(alpha_colors[j])
            endfor 
            if i eq 0 then al_legend, string(alpha_uniq,format='(f5.2)'),psym=cgsymcat(zmet_symbols[i]),colors=alpha_colors,/bottom,/right,box=0,charsize=1
   
         endfor
      device,/close
   
      psname='plots/'+name+'deltazmet_sn_fixzmet.eps'
      device, filename = psname,xsize = 30,ysize = 8, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         !p.multi = [0,3,1]
         for i=0,n_elements(zmet_uniq)-1 do begin
            sel = where(izmet eq zmet_uniq(i),csel)
            if csel eq 0 then stop,'oops1'
            plot,isn(sel),dz(sel),/nodata,xtitle='SN per A',ytitle= delta+'[Z/H]',$
                xrange=[min(isn)-5,max(isn)+5],title='input zmet: '+strtrim(string(zmet_uniq(i),format='(f5.2)'),2)
            oplot,!x.crange,[0,0],color=fsc_color('gray')
            for j=0,n_elements(alpha_uniq)-1 do begin
               sel = where(izmet eq zmet_uniq(i) and ialpha eq alpha_uniq(j), csel)
               if csel eq 0 then stop,'oop2'
               cgerrplot,isn(sel)+(j-2)*0.5,dzl(sel),dzu(sel),color=alpha_colors[j]
               oplot,isn(sel)+(j-2)*0.5,dz(sel),psym=cgsymcat(zmet_symbols[i]),color=fsc_color(alpha_colors[j])
            endfor
            if i eq 0 then al_legend, string(alpha_uniq,format='(f5.2)'),psym=cgsymcat(zmet_symbols[i]),colors=alpha_colors,/bottom,/right,box=0,charsize=1
         endfor
      device,/close

      psname='plots/'+name+'deltaalpha_sn_fixzmet.eps'
      device, filename = psname,xsize = 30,ysize = 8, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         !p.multi = [0,3,1]
         for i=0,n_elements(zmet_uniq)-1 do begin
            sel = where(izmet eq zmet_uniq(i),csel)
            if csel eq 0 then stop,'oops1'
            plot,isn(sel),dalpha(sel),/nodata,xtitle='SN per A',ytitle= delta+'[alpha/H]',$
                xrange=[min(isn)-5,max(isn)+5],title='input zmet: '+strtrim(string(zmet_uniq(i),format='(f5.2)'),2)
            oplot,!x.crange,[0,0],color=fsc_color('gray')
            for j=0,n_elements(alpha_uniq)-1 do begin
               sel = where(izmet eq zmet_uniq(i) and ialpha eq alpha_uniq(j), csel)
               if csel eq 0 then stop,'oop2'
               cgerrplot,isn(sel)+(j-2)*0.5,dalphal(sel),dalphau(sel),color=alpha_colors[j]
               oplot,isn(sel)+(j-2)*0.5,dalpha(sel),psym=cgsymcat(zmet_symbols[i]),color=fsc_color(alpha_colors[j])
            endfor
            if i eq 0 then al_legend, string(alpha_uniq,format='(f5.2)'),psym=cgsymcat(zmet_symbols[i]),colors=alpha_colors,/bottom,/right,box=0,charsize=1
         endfor
      device,/close
   endfor
end
