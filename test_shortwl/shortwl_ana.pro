pro shortwl_ana
   ;read data
   short = mrdfits('/scr2/nichal/workspace4/test_shortwl/datasps/test_shortwl_cl1604/sps_fitall.fits.gz',1)
   full = mrdfits('/scr2/nichal/workspace4/test_shortwl/datasps/test_fullwl_cl1604/sps_fitall.fits.gz',1)

   for ff=0,1 do begin         
      if ff eq 0 then begin
         str = short
         name='short'
      endif
      if ff eq 1 then begin
         str = full
         name='full'
      endif
      ;get correct/input values
      nstr = n_elements(str)
      izmet = fltarr(nstr)
      iage = fltarr(nstr)
      isn = fltarr(nstr)
      for i=0,nstr-1 do begin
          exts = strsplit(str[i].objname,'_',/extract)
          mop = (strmid(exts[0],4,1) eq 'm') ? -1. : 1.
          izmet[i] = mop*(float(strmid(exts[0],5,1))+0.01*float(strmid(exts[0],7,2)))
          iage[i] = float(strmid(exts[1],0,1))
          isn[i] = float(strmid(exts(2),2))
      endfor
      ;set colors and symbols in the plots
      zmet_uniq = izmet[uniq(izmet,sort(izmet))]
      age_uniq  = iage[uniq(iage,sort(iage))]
      sn_uniq = isn[uniq(isn,sort(isn))]
      age_symbols = [14,16,15,46,24,22] ;for cgsymcat
      age_symbols = age_symbols[0:n_elements(age_uniq)-1]
      zmet_colors = ['blu6','teal','org4','red6','deeppink','gray','dodgerblue'] ;for fsc_colors
      zmet_colors = zmet_colors[0:n_elements(zmet_uniq)-1]
      ;calculate differences
      dage = alog10(str.age/iage)
      dageu = alog10(str.ageupper/iage)
      dagel = alog10(str.agelower/iage)
      dz  = str.feh-izmet 
      dzu = str.fehupper-izmet
      dzl = str.fehlower-izmet
     
      ;symbols
      Delta = '!9'+string("104B)+'!x'
      Delta = '!9'+string("104B)+'!x'
    
      ;plotting
      set_plot,'ps'
      !p.font=0
      !p.charsize=2
      psname='plots/'+name+'wl_deltaage_sn_fixage.eps'
      device, filename = psname,xsize = 30,ysize = 8, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         !p.multi = [0,3,1]
         for i=0,n_elements(age_uniq)-1 do begin
            sel = where(iage eq age_uniq(i),csel)
            if csel eq 0 then stop,'oops1'
            plot,isn(sel),dage(sel),/nodata,xtitle='SN per A',ytitle= delta+'log(Age/Gyr)',xrange=[min(isn)-5,max(isn)+5],$
                 title='input age: '+strtrim(string(age_uniq(i),format='(I)'),2)+' Gyr'
            oplot,!x.crange,[0,0],color=fsc_color('gray')
            for j=0,n_elements(zmet_uniq)-1 do begin
               sel = where(iage eq age_uniq(i) and izmet eq zmet_uniq(j), csel)
               if csel eq 0 then stop,'oop2'
               cgerrplot,isn(sel)+(j-2)*0.5,dagel(sel),dageu(sel),color=zmet_colors[j]
               oplot,isn(sel)+(j-2)*0.5,dage(sel),psym=cgsymcat(age_symbols[i]),color=fsc_color(zmet_colors[j])
            endfor 
            if i eq 0 then al_legend, string(zmet_uniq,format='(f5.2)'),psym=cgsymcat(age_symbols[i]),colors=zmet_colors,/bottom,/right,box=0,charsize=1
   
         endfor
      device,/close
   
      psname='plots/'+name+'wl_deltazmet_sn_fixage.eps'
      device, filename = psname,xsize = 30,ysize = 8, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         !p.multi = [0,3,1]
         for i=0,n_elements(age_uniq)-1 do begin
            sel = where(iage eq age_uniq(i),csel)
            if csel eq 0 then stop,'oops1'
            plot,isn(sel),dz(sel),/nodata,xtitle='SN per A',ytitle= delta+'[Z/H]',$
                xrange=[min(isn)-5,max(isn)+5],title='input age: '+strtrim(string(age_uniq(i),format='(I)'),2)+' Gyr'
            oplot,!x.crange,[0,0],color=fsc_color('gray')
            for j=0,n_elements(zmet_uniq)-1 do begin
               sel = where(iage eq age_uniq(i) and izmet eq zmet_uniq(j), csel)
               if csel eq 0 then stop,'oop2'
               cgerrplot,isn(sel)+(j-2)*0.5,dzl(sel),dzu(sel),color=zmet_colors[j]
               oplot,isn(sel)+(j-2)*0.5,dz(sel),psym=cgsymcat(age_symbols[i]),color=fsc_color(zmet_colors[j])
            endfor
            if i eq 0 then al_legend, string(zmet_uniq,format='(f5.2)'),psym=cgsymcat(age_symbols[i]),colors=zmet_colors,/bottom,/right,box=0,charsize=1
         endfor
      device,/close
   endfor
end
