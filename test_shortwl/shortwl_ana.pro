pro shortwl_ana
   ;read data
   short = mrdfits('/scr2/nichal/workspace4/test_shortwl/datasps/test_shortwl_cl1604/sps_fit.fits.gz',1)

   ;get correct/input values
   nshort = n_elements(short)
   izmet = fltarr(nshort)
   iage = fltarr(nshort)
   isn = fltarr(nshort)
   for i=0,nshort-1 do begin
       exts = strsplit(short[i].objname,'_',/extract)
       mop = (strmid(exts[0],4,1) eq 'm') ? -1. : 1.
       izmet[i] = mop*float(strmid(exts[0],5,1))+0.01*float(strmid(exts[0],7,2))
       iage[i] = float(strmid(exts[1],0,1))
       isn[i] = float(strmid(exts(3),3))
   endfor

   ;set colors and symbols in the plots
   zmet_uniq = izmet[uniq(izmet,sort(izmet))]
   age_uniq  = iage[uniq(iage,sort(iage))]
   sn_uniq = isn[uniq(isn,sort(isn))]
   age_symbols = [14,16,46,24,22] ;for cgsymcat
   age_symbols = age_symbols[0:n_elements(age_uniq)-1]

   ;calculate differences
   dage = short.age-iage
   dageu = short.ageupper-iage
   dagel = short.agelower-iage
   dz  = short.feh-izmet 
   dzu = short.fehupper-izmet
   dzl = short.fehlower-izmet
  
   ;symbols
   Delta = '!9'+string("104B)+'!x'
   Delta = '!9'+string("104B)+'!x'
 
   ;plotting
   set_plot,'ps'
   !p.font=0
   psname='shortwl_deltaage_sn_fixage.eps'
   distinct_colors, N_COLORS = n_elements(zmet_uniq) ;0 is black, 255 = white and n colors in between
   device, filename = psname,xsize = 30,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,3,1]
      for i=0,n_elements(age_uniq)-1 do begin
         sel = where(iage eq age_uniq(i),csel)
         if csel eq 0 then stop,'oops1'
         plot,isn(sel),dage(sel),/nodata,xtitle='SN per A',ytitle= delta+'Age (Gyr)'
         for j=0,n_elements(zmet_uniq)-1 do begin
            sel = where(iage eq age_uniq(i) and izmet eq zmet_uniq(j), csel)
            if csel eq 0 then stop,'oop2'
            cgerrplot,isn(sel),dagel(sel),dageu(sel),color=fix(255/(j+1))
            oplot,isn(sel),dage(sel),psym=cgsymcat(age_symbols[i]),color=fix(255/(j+1))
            xyouts,0.6,0.9,'input age: '+strtrim(string(age_uniq(i),format='(I)'),2)+' Gyr'
         endfor 
         if i eq 1 then al_legend,zmet_unique,psym=cgsymcat(age_symbols[i]),colors=fix(255/(indgen(n_Elements(zmet_uniq))+1)),/top,/right
      endfor
   device,/close

   psname='shortwl_deltazmet_sn_fixage.eps'
   distinct_colors, N_COLORS = n_elements(zmet_uniq) ;0 is black, 255 = white and n colors in between
   device, filename = psname,xsize = 30,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,3,1]
      for i=0,n_elements(age_uniq)-1 do begin
         sel = where(iage eq age_uniq(i),csel)
         if csel eq 0 then stop,'oops1'
         plot,isn(sel),dz(sel),/nodata,xtitle='SN per A',ytitle=delta+'[z/H] (Gyr)'
         for j=0,n_elements(zmet_uniq)-1 do begin
            sel = where(iage eq age_uniq(i) and izmet eq zmet_uniq(j), csel)
            if csel eq 0 then stop,'oop2'
            cgerrplot,isn(sel),dzl(sel),dzu(sel),color=fix(255/(j+1))
            oplot,isn(sel),dz(sel),psym=cgsymcat(age_symbols[i]),color=fix(255/(j+1))
         endfor 
         if i eq 1 then al_legend,zmet_unique,psym=cgsymcat(age_symbols[i]),colors=fix(255/(indgen(n_Elements(zmet_uniq))+1)),/top,/right
      endfor
   device,/close

end
