pro compare_alf_ms0451
;read alf
  alf_files = file_search('/scr2/nichal/alfresults/default_config/*.sum',count=nfile)
;read sps
  sciall = mrdfits('/scr2/nichal/workspace4/ana/all/all_sci_ms05.fits',1)

  alf_results = replicate({objname:'',feh:0.,feh84:0.,feh16:0.,feh98:0.,feh03:0.,$
                           logage:0.,logage84:0.,logage16:0.,logage98:0.,logage03:0.,$
                           mgfe:0.,mgfe84:0.,mgfe16:0.,mgfe98:0.,mgfe03:0.},nfile)
  nl_results = replicate({objname:'',feh:0.,feh84:0.,feh16:0.,$
                          logage:0.,logage84:0.,logage16:0.,$
                          mgfe:0.,mgfe84:0.,mgfe16:0.},nfile)  
  for i=0,nfile-1 do begin
     filenow = file_basename(alf_files[i])
     objname = (strsplit(filenow,'_',/extract))[0]
     alf = read_alf(filenow)
     print, objname
     loc = where(strmatch(sciall.objname,objname+'*') eq 1,cmatch)
     if cmatch ne 1 then stop,'something is wrong'
     nl = sciall[loc]
     alf_results[i].objname = objname
     nl_results[i].objname = objname

     alf_results[i].feh = alf[1].feh
     alf_results[i].feh84 = alf[6].feh
     alf_results[i].feh16= alf[4].feh
     alf_results[i].feh98 = alf[7].feh
     alf_results[i].feh03= alf[3].feh

     alf_results[i].logage = alf[1].logage
     alf_results[i].logage84 = alf[6].logage
     alf_results[i].logage16= alf[4].logage
     alf_results[i].logage98 = alf[7].logage
     alf_results[i].logage03= alf[3].logage

     alf_results[i].mgfe = alf[0].mgfe
     alf_results[i].mgfe84 = alf[1].mgfe+alf[2].mgfe;alf[6].mgfe
     alf_results[i].mgfe16= alf[1].mgfe-alf[2].mgfe;alf[4].mgfe
     alf_results[i].mgfe98 = alf[7].mgfe
     alf_results[i].mgfe03= alf[3].mgfe
     print, alf.mgfe
     nl_results[i].feh = nl.feh
     nl_results[i].feh84 = nl.fehupper
     nl_results[i].feh16= nl.fehlower
     nl_results[i].logage = alog10(nl.age)
     nl_results[i].logage84 = alog10(nl.ageupper)
     nl_results[i].logage16= alog10(nl.agelower)
     nl_results[i].mgfe = nl.alphafe
     nl_results[i].mgfe84 = nl.alphafeupper
     nl_results[i].mgfe16= nl.alphafelower
 
  endfor
   set_plot,'ps'
   !p.multi = [0,2,2]
   !p.font = 0
   !p.charsize = 1
   sunsym = sunsymbol()
   psname= 'compare_alf.eps'
   sym = [14,15,16,17,22,24,18]
   col = ['saddlebrown', 'firebrick','darkorchid','ryb8','goldenrod','pink','orangered']
   device, filename = psname,xsize = 20,ysize = 20, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,alf_results.feh,nl_results.feh,psym=1,/nodata,xtitle='[Fe/H] ALF',ytitle='[Fe/H] MPFIT',xrange=[-0.3,0.3],yrange=[-0.3,0.3],xstyle=1,ystyle=1
      cgerrplot,alf_results.feh,nl_results.feh16,nl_results.feh84
      cgerrplot,nl_results.feh,alf_results.feh16,alf_results.feh84,/horizontal
      ;cgerrplot,nl_results.feh,alf_results.feh03,alf_results.feh98,/horizontal,color=fsc_color('gray')
      for i=0,nfile-1 do oplot,[alf_results[i].feh],[nl_results[i].feh],psym=cgsymcat(sym[i]),color=fsc_color(col[i])
      oplot,!x.crange,!x.crange,color=fsc_color('red')
     
      plot,alf_results.logage,nl_results.logage,psym=1,/nodata,xtitle='log Age ALF',ytitle='log Age MPFIT',xrange=[0,1],yrange=[0,1]
      cgerrplot,alf_results.logage,nl_results.logage16,nl_results.logage84
      cgerrplot,nl_results.logage,alf_results.logage16,alf_results.logage84,/horizontal
      for i=0,nfile-1 do oplot,[alf_results[i].logage],[nl_results[i].logage],psym=cgsymcat(sym[i]),color=fsc_color(col[i])
      oplot,!x.crange,!x.crange,color=fsc_color('red')

      plot,alf_results.mgfe,nl_results.mgfe,psym=1,/nodata,xtitle='[Mg/Fe] ALF',ytitle='[Mg/Fe] MPFIT',xrange=[-0.1,0.5],yrange=[-0.1,0.5],xstyle=1,ystyle=1
      cgerrplot,alf_results.mgfe,nl_results.mgfe16,nl_results.mgfe84
      cgerrplot,nl_results.mgfe,alf_results.mgfe16,alf_results.mgfe84,/horizontal
      for i=0,nfile-1 do oplot,[alf_results[i].mgfe],[nl_results[i].mgfe],psym=cgsymcat(sym[i]),color=fsc_color(col[i])
      oplot,!x.crange,!x.crange,color=fsc_color('red')
   device,/close
stop
end
