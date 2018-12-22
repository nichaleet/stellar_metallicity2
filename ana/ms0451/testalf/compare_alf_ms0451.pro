pro compare_alf_ms0451
;read alf
  diralf = '/scr2/nichal/alfresults/default_config/'
  alf_files = file_search(diralf+'*.sum',count=nfile)
;read sps
  mgn = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
  mgo = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit06.fits.gz',1)
; fix mgh with library enhancement
; MILES Library abundance (from Conroy17)
  miles = {feh:[-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],$
           ofe:[0.6,0.5,0.5,0.4,0.3,0.2,0.2,0.1,0.0,0.0],$
           mgfe:[0.4,0.4,0.4,0.4,0.34,0.22,0.14,0.11,0.05,0.04],$
           cafe:[0.32,0.30,0.28,0.26,0.26,0.17,0.12,0.06,0.00,0.00]}
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  ;fix the library abundance
   tofix = interpol(miles.mgfe,miles.feh,mgn.feh)
   mgn.alphafe = mgn.alphafe+tofix
   mgn.alphafeupper = mgn.alphafeupper+tofix
   mgn.alphafelower = mgn.alphafelower+tofix
   tofix = interpol(miles.mgfe,miles.feh,mgo.feh)
   mgo.alphafe = mgo.alphafe+tofix
   mgo.alphafeupper = mgo.alphafeupper+tofix
   mgo.alphafelower = mgo.alphafelower+tofix
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  alf_results = replicate({objname:'',feh:0.,feh84:0.,feh16:0.,feh98:0.,feh03:0.,$
                           logage:0.,logage84:0.,logage16:0.,logage98:0.,logage03:0.,$
                           mgfe:0.,mgfe84:0.,mgfe16:0.,mgfe98:0.,mgfe03:0.},nfile)
  nl_mgn = replicate({objname:'',feh:0.,feh84:0.,feh16:0.,$
                          logage:0.,logage84:0.,logage16:0.,$
                          mgfe:0.,mgfe84:0.,mgfe16:0.},nfile)  
  nl_mgo = nl_mgn

  for i=0,nfile-1 do begin
     filenow = file_basename(alf_files[i])
     objname = (strsplit(filenow,'_',/extract))[0]
     
     alf = read_alf(filenow,dir=diralf)
     print, objname
     loc = where(strmatch(mgn.objname,objname+'*') eq 1,cmatch)
     if cmatch ne 1 then stop,'something is wrong'
     if mgn(loc).ra-mgo(loc).ra ne 0 then stop, 'the two sps_fits filesdo not match. stopped'

     nl1 = mgn[loc]
     nl2 = mgo[loc]
     alf_results[i].objname = objname
     nl_mgn[i].objname = objname
     nl_mgo[i].objname = objname     

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
     nl_mgn[i].feh = nl1.feh
     nl_mgn[i].feh84 = nl1.fehupper
     nl_mgn[i].feh16= nl1.fehlower
     nl_mgn[i].logage = alog10(nl1.age)
     nl_mgn[i].logage84 = alog10(nl1.ageupper)
     nl_mgn[i].logage16= alog10(nl1.agelower)
     nl_mgn[i].mgfe = nl1.alphafe
     nl_mgn[i].mgfe84 = nl1.alphafeupper
     nl_mgn[i].mgfe16= nl1.alphafelower
 
     nl_mgo[i].feh = nl2.feh
     nl_mgo[i].feh84 = nl2.fehupper
     nl_mgo[i].feh16= nl2.fehlower
     nl_mgo[i].logage = alog10(nl2.age)
     nl_mgo[i].logage84 = alog10(nl2.ageupper)
     nl_mgo[i].logage16= alog10(nl2.agelower)
     nl_mgo[i].mgfe = nl2.alphafe
     nl_mgo[i].mgfe84 = nl2.alphafeupper
     nl_mgo[i].mgfe16= nl2.alphafelower
 
  endfor
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1
   sunsym = sunsymbol()
   sym = [14,15,16,17,22,24,18]
   col = ['saddlebrown', 'firebrick','darkorchid','ryb8','goldenrod','pink','orangered']
   colmgn = 'slateblue'
   colmgo = 'crimson'
   psname= 'compare_alf_feh.eps'
   device, filename = psname,xsize = 12,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,alf_results.feh,nl_mgn.feh,psym=1,/nodata,xtitle='[Fe/H] ALF',ytitle='[Fe/H] MPFIT',xrange=[-0.3,0.3],yrange=[-0.3,0.3],xstyle=1,ystyle=1
      cgerrplot,alf_results.feh,nl_mgn.feh16,nl_mgn.feh84
      cgerrplot,nl_mgn.feh,alf_results.feh16,alf_results.feh84,/horizontal
      ;cgerrplot,nl_mgn.feh,alf_results.feh03,alf_results.feh98,/horizontal,color=fsc_color('gray')
      cgerrplot,alf_results.feh,nl_mgo.feh16,nl_mgo.feh84
      cgerrplot,nl_mgo.feh,alf_results.feh16,alf_results.feh84,/horizontal
      for i=0,nfile-1 do begin
         oplot,[alf_results[i].feh],[nl_mgn[i].feh],psym=cgsymcat(sym[i]),color=fsc_color(colmgn)
         oplot,[alf_results[i].feh],[nl_mgo[i].feh],psym=cgsymcat(sym[i]),color=fsc_color(colmgo)
      endfor

      oplot,!x.crange,!x.crange,color=fsc_color('darkgray'),linestyle=2
   device,/close

   psname= 'compare_alf_age.eps'
   device, filename = psname,xsize = 12,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,alf_results.logage,nl_mgn.logage,psym=1,/nodata,xtitle='log Age ALF',ytitle='log Age MPFIT',xrange=[0,1],yrange=[0,1]
      cgerrplot,alf_results.logage,nl_mgn.logage16,nl_mgn.logage84
      cgerrplot,nl_mgn.logage,alf_results.logage16,alf_results.logage84,/horizontal
      cgerrplot,alf_results.logage,nl_mgo.logage16,nl_mgo.logage84
      cgerrplot,nl_mgo.logage,alf_results.logage16,alf_results.logage84,/horizontal
      for i=0,nfile-1 do begin
         oplot,[alf_results[i].logage],[nl_mgn[i].logage],psym=cgsymcat(sym[i]),color=fsc_color(colmgn)
         oplot,[alf_results[i].logage],[nl_mgo[i].logage],psym=cgsymcat(sym[i]),color=fsc_color(colmgo)
      endfor

      oplot,!x.crange,!x.crange,color=fsc_color('darkgray'),linestyle=2
   device,/close

   psname= 'compare_alf_mgh.eps'
   device, filename = psname,xsize = 12,ysize = 10, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,alf_results.mgfe,nl_mgn.mgfe,psym=1,/nodata,xtitle='[Mg/Fe] ALF',ytitle='[Mg/Fe] MPFIT',xrange=[-0.1,0.5],yrange=[-0.1,0.5],xstyle=1,ystyle=1
      cgerrplot,alf_results.mgfe,nl_mgn.mgfe16,nl_mgn.mgfe84
      cgerrplot,nl_mgn.mgfe,alf_results.mgfe16,alf_results.mgfe84,/horizontal
      cgerrplot,alf_results.mgfe,nl_mgo.mgfe16,nl_mgo.mgfe84
      cgerrplot,nl_mgo.mgfe,alf_results.mgfe16,alf_results.mgfe84,/horizontal
      for i=0,nfile-1 do begin
         oplot,[alf_results[i].mgfe],[nl_mgn[i].mgfe],psym=cgsymcat(sym[i]),color=fsc_color(colmgn)
         oplot,[alf_results[i].mgfe],[nl_mgo[i].mgfe],psym=cgsymcat(sym[i]),color=fsc_color(colmgo)
      endfor
      oplot,!x.crange,!x.crange,color=fsc_color('darkgray'),linestyle=2
   device,/close
stop
end
