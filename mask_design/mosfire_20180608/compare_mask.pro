pro compare_mask
   num_bright1= [102,354,500,501,515,214,394,511,392,512,510,173,509]
   num_bright2= [233,244,147,266,149,101,482,425,66,423]
   num_faints = [403,133,129,337,371,370,120]
   num_faintn = [244,248,147,149,428,66]
   num_faintc = [173,391,393,122,512,214,342,500,519]
   num_faint = {s:num_faints,n:num_Faintn,c:num_faintc}
   num_Bright = {one:num_bright1,two:num_bright2}
   zlimit = 0.8793 ;below this Mgb is not in MOSFIRE
   
   datafile = 'cl1604_candidate_list.fits' 
  ;    list = replicate({num:0,objname:'',bluesn:0.,redsn:0.,zmag:0.,f814wmag:0.,$
  ;                      logmass:0.,z:0.,oiiew:0.,ra:0d,dec:0d},ngals)
   list = mrdfits(datafile,1)

  ;faint stuff
   color = ['darkred','seagreen','darkorchid']
   fcolor = ['indianred','lightseagreen','thistle']
   forientation=[45,135,0.]
   psname='faint_masks.eps'
   massbin = 0.1
   massrange=[10.,11.5]
   symsize = [1.5,1.5,1.5]
   !p.font = 0
   set_plot,'ps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   for i=0,2 do begin
      fnow = num_Faint.(i)
      match,list.num,fnow,sublist,subfnow,count=count
      if count ne n_Elements(fnow) then stop,'cannot find some elements'
      curlist = list(sublist)
      goodz = where(curlist.z gt zlimit,ngoodz,complement=badz,ncomplement=nbadz)
      listnow = curlist(goodz)
      goodmag = where(listnow.zmag le 21.7,ngoodmag,complement=badmag,ncomplement=nbadmag)
      print, color(i),ngoodmag
      print,fnow
      print,curlist.zmag
      print,curlist.logmass
      if i eq 0 then plot,listnow(goodmag).logmass,listnow(goodmag).zmag,xrange=massrange,$
                     xtitle='log mass',/nodata,psym=1,yrange=[20,23]
      oplot,listnow(goodmag).logmass,listnow(goodmag).zmag,psym=cgsymcat(14),color=fsc_color(color[i]),symsize=symsize[i]
      oplot,listnow(badmag).logmass,listnow(badmag).zmag,psym=cgsymcat(4),color=fsc_color(color[i]),symsize=symsize[i]
      if nbadz gt 0 then begin
         listnow = curlist(badz)
         oplot,listnow.logmass,listnow.zmag,psym=cgsymcat(16),color=fsc_color(color[i]),symsize=0.8
      endif

;      if i eq 0 then plothist,listnow(goodmag).logmass,bin=massbin,xrange=massrange,$
;                     xtitle='log mass',yrange=[0,3]
;      plothist,listnow(goodmag).logmass,bin=massbin,/overplot,fcolor=fsc_color(fcolor(i)),$
;               color=fsc_color(color[i]),/fline,forientation=forientation[i]
;      plothist,listnow(badmag).logmass,bin=massbin,/overplot,fcolor=fsc_color(fcolor(i)),$
;               color=fsc_color(color[i]),/fline,forientation=forientation[i]
    endfor
   device,/close
   
   ;bright stuff
   psname='bright_masks.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   for i=0,1 do begin
      fnow = num_bright.(i)
      match,list.num,fnow,sublist,subfnow,count=count
      if count ne n_Elements(fnow) then stop,'cannot find some elements'
      curlist = list(sublist)
      goodz = where(curlist.z gt zlimit,ngoodz,complement=badz,ncomplement=nbadz)
      listnow = curlist(goodz)
      goodmag = where(listnow.zmag lt 21.,ngoodmag,complement=badmag,ncomplement=nbadmag)
      print, color(i),ngoodmag
      if i eq 0 then plot,listnow(goodmag).logmass,listnow(goodmag).zmag,xrange=massrange,$
                     xtitle='log mass',/nodata,psym=1,yrange=[20,23]
      oplot,listnow(goodmag).logmass,listnow(goodmag).zmag,psym=cgsymcat(14),color=fsc_color(color[i]),symsize=symsize[i]
      oplot,listnow(badmag).logmass,listnow(badmag).zmag,psym=cgsymcat(4),color=fsc_color(color[i]),symsize=symsize[i]
      if nbadz gt 0 then begin
         listnow = curlist(badz)
         oplot,listnow.logmass,listnow.zmag,psym=cgsymcat(16),color=fsc_color(color[i]),symsize=0.8
      endif
    endfor
   device,/close
  stop 
end
