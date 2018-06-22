pro make_targetlist

;From visually inspecting current spectra, (see page 55-57 in NL's logbook
;The list are then created and the priorities are set by current SN and magnitude

   listnum = [49,58,66,86,94,101,102,113,120,122,126,127,129,133,147,148,149,155,162,173,187,$
           188,197,202,207,214,230,233,239,244,248,266,278,283,286,287,289,291,296,298,$
           307,312,328,331,337,342,343,345,354,357,370,371,382,385,389,391,392,393,394,$
           395,399,403,412,417,423,425,428,436,482,497,500,501,509,510,511,512,515,519]
   highpnum = [497,501,509,510,515]
   
   zlimit = 0.8793 ;below this Mgb is not in MOSFIRE
   
   datafile = 'cl1604_candidate_list.fits' 
   if file_test(datafile) eq 0 then begin
      ;get data in
      ngals = n_Elements(listnum)
      list = replicate({num:0,objname:'',bluesn:0.,redsn:0.,zmag:0.,f814wmag:0.,$
                        logmass:0.,z:0.,oiiew:0.,ra:0d,dec:0d},ngals)
      for i=0,ngals-1 do begin
         file = file_search('/scr2/nichal/workspace4/prepspec_cl1604/data/spline/*_'+$
                             strtrim(string(listnum[i]),2)+'_*.fits.gz',count=cfile)
         if cfile ne 1 then stop,'found more than 1 file',file
         aa = mrdfits(file[0],1,/silent)
         list[i].num = listnum[i]
         list[i].objname = aa.phot_id
         wlris = where(strmatch(aa.indiv_instru,'*lris*'),clris)
         wdeimos = where(strmatch(aa.indiv_instru,'*deimos*'),cdeimos)
         wblue = []
         if clris gt 0 then wblue=[wblue,wlris]
         if cdeimos gt 0 then wblue = [wblue,wdeimos]
         cblue = clris+cdeimos
         if cblue gt 0 then list[i].bluesn = sqrt(total((aa(wblue).indiv_sn)^2)) 
         wmosfire = where(strmatch(aa.indiv_instru,'*mosfire*'),cmosfire)
         if cmosfire gt 0 then list[i].redsn = sqrt(total((aa(wmosfire).indiv_sn)^2))
         list[i].zmag = aa.zmag
         list[i].f814wmag = aa.f814wmag
         list[i].z = aa.z
         list[i].logmass = aa.LOGMSTAR_SED_BL
         list[i].oiiew = aa.ew_oii
         list[i].ra = aa.ra
         list[i].dec = aa.dec
      endfor
      mwrfits,list,datafile,/silent,/create
  endif
  list = mrdfits(datafile,1)

  goodz = where(list.z gt zlimit,complement=badz,ncomplement=nlowz)
  listlowz = list(badz)
  list = list(goodz)
  radec,list.ra,list.dec,ihr,imin,xsec,ideg,imn,xsc
  radec,listlowz.ra,listlowz.dec,ihrlowz,iminlowz,xseclowz,ideglowz,imnlowz,xsclowz
  ;shorter exposure
  bright = where(list.zmag lt 21.,cbright,complement=faint,ncomplement=cfaint)
  brightpri = fix(list(bright).bluesn*2000)
  match,list(bright).num,highpnum,suba,subb,count=chighp
  if chighp gt 0 then brightpri(suba) = max(brightpri)+500
  write_ds9_regionfile, list(bright).ra, list(bright).dec,$
       comment=strtrim(string(list(bright).num),2), filename='cl1604_brightobj.reg',color='red'
  writecol,'cl1604_brightobj.txt',strtrim(string(list(bright).num),2)+'_'+$
           list(bright).objname,brightpri,list(bright).zmag,$
           ihr(bright),imin(bright),xsec(bright),ideg(bright),imn(bright),xsc(bright),$
           fltarr(cbright)+2000.,fltarr(cbright)+2000.,fltarr(cbright),fltarr(cbright),$
           fmt='(A,1x,I5,2x,F4.1,2x,I2,2x,I2,2x,F5.2,2x,I2,2x,I2,2x,F5.2,2x,F6.1,2x,F6.1,2x,F3.1,2x,F3.1)'
  ;longer exposure
  faintpri = fix((25.-list(faint).zmag)*1000.+(list(faint).bluesn*100))
  match,list(faint).num,highpnum,suba,subb,count=chighp
  if chighp gt 0 then faintpri(suba) = max(brightpri)
  write_ds9_regionfile, list(faint).ra, list(faint).dec,$
       comment=strtrim(string(list(faint).num),2), filename='cl1604_faintobj.reg',color='green'
  writecol,'cl1604_faintobj.txt',strtrim(string(list(faint).num),2)+'_'+$
           list(faint).objname,faintpri,list(faint).zmag,$
           ihr(faint),imin(faint),xsec(faint),ideg(faint),imn(faint),xsc(faint),$
           fltarr(cfaint)+2000.,fltarr(cfaint)+2000.,fltarr(cfaint),fltarr(cfaint),$
           fmt='(A,1x,I5,2x,F4.1,2x,I2,2x,I2,2x,F5.2,2x,I2,2x,I2,2x,F5.2,2x,F6.1,2x,F6.1,2x,F3.1,2x,F3.1)'
  ;filler galaxies
  lowzpri = fix((25.-listlowz.zmag)*100)
  write_ds9_regionfile, listlowz.ra, listlowz.dec,$
       comment=strtrim(string(listlowz.num),2), filename='cl1604_lowzobj.reg',color='blue'
  writecol,'cl1604_lowzobj.txt',strtrim(string(listlowz.num),2)+'_'+$
           listlowz.objname,lowzpri,listlowz.zmag,$
           ihrlowz,iminlowz,xseclowz,ideglowz,imnlowz,xsclowz,$
           fltarr(nlowz)+2000.,fltarr(nlowz)+2000.,fltarr(nlowz),fltarr(nlowz),$
           fmt='(A,1x,I4,2x,F4.1,2x,I2,2x,I2,2x,F5.2,2x,I2,2x,I2,2x,F5.2,2x,F6.1,2x,F6.1,2x,F3.1,2x,F3.1)'

  stop 
end
