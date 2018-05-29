pro addobject,matchpos,cmatchpos,file, instru, mask, slit
common catfile, cat, objfound
common outfile, objall, maxdup, minz, maxz
   if cmatchpos gt 1 then begin
       print, file+' has '+ string(cmatchpos,format='(I3)')+' matches in the catalog'
       print, 'mask: ', cat(matchpos).mask
       print, 'slit: ', cat(matchpos).slit
       print, 'Choose the first match....'
       matchpos = matchpos[0]
       cmatchpos = 1
   endif
   if cmatchpos eq 0 then print,file+' has no matched object in the catalog'
   if cmatchpos eq 1 and cat(matchpos).z ge minz and cat(matchpos).z le maxz then begin
      if objfound(matchpos) eq -999 then begin
         objnum = n_elements(objall)
         objfound(matchpos) = objnum
         obj = create_struct(name='obj',cat(matchpos),'files',strarr(maxdup),'instrus',strarr(maxdup),$
                             'masks',strarr(maxdup),'slits',strarr(maxdup),'ndup',0)
         obj.files[obj.ndup] = file
         obj.instrus[obj.ndup] = instru
         obj.masks[obj.ndup] = mask
         obj.slits[obj.ndup] = slit
         obj.ndup += 1
         objall = [objall,obj]
      endif else begin
         objnum = objfound(matchpos)
         obj = objall[objnum]
         obj.files[obj.ndup] = file
         obj.instrus[obj.ndup] = instru
         obj.masks[obj.ndup] = mask
         obj.slits[obj.ndup] = slit
         obj.ndup += 1
         objall[objnum] = obj
      endelse
   endif
end

pro find_matching1604
common catfile, cat, objfound
common outfile, objall, maxdup, minz, maxz
;Look at all existing spectra, match them based on ID/slit name/mask name with the catalog, and list all spectra of the same object under one object name. Only save those with redshfit 0.8<z<1.0
;A lot of this is copied from and old code in /scr2/nichal/workspace2/catalogs/orelse_sample_ana.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   maxdup = 5
   minz = 0.8
   maxz = 1.0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;reading in data
   ;catalog from the ORELSE group
   cat = mrdfits('/scr2/nichal/workspace2/catalogs/Cl1604.fits.gz',1,hdr)
   ;old deimos and lris spectra
   odfile = file_search('/scr2/nichal/keck/deimos/Cl1604/deimos/spec1d*.fits.gz',count=codfile) 
   olfile = file_search('/scr2/nichal/keck/deimos/Cl1604/lris/spec1d*.fits.gz',count=colfile)  
   ;new lris spectra 
   nlfile1 = file_search('/scr2/nichal/WMKO/LRIS/2017may_cl1604/reduce/cl1604_1*.r.*.msdc.fits*',count=cnlfile1)
   nlfile2 = file_search('/scr2/nichal/WMKO/LRIS/2018may_cl1604/reduce/cl1604_1*.r.*.msdc.fits*',count=cnlfile2)
   nlfile = nlfile1+nlfile2
   cnlfile = cnlfile1+cnlfile2
   ;new mosfire spectra
   mfile =  file_search('/scr2/nichal/keck/mosfire/2014jun/Cl1604_*/onedspec/spec1d*.fits',count=cmfile)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Match each spectrum with the catalog
   catmask = strtrim(cat.mask,2)
   catslit = strtrim(cat.slit,2)
   catobjname = strtrim(cat.phot_id,2)
   ncat = n_elements(cat)
   objfound = lonarr(ncat)-999
   objall = []
;1) old deimos
   for i=0,codfile-1 do begin
      filebase = file_basename(odfile[i])
      extensions = strsplit(filebase,'.',/extract)
      slit = extensions[2]
      mask = extensions[1]
      matchpos = where(catmask eq mask and catslit eq slit,cmatchpos)
      addobject,matchpos,cmatchpos,odfile[i],'deimos_old',mask,slit
   endfor
;2) old lris
  for i=0,colfile-1 do begin
     filebase = file_basename(olfile[i])
     extensions = strsplit(filebase,'.',/extract)
     if n_elements(extensions) eq 6 then slit = extensions[2]
     if n_elements(extensions) eq 7 then slit = extensions[2]+'.'+extensions[3]
     mask = extensions[1]
     if strmid(mask,0,3) ne 'old' then mask='LRIS.'+mask
     if mask eq 'oldLRISd' then mask = 'oldLRIS'
     if strmid(slit,0,1) eq '0' then slit =  strmid(slit,1)
     matchpos = where(catmask eq mask and catslit eq slit,cmatchpos)
     addobject,matchpos,cmatchpos,olfile[i],'lris_old',mask,slit
  endfor
;3 new lris
   for i=0,cnlfile-1 do begin
      filebase = file_basename(nlfile[i])
      extensions = strsplit(filebase,'.',/extract)
      objname = extensions[2]
      mask = extensions[0]
      slit = ''
      matchpose = where(catobjname eq objname,cmatchpos)
      addobject,matchpose,cmatchpos,nlfile[i],'lris',mask,slit
   endfor
;4) new mosfire
   for i=0, cmfile-1 do begin
      filebase = file_basename(mfile[i])
      extensions = strsplit(filebase,'.',/extract)
      mask = extensions[1]
      objname = extensions[3]
      matchpose = where(catobjname eq objname, cmatchpos)
      slit = ''
      addobject,matchpose,cmatchpos,mfile[i],'mosfire',mask,slit
   endfor
   print, 'Total found ', n_Elements(objall), ' members.' 
   mwrfits,objall,'matched_cl1604obj_members.fits',/create,/silent

end
