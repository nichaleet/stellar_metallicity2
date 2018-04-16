pro find_matching0451
   ;search reduced deimos spectra and match them according to both name and RA and DEC (distance lt 2 arcsec)
   ;output is structure matched_ms0451obj.fits

   files = file_search('/scr2/nichal/keck/deimos/Cl0024MS0451/nicha_reduced/ms0451/masks/*/spec1d*.{fits,fits.gz}', count=cspec)
   ;first list all the files and sort according to RA 
   mask = []
   slit = []
   ra = []
   dec = []
   objname = []
   filename = []
   exptime = []
   for i=0,cspec-1 do begin
      hdr = headfits(files[i],exten=1)
      ira = sxpar(hdr,'ra_obj')
      idec = sxpar(hdr,'dec_obj')
      if strmatch(ira,'0.0*') eq 0 then begin
      ;only keeps those with ra and dec, which means not including serendips
         ra = [ra,ira]
         if strmatch(idec,'-\**') eq 1 then begin
            decobj = strsplit(idec,':',/extract)
            decmask = strsplit(sxpar(hdr,'dec'),':',/extract)
            idec = strjoin([decmask[0],decobj[1:2]],':')
         endif
         dec = [dec,idec]
         basefile = file_basename(files[i])
         extensions = strsplit(basefile, '.', /extract)
         objname = [objname,extensions[3]]
         mask = [mask,extensions[1]]
         slit = [slit,extensions[2]]
         filename = [filename,files[i]]
         curexptime = double(sxpar(hdr,'exptime'))
         exptime = [exptime,curexptime]
      endif ;else print, sxpar(hdr,'objno'),' is not included'
   endfor 
    
   order = sort(ra)
   mask = mask(order)
   slit = slit(order)
   rastr = ra(order)
   decstr = dec(order)
   objname = objname(order)
   filename = filename(order)
   exptime = exptime(order)

   ;convert ra and dec into decimal
   nspec = n_Elements(rastr)
   ra = fltarr(nspec)
   dec = fltarr(nspec)
   for i=0,nspec-1 do begin
     raarr = float(strsplit(rastr(i),':',/extract))
     decarr = float(strsplit(decstr(i),':',/extract))
     ra(i) = (raarr(0)+raarr(1)/60.+raarr(2)/3600.)*15.
     dec(i) = signum(decarr(0))*(abs(decarr(0))+decarr(1)/60.+decarr(2)/3600.)
   endfor  
   
   ;find matching
   ndup = 7
   objlist = {objname:strarr(ndup),rastring:strarr(ndup),decstring:strarr(ndup),$
              ra:fltarr(ndup),dec:fltarr(ndup),mask:strarr(ndup),slit:strarr(ndup),$
              filename:strarr(ndup),exptime:dblarr(ndup),count:0}
   objlist = replicate(objlist,1500) 
   foundmatch = bytarr(nspec)
   iobj = 0L
   for i=0,nspec-1 do begin
      if foundmatch(i) eq 1 then continue ;only do those with no match yet
      ;1) match by name
      wnamematch = where(strmatch(objname,'*'+objname(i)+'*'),cnamematch)
      ;2) match by ra dec
      gcirc,2,ra(i),dec(i),ra,dec,dis
      wlocmatch = where(dis lt 2.,clocmatch)
      ;3) union the two matches
      wmatch = setunion(wnamematch,wlocmatch)
      cmatch = n_elements(wmatch)
      foundmatch(wmatch)=1 ;flag the match
      ;keep the matches
      objlist(iobj).objname[0:cmatch-1]=objname(wmatch)
      objlist(iobj).ra[0:cmatch-1]=ra(wmatch)
      objlist(iobj).dec[0:cmatch-1]=dec(wmatch)
      objlist(iobj).rastring[0:cmatch-1]=rastr(wmatch)
      objlist(iobj).decstring[0:cmatch-1]=decstr(wmatch)
      objlist(iobj).mask[0:cmatch-1]=mask(wmatch)
      objlist(iobj).slit[0:cmatch-1]=slit(wmatch)
      objlist(iobj).filename[0:cmatch-1]=filename(wmatch)
      objlist(iobj).exptime[0:cmatch-1]=exptime(wmatch)
      objlist(iobj).count = cmatch
      iobj +=1
   endfor   
    
   objlist = objlist[0:iobj-1]
   mwrfits,objlist,'matched_ms0451obj.fits',/create,/silent
   objlist = mrdfits('matched_ms0451obj.fits',1)
   wmult = where(objlist.count gt 1,cmult)
   for i=0,cmult-1 do begin
      nspec = objlist(wmult(i)).count
      print, objlist(wmult(i)).objname[0:nspec-1]
      print, objlist(wmult(i)).rastring[0:nspec-1]
      print, objlist(wmult(i)).decstring[0:nspec-1]
    endfor

   stop
end
