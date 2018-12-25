pro cl1604_spec_info
;for now, this is to reply email from Brian Lemaux asking to trade MOSFIRE data
   redodata = 0

   if redodata eq 1 then begin
      files = file_search('/scr2/nichal/workspace4/prepspec_cl1604/data/spline/spline*.fits.gz',count=cfiles)
      ndup=4
      mostemp={objnum:'', $
               objname:'', $
               sn:-999d, $
               ra:-999d, $
               dec:-999d, $
               indiv_sn:fltarr(ndup),$
               indiv_ra:dblarr(ndup), $
               indiv_dec:dblarr(ndup), $
               indiv_mask:strarr(ndup),$
               indiv_slit:strarr(ndup),$
               indiv_filename:strarr(ndup),$
               indiv_good:bytarr(ndup), $
               indiv_count:-999L,$
               indiv_instru:strarr(ndup),$
               Z:-999d,$
               LOGMSTAR_SED_BL:-999d,$
               LOGMSTAR_SED_BLERR:-999d,$
               EW_OII:-999d,$
               EW_OIIERR:-999d, $
               SNMOSFIRE:0.}
      mosdataarr = replicate(mostemp,cfiles)
      goodarr = bytarr(cfiles)
 
      for i=0,cfiles-1 do begin
         spec = mrdfits(files[i],1,/silent)
         if spec.indiv_count gt 1 and total(strmatch(spec.indiv_instru,'mosfire*')) ge 1 then begin
            locmos = where(strmatch(spec.indiv_instru,'mosfire*') eq 1, clocmos)
            if clocmos gt 1 then stop
            mosdata = mostemp
            struct_assign,spec,mosdata
            mosdata.snmosfire = spec.indiv_sn[locmos]
            mosdata.ra = spec.indiv_ra[locmos]
            mosdata.dec = spec.indiv_dec[locmos]
            mosdataarr[i] = mosdata
            goodarr[i] = 1
         endif 
      endfor
      good = where(goodarr eq 1,cgood)
      mosdataarr = mosdataarr(good)
      mwrfits,mosdataarr,'mosfire_cl1604_info.fits',/silent,/create
   endif else mosdataarr = mrdfits('mosfire_cl1604_info.fits',1)
   stop
   writecol,'mosfire_cl1604_info.txt',mosdataarr.objname,mosdataarr.ra, mosdataarr.dec,mosdataarr.logmstar_sed_bl,mosdataarr.snmosfire,mosdataarr.sn,mosdataarr.ew_oii
   stop
end
