pro combinespec_ms0451,rematch=rematch,plot=plot
;combine all the spectra of MS0451 observed with DEIMOS
;only take the possible members (z=0.5 to 0.6)
;combine continuum divided spectra with weighted mean
;outputs are in outdir
;there are likely 6 extensions per files
;ext 1: combined continuum divided spectra of boxcar extraction and its lambda, contdivivar, calculated sn
;ext 2: structures of original boxcar extracted spectra prior to combination.
;ext 3: same as ext 1 but with gaussian phot extraction
;ext 4: same as ext 2 but wth gaussian phot extraction
;ext 5; list file containing the matches between different spectra of the objects
;ext 6; catalog info of the object from Moran's catalog, alsi add a field of stellar mass (kcorrect)

   outdir='/scr2/nichal/workspace4/prepspec/combined_uniq_ms0451/'
   if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
   npix = 4096

   if keyword_set(rematch) then begin
      cat = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451master_specz.fits.gz',1)
      badcat = where(cat.ra eq 99., complement=goodcat)
      cat = cat(goodcat)
      list = mrdfits('matched_ms0451obj.fits',1)
      ;match object in the list and catalog 
      ras = list.ra[0]
      decs = list.dec[0]
      spherematch, cat.ra,cat.dec,ras,decs,2./3600.,w1,w2,distance
      cat = cat(w1)
      list = list(w2)
      ;only take the galaxies with z from 0.52 to 0.56
      goodz = where(cat.z gt 0.5 and cat.z lt 0.6,cgoodz)
      print,'total of',cgoodz, ' objects between z=0.5 and z=0.6'
      cat = cat(goodz)
      list = list(goodz)
      distance = distance(goodz)
      ;calculate mass
      mass = getms0451mass(cat)
      stop
      ;writefits
      mwrfits,list,'matched_ms0451obj_members.fits',/create,/silent 
      mwrfits,cat,'matched_ms0451obj_members.fits',/silent
      mwrfits,distance,'matched_ms0451obj_members.fits',/silent
      mwrfits,mass,'matched_ms0451obj_members.fits',/silent
   endif
   list = mrdfits('matched_ms0451obj_members.fits',1)
   cat  = mrdfits('matched_ms0451obj_members.fits',2)  
   distance = mrdfits('matched_ms0451obj_members.fits',3) 
   mass = mrdfits('matched_ms0451obj_members.fits',4)
   nobjs = n_elements(list)
   for iobj=0,nobjs-1 do begin ;loop over the list of uniq gals
      obj = list[iobj]
      ndup = obj.count
      ;deal with objname (if there are two names, merge them)
      if ndup gt 1 then begin
          pos=uniq(strtrim(obj.objname[0:ndup-1],2))
          if n_elements(pos) gt 1 then name =  strjoin(strtrim(obj.objname(pos),2),'_') $
          else name = strtrim(obj.objname[0],2)
      endif else name = strtrim(obj.objname[0],2)
      print,iobj,nobjs-1,format='("now doing ",I3,"/",I3)'
      print, 'specname    catname   redshift   distance'
      print,name,' ',strtrim(cat[iobj].SPECNAME,2),cat[iobj].z,distance[iobj]
      ;;;;;;;;;;;
      fileout = outdir+'spec1d.combinedmasks.'+string(iobj,format='(I03)')+'.'+name+'.fits'

      for iframe=1,3,2 do begin ;loop over extraction type in deimos
         strall=[]
         hdrall=[]
         for idup=0,ndup-1 do begin
             strblu = mrdfits(strtrim(obj.filename[idup],2),iframe,hdrblu,/silent,status=status1)
             strred = mrdfits(strtrim(obj.filename[idup],2),iframe+1,hdrred,/silent,status=status2)
             if status1 ne 0 or status2 ne 0 then stop,'stopped: raw spectrum error'
             if status1 ne 0 or status2 ne 0 and iframe eq 3 then goto, jump1
             if typename(strblu) eq 'INT' then goto, jump1
             strnow = join_deimos_blu_red(strblu,strred)
             hdrnow = hdrblu
             sxaddpar,hdrnow,'extname',strtrim(sxpar(hdrblu,'extname'),2)+' '+strtrim(sxpar(hdrred,'extname'),2)
             if iframe eq 1 and keyword_set(plot) then plotcont=1 else plotcont=0
             continuumnormalize_deimos,strnow,cat[iobj].z,hdrnow,/splinecon,plot=plotcont ;optional ,/polycon or ,/splinecon(default)
             strall=[strall,strnow]
             if idup eq 0 then sizehdr = (size(hdrnow,/dimension))[0]
             if idup gt 0 and (size(hdrnow,/dimension))[0] gt sizehdr then  hdrnow = hdrnow[0:sizehdr-1]
             if idup gt 0 and (size(hdrnow,/dimension))[0] lt sizehdr then  hdrnow = [hdrnow,strarr(sizehdr-(size(hdrnow,/dimension))[0])]
             hdrall=[[hdrall],[hdrnow]]
         endfor

         if ndup eq 1 then begin
             strout={lambda:strall.lambda,dlam:strall.dlam,contdiv:strall.contdiv,contdivivar:strall.contdivivar,sn:strall.sn,exptime:total(obj.exptime)} 
         endif         

         if ndup gt 1 then begin ;stack continuum and make output structure (strout)
             sn = strall.sn
             ;shift velocity to heliocentric frame
             for idup=0,ndup-1 do begin   
                 hdrnow = hdrall[*,idup]
                 lambda = strnow.lambda
                 jd = double(sxpar(hdrnow,'mjd-obs'))+2400000.5D 
                 vhelio = helio_deimos(obj.ra[idup],obj.dec[idup],2000,jd=jd)
                 lambda = (1.+vhelio/3e5)*lambda
                 strall[idup].lambda = lambda  ;ok it's shifted
                 print, '   vhelio, sn, exptime:', vhelio, sn(idup), obj.exptime(idup)
             endfor   

             ;find the total wl range (lambdafull)
             lambdaarr = strall.lambda ;[npix,ndup]
             minwl = min(lambdaarr,wminwl)
             maxwl = max(lambdaarr,wmaxwl)
             wminwl = array_indices(lambdaarr,wminwl)
             wmaxwl = array_indices(lambdaarr,wmaxwl)
             maxsn = max(sn,ref)

             ;if overlapped, use the lambda of the spec with max sn as ref
             ;attach the extra lambdas to the main lambdas
             lambdafull = lambdaarr[*,ref]
             if wmaxwl[1] ne ref then begin 
                 wadd = where(lambdaarr[*,wmaxwl[1]] gt max(lambdafull),cwadd)
                 if cwadd gt 0 then lambdafull = [lambdafull,lambdaarr[wadd,wmaxwl[1]]]
             endif
             if wminwl[1] ne ref then begin
                 wadd = where(lambdaarr[*,wminwl[1]] lt min(lambdafull),cwadd)
                 if cwadd gt 0 then lambdafull = [lambdaarr[wadd,wminwl[1]],lambdafull]
             endif
             ;get data in and interpolate
             npix = n_elements(lambdafull)
             strout={lambda:lambdafull,dlam:fltarr(npix),contdiv:fltarr(npix),contdivivar:fltarr(npix),sn:0.,exptime:0.}
             contdivarr=fltarr(npix,ndup)+1./0.
             contdivivararr = fltarr(npix,ndup)+1./0.
             dlamarr = fltarr(npix,ndup)+1./0.
             framearr = intarr(npix,ndup)
             for k=0,ndup-1 do begin
                strnow = strall[k]
                loc = value_locate(lambdafull,strnow.lambda)
                if loc[0] eq -1 then stop,'this should not have happened'
                contdivarr[loc,k]=interpol(strnow.contdiv,strnow.lambda,lambdafull[loc])
                contdivivararr[loc,k]=interpol(strnow.contdivivar,strnow.lambda,lambdafull[loc])
                dlamarr[loc,k]=interpol(strnow.dlam,strnow.lambda,lambdafull[loc])
                framearr[loc,k] +=1
             endfor

             ;weighted average
             nframearr = total(framearr,2)
             wnodup = where(nframearr eq 1, cwnodup)
             if cwnodup gt 0 then begin
                 strout.dlam(wnodup) = total(dlamarr[wnodup,*],2,/nan)
                 strout.contdiv(wnodup) = total(contdivarr[wnodup,*],2,/nan)
                 strout.contdivivar(wnodup) = total(contdivivararr[wnodup,*],2,/nan)
             endif
             wdup = where(nframearr gt 1, cwdup)
             for ii=0,cwdup-1 do begin                 
                 wnan = where(~finite(contdivivararr[wdup[ii],*]),cnan)
                 dcontdiv = 1./sqrt(abs(contdivivararr[wdup[ii],*]))
                 if cnan gt 0 then dcontdiv(wnan) = 1./0.
                 err = abs(dcontdiv/contdivarr[wdup[ii],*])
                 strout.dlam[wdup[ii]] = wmean(dlamarr[wdup[ii],*],abs(err*dlamarr[wdup[ii],*]),/nan) 
                 strout.contdiv[wdup[ii]]=wmean(contdivarr[wdup[ii],*],dcontdiv,error=contdivivar,/nan)
                 strout.contdivivar[wdup[ii]]=contdivivar
             endfor
             if keyword_set(plot) then begin
                colors = ['purple','cyan','green','yellow','orange','red','brown']
                z=cat[iobj].z
                plot,strout.lambda/(1.+z),strout.contdiv,xtitle='lambda',ytitle='contdiv',title=name,yrange=[-4,4]
                for iplot=ndup-1,0,-1 do oplot,strall[iplot].lambda/(1.+z),strall[iplot].contdiv,color=fsc_color(colors[iplot])
                oplot,strout.lambda/(1.+z),strout.contdiv,color=fsc_color('white')
                wait,1
             endif    
         endif ;if ndup gt 1

         ;write data
         if iframe eq 1 then create=1 else create=0
         if ndup gt 1 then hdrout=hdrall[*,0] else hdrout=hdrall
         sxaddpar,hdrout,'exptime',total(obj.exptime)
         strout.exptime = total(obj.exptime)
         mwrfits,strout,fileout,hdrout,create=create,/silent
         mwrfits,strall,fileout
         jump1: continue
      endfor ;for each frame 
      mwrfits,obj,fileout
      catout = create_struct(cat[iobj],'mass_kcorrect',mass[iobj])
      mwrfits,catout,fileout
  endfor ;loop each object
end
