function make_chisqarr_worm2,refa,refz,refmg,sn
   common degen, fehgrid, agegrid, mgfegrid, degen_file ,nfeh,nage,nmgfe 
   common misc, feharr,agearr,mgfearr,deltafeh,deltaage,deltamgfe
   common chistuff, dirchi
   ageref = agegrid(refa)
   fehref = fehgrid(refz)
   mgferef = mgfegrid(refmg)

   ;find where is the spec closest to the input age and feh and alpha
   min_dist = min(sqrt((fehgrid-fehref)^2+(agegrid-ageref)^2+(mgfegrid-mgferef)^2),iref)
   if min_dist ne 0. then print,strtrim(string(min_dist))+' matching might be off grid'
   loc = array_indices(fehgrid,iref)

   ;make noisy ref spectra
   specarr = mrdfits(degen_file(loc[1]),loc[2]+1,/silent)
   refspecstr = specarr(loc[0])
   if refspecstr.feh ne fehgrid(iref) or refspecstr.age ne agegrid(iref) or refspecstr.alpha ne mgfegrid(iref) $
        then stop,'grid does not match with file'
   npix = n_elements(refspecstr.lambda)
   refspec_err = abs(refspecstr.spec/sn)
   refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

   ;make chisqarr
   gridsize = size(fehgrid,/dimensions)
   chisqarr = fltarr(gridsize)

   for ia = 0,n_elements(degen_file)-1 do for iap = 1,gridsize[2] do begin
      specarr = mrdfits(degen_file(ia),iap,/silent)
      if gridsize[0] ne n_elements(specarr) then stop,'grid does not match size'
      for i=0,n_elements(specarr)-1 do begin
         loc = where(fehgrid eq specarr[i].feh and agegrid eq specarr[i].age and $
               mgfegrid eq specarr[i].alpha,cloc)
         if cloc ne 1 then stop,'oops'
         ind = array_indices(fehgrid,loc)
         chisqarr[ind[0],ind[1],ind[2]] = total((refspec-specarr[i].spec)^2/refspec_err)
      endfor
   endfor
   
   ;write fits file
   mkhdr,header,chisqarr
   paraname = ['CTYPE1','CTYPE2','CTYPE3']
   paraval =['[Fe/H]','Age','[Mg/Fe]']
   for ii=0,2 do sxaddpar,header,paraname(ii),paraval(ii)
   paraname = ['CRVAL1','CDELT1','CRPIX1','CROTA1',$
	   'CRVAL2','CDELT2','CRPIX2','CROTA2',$
	   'CRVAL3','CDELT3','CRPIX3','CROTA3']
   paraval = [fehgrid(0),fehgrid(1)-fehgrid(0),0.,0.,$
	      agegrid(0),agegrid(1)-agegrid(0),0.,0.,$
	      mgfegrid(0),mgfegrid(1)-mgfegrid(0),0.,0.]
   for ii =0, 11 do sxaddpar,header,paraname(ii),paraval(ii)
   sxaddpar,header,'AGE',ageref
   sxaddpar,header,'FEH',fehref
   sxaddpar,header,'MGFE',mgferef

   if file_test(dirchi,/directory) eq 0 then file_mkdir,dirchi
   filechi = dirchi+'/chisq_age'+strtrim(string(fix(refa)),2)+'_feh'+$
           strtrim(string(fix(refz)),2)+'_mgfe'+strtrim(string(fix(refmg)),2)+$
           '.fits'
   writefits,filechi,chisqarr,header
   print,'wrote ',filechi
   return,chisqarr
end

pro initialize_degenfiles, degendir
   common degen, fehgrid, agegrid, mgfegrid, degen_file,nfeh,nage,nmgfe 
   common misc, feharr,agearr,mgfearr,deltafeh,deltaage,deltamgfe

   degenfile = file_search(degendir+'age*_sspdegen_gridspec_ms0451.fits',count=cdegenfile)
   for i=0,cdegenfile-1 do begin
      degen_sps = mrdfits(degenfile(i),1,/silent)
      if i eq 0 then begin
           fits_info,degenfile(i),/silent,n_ext=nmgfe
           nfeh = n_elements(degen_sps)
           nage = cdegenfile
           fehgrid = fltarr(nfeh,nage,nmgfe)
           agegrid = fltarr(nfeh,nage,nmgfe)
           mgfegrid = fltarr(nfeh,nage,nmgfe)
           degen_file = strarr(nage)
           for j=0,nmgfe-1 do begin
              degen_sps = mrdfits(degenfile(i),j+1,/silent)
              mgfegrid[*,*,j] = rebin(degen_sps.alpha,nfeh,nage)
           endfor
      endif
      filebase = file_basename(degenfile(i),'.fits')
      agei = fix(strmid(filebase,3,2))  ;0 to 50 

      degen_file[agei] = degenfile(i)
      fehgrid[*,agei,*] = rebin(degen_sps.feh,nfeh,nmgfe)
      agegrid[*,agei,*] = rebin(degen_sps.age,nfeh,nmgfe)
   endfor
   feharr = reform(fehgrid[*,0,0])
   agearr = reform(agegrid[0,*,0])
   mgfearr = reform(mgfegrid[0,0,*]) 

   deltafeh = feharr[1]-feharr[0]
   deltaage = agearr[1]-agearr[0]
   deltamgfe = mgfearr[1]-mgfearr[0]
end

pro plotworm_feh_age,posmg,outname
   common degen, fehgrid, agegrid, mgfegrid, degen_file ,nfeh,nage,nmgfe 
   common misc, feharr,agearr,mgfearr,deltafeh,deltaage,deltamgfe
   common ref, ref_ia, ref_iz, ref_img, ageref, fehref, mgferef,snset
   common chistuff, dirchi

   set_plot,'ps'
   !p.font=0
   burd_colors
   device, filename = outname,xsize = 15,ysize = 10, $
   	xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,findgen(100)/125.-0.6,alog10(findgen(126)/10.+0.5)+9.,xrange=[-0.6,0.2],yrange=alog10([2,12.])+9,/nodata,$
   	xtitle='[Fe/H]',xstyle=1,ystyle=5
   axis, yaxis=0,ytitle='log(Age/yr)',yrange=alog10([2,12.])+9,ystyle=1
   axis, yaxis=1,yticks=6,ytickv=alog10([2,3,4,5,6,8,11.])+9.,ytickn=['2','3','4','5','6','8','11'],ycharsize=0.7
   fehcolor=[10,60,100,190,254]

   for i=0,n_Elements(ref_ia)-1 do for j=0,n_elements(ref_iz)-1 do begin

      filechi = dirchi+'/chisq_age'+strtrim(string(fix(ref_ia(i))),2)+'_feh'+$
           strtrim(string(fix(ref_iz(j))),2)+'_mgfe'+strtrim(string(fix(posmg)),2)+$
           '.fits'
      if file_test(filechi) then chisqarr = readfits(filechi) else $ 
           chisqarr = make_chisqarr_worm2(ref_ia(i),ref_iz(j),posmg,snset)
      
      ;calculate probability from chisqarr
      Lgrid = -0.5*(chisqarr-min(chisqarr))
      probgrid = exp(double(Lgrid))

      ;marginalized over mgfe
      probfehage = fltarr(nfeh,nage)
      for ii=0,nfeh-1 do for jj=0,nage-1 do $
	      probfehage[ii,jj] = int_tabulated(mgfegrid[ii,jj,*],probgrid[ii,jj,*])

      ;check volume
      volume = total(probfehage)*deltafeh*deltaage
      probfehage = probfehage/volume 
      print, volume 
      ;shift midage and midfeh to the center of the fehref and ageref
     ; shiftfeh = fehref-midfeh
     ; shiftfeh = fehref-maxprobfeh
     ; shiftage = ageref-midage
     ; shiftage = ageref-maxprobage

      ;contour,probgrid,fehgrid+shiftfeh,alog10(agegrid+shiftage)+9.,levels=[maxprob/1.65],c_linestyle=0,/overplot,color=fehcolor(i)
      maxprob = max(probfehage)
      contour,probfehage,fehgrid,alog10(agegrid)+9.,levels=[maxprob/1.65],c_linestyle=0,/overplot,color=fehcolor(i) ;sqrt(e)=1.65 
      oplot,[fehref[j]],alog10([ageref[i]])+9.,psym=7,color=fehcolor(j)   	
      print,'plotted ',fehref[j],alog10([ageref[i]])+9.
   endfor
   device,/close
   stop
end

pro plotworm_feh_alpha,posage,outname
   common degen, fehgrid, agegrid, mgfegrid, degen_file ,nfeh,nage,nmgfe 
   common misc, feharr,agearr,mgfearr,deltafeh,deltaage,deltamgfe
   common ref, ref_ia, ref_iz, ref_img, ageref, fehref, mgferef,snset
   common chistuff, dirchi

   print,'Making worm plot for feh vs Mgfe'
   print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
   set_plot,'ps'
   !p.font=0
   burd_colors
   device, filename = outname,xsize = 15,ysize = 10, $
   	xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,findgen(100)/125.-0.6,alog10(findgen(126)/10.+0.5)+9.,xrange=[-0.6,0.2],yrange=alog10([2,12.])+9,/nodata,$
   	xtitle='[Fe/H]',xstyle=1,ystyle=5
   axis, yaxis=0,ytitle='log(Age/yr)',yrange=alog10([2,12.])+9,ystyle=1
   axis, yaxis=1,yticks=6,ytickv=alog10([2,3,4,5,6,8,11.])+9.,ytickn=['2','3','4','5','6','8','11'],ycharsize=0.7
   fehcolor=[10,60,100,190,254]

   for i=0,n_Elements(ref_img)-1 do for j=0,n_elements(ref_iz)-1 do begin

      chisqarr = make_chisqarr_worm2(posage,ref_iz(j),ref_img(i),snset)
      
      ;calculate probability from chisqarr
      Lgrid = -0.5*(chisqarr-min(chisqarr))
      probgrid = exp(double(Lgrid))

      ;marginalized over age
      probfehmgfe = fltarr(nfeh,nmgfe)
      for ii=0,nfeh-1 do for jj=0,nage-1 do $
	      probfehmgfe[ii,jj] = int_tabulated(agegrid[ii,*,jj],probgrid[ii,*,jj])

      ;check volume
      volume = total(probfehmgfe)*deltafeh*deltamgfe
      probfehmgfe = probfehmgfe/volume 
      print, volume 
      ;shift midage and midfeh to the center of the fehref and ageref
     ; shiftfeh = fehref-midfeh
     ; shiftfeh = fehref-maxprobfeh
     ; shiftage = ageref-midage
     ; shiftage = ageref-maxprobage

      ;contour,probgrid,fehgrid+shiftfeh,alog10(agegrid+shiftage)+9.,levels=[maxprob/1.65],c_linestyle=0,/overplot,color=fehcolor(i)
      maxprob = max(probfehmgfe)
      contour,probfehmgfe,fehgrid,mgfegrid,levels=[maxprob/1.65],c_linestyle=0,/overplot,color=fehcolor(i) ;sqrt(e)=1.65 
      oplot,[fehref[j]],alog10([ageref[i]])+9.,psym=7,color=fehcolor(j)   	
   endfor
   device,/close
end

pro make_worm_plot2,sn=sn
   common degen, fehgrid, agegrid, mgfegrid, degen_file ,nfeh,nage,nmgfe 
   common misc, feharr,agearr,mgfearr,deltafeh,deltaage,deltamgfe
   common ref, ref_ia, ref_iz, ref_img, ageref, fehref, mgferef,snset
   common chistuff, dirchi

   ;READ THE SSP FILES AND SET THE GRIDS
   degendir = '/scr2/nichal/workspace4/sps_fit/sspdegen/sspdegen_ms0451/'
   initialize_degenfiles,degendir

   ;SETTING
   ;find chisquare at each reference point
   ref_ia = [8,12,16,21,25] ;~2,3,4,5,6 Gyr 
   ref_iz = [22,28,34,40,46] ; -0.45,-0.3,-0.15,0.0,0.15 dex
   ref_img = [3,10,17,23,30] ; -0.3,-0.1,0.1,0.3,0.5 dex

   ageref = agearr(ref_ia)
   fehref = feharr(ref_iz)
   mgferef = mgfearr(ref_img) 
 
   outdir = '/scr2/nichal/workspace4/ana/wormplot/output/'
   dirchi='chisqoutputsn'+strtrim(string(fix(sn)),2)
   if ~keyword_set(sn) then sn = 15.   
   snset = float(sn)

   ;1)plot between feh and age at fixed alpha = 0.1 
   outname=outdir+'/worm_plot_feh_age_sn'+strtrim(string(sn),2)+'.eps'
   posmg = ref_img[2] ; should be where [mg/fe] = 0.1
   plotworm_feh_age,posmg,outname
end
