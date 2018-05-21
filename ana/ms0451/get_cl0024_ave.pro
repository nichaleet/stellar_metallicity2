pro get_cl0024_ave

;GETTING SCI DATA
   sci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_cl0024_ana_afterrev_final.fits',1)
   probsci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/cl0024_feh_age_probdist_final.fits',1)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Averaging
   massbins = [8.9,9.3,9.6,10.,10.25,10.5,10.9]+0.6
   nbins = n_elements(massbins)-1
   ave_mass    = fltarr(nbins)
   ave_mass_dev= fltarr(nbins)
   ave_feh     = fltarr(nbins) 
   ave_feh_dev = fltarr(nbins)
   bndry_mass = fltarr(nbins*2)
   for i=0,nbins-1 do begin
      if i eq nbins-1 then msel = where(sci.logmstar gt massbins[i],cmsel) else $
      msel = where(sci.logmstar gt massbins[i] and sci.logmstar lt massbins[i+1],cmsel)
      meanmc,sci(msel).feh,probsci(msel).dfeh,probsci(msel).probdfeh,fehmean,sigmam,sigmad,sigmas;,/plot
      ;meanerr,sci(msel).feh,(sci(msel).fehupper-sci(msel).fehlower)*0.5,fehmean,sigmam,sigmad,sigmas		
      ave_feh(i) = fehmean
      ave_feh_dev(i) = sigmas
      ave_mass(i) = mean(sci(msel).logmstar)
      ave_mass_dev(i) = stdev(sci(msel).logmstar)
      bndry_mass[i*2:i*2+1] = [massbins[i],massbins[i+1]]
   endfor
   hifeh = interpol(ave_feh+ave_feh_dev,ave_mass,bndry_mass)
   lofeh = interpol(ave_feh-ave_feh_dev,ave_mass,bndry_mass)
   wtoolow = where(lofeh lt -1.,cwtoolow)
   if cwtoolow ge 1 then lofeh(wtoolow) = -1.

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Get SDSS data
   sdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_sdss_ana_afterev_final.fits',1)
   catsdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sci_sdss_ana.fits',2)
   probsdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sdss_feh_age_probdist.fits',1)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Average the SDSS data
   massbins_sdss = [9,9.8,10.3,10.5,10.7,10.9,11.2,11.5]
   massbins_sdss = [9,9.4,9.8,10.1,10.4,10.7,10.9,11.2,11.5]
   nbins_sdss = n_elements(massbins_sdss)-1
   ave_mass_sdss    = fltarr(nbins_sdss)
   ave_mass_dev_sdss= fltarr(nbins_sdss)
   ave_feh_sdss     = fltarr(nbins_sdss)
   ave_feh_dev_sdss = fltarr(nbins_sdss)
   
   for i=0,nbins_sdss-1 do begin
      msel = where(sdss.logmstar gt massbins_sdss[i] and sdss.logmstar lt massbins_sdss[i+1],cmsel)
      meanmc,sdss(msel).feh,probsdss(msel).dfeh,probsdss(msel).probdfeh,fehmean,sigmam,sigmad,sigmas;,/plot
      ;meanerr,sdss(msel).feh,sdss(msel).feherr,fehmean,sigmam,sigmad,sigmas
      ave_feh_sdss(i) = fehmean
      ave_feh_dev_sdss(i) = sigmas
      ave_mass_sdss(i) = mean(sdss(msel).logmstar)
      ave_mass_dev_sdss(i) = stdev(sdss(msel).logmstar)
      ;if i eq nbins_sdss-1 then stop
   endfor
   bndry_mass_sdss = massbins_sdss
   hifeh_sdss = interpol(ave_feh_sdss+ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
   lofeh_sdss = interpol(ave_feh_sdss-ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
   mtoohi = where(bndry_mass_sdss gt 11.5,cmtoohi)
   if cmtoohi gt 0 then remove, mtoohi, bndry_mass_sdss, hifeh_sdss,lofeh_sdss
   
   ;Average the SDSS catalog data (measurements from G05)
   ave_mass_catsdss    = fltarr(nbins_sdss)
   ave_mass_dev_catsdss= fltarr(nbins_sdss)
   ave_feh_catsdss     = fltarr(nbins_sdss)
   ave_feh_dev_catsdss = fltarr(nbins_sdss)
   
   for i=0,nbins_sdss-1 do begin
   	msel = where(catsdss.logm50 gt massbins_sdss[i] and catsdss.logm50 lt massbins_sdss[i+1],cmsel)
   	meanerr,catsdss(msel).z50,(catsdss(msel).z84-catsdss(msel).z16)/2.,fehmean,sigmam,sigmad,sigmas
   	ave_feh_catsdss(i) = fehmean
   	ave_feh_dev_catsdss(i) = sigmas
   	ave_mass_catsdss(i) = mean(catsdss(msel).logm50)
   	ave_mass_dev_catsdss(i) = stdev(catsdss(msel).logm50)
   	;if i eq nbins_sdss-1 then stop
   endfor
   bndry_mass_catsdss = massbins_sdss
   hifeh_catsdss = interpol(ave_feh_catsdss+ave_feh_dev_catsdss,ave_mass_catsdss,bndry_mass_catsdss)
   lofeh_catsdss = interpol(ave_feh_catsdss-ave_feh_dev_catsdss,ave_mass_catsdss,bndry_mass_catsdss)
   mtoohi = where(bndry_mass_catsdss gt 11.5,cmtoohi)
   if cmtoohi gt 0 then remove, mtoohi, bndry_mass_catsdss, hifeh_catsdss,lofeh_catsdss
  
   ;SAVE THE POLYFILL VALUES
   save,hifeh,lofeh,bndry_mass,hifeh_sdss,lofeh_sdss,bndry_mass_sdss,filename='leetho18a_avevals.sav'

end 
