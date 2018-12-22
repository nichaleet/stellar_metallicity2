
pro leetho2018_cl0024ana_forjobtalk
;mainly for plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;GETTING SCI DATA
   sci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_cl0024_ana_afterrev_final.fits',1)
   probsci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/cl0024_feh_age_probdist_final.fits',1)
   ageform = (galage(sci.zfit,1000.)/1.e9-sci.age)>0. ;age of universe when it was formed
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Get SDSS data
   sdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_sdss_ana_afterev_final.fits',1)
   catsdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sci_sdss_ana.fits',2)
   probsdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sdss_feh_age_probdist.fits',1)
   sdss_ageform = (galage(sdss.zfit,1000.)/1.e9-sdss.age)>0. ;age of universe when it was formed	
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;restore average values
   restore,'/scr2/nichal/workspace4/ana/ms0451/leetho18a_avevals.sav'
   ;bndry_mass,hifeh,lowfeh, bndry_mass_sdss,hifeh_sdss,lofeh_sdss
   wtoolow = where(lofeh lt -1.,cwtoolow)
   if cwtoolow ge 1 then lofeh(wtoolow) = -1.
   mtoohi = where(bndry_mass_sdss gt 11.5,cmtoohi)
   if cmtoohi gt 0 then remove, mtoohi, bndry_mass_sdss, hifeh_sdss,lofeh_sdss
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;Average the SDSS catalog data (measurements from G05)
   massbins_sdss = [9,9.4,9.8,10.1,10.4,10.7,10.9,11.2,11.5]
   nbins_sdss = n_elements(massbins_sdss)-1
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
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;GETTING LITERATURE VALUES
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Get the average values measured in Gallazzi05
   g05_mass = [9.00,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5];the first time is actually 8.91
   g05_feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
   ;below are the stars in the figure 8 of Gallazzi05
   g05_feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
   g05_fehlo = g05_feh-g05_feherr
   g05_fehhi = g05_feh+g05_feherr
   
   ;below are the diamonds in figure 8 of Gallazzi05
   ;g05_fehlo= [-1.11,-1.07,-1.1,-1.03,-0.97,-0.9,-0.8,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04]
   ;g05_fehhi= [0.0,0.0,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28]
   toolow = where(g05_fehlo lt -1.,ctoolow)
   if ctoolow gt 0 then g05_fehlo(toolow) = -1

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Choi's Data
   z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01],mgfe:[0.04,0.13,0.16,0.21,0.23],mgfeerr:[0.05,0.03,0.01,0.01,0.02]}
   z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02],mgfe:[0.17,0.19,0.20,0.22,0.23],mgfeerr:[0.06,0.02,0.01,0.02,0.04]}
   z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02],$
          mgfe:[0.18,0.22,0.21,0.29],mgfeerr:[0.04,0.02,0.01,0.03]}
   z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02],$
          mgfe:[0.23,0.23,0.30],mgfeerr:[0.04,0.02,0.03]}
   z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03],$
          mgfe:[0.05,0.09,0.19],mgfeerr:[0.13,0.05,0.04]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Conroy2014 data (SDSS)
   conroy = {mass:[9.63,9.75,9.80,10.08,10.55,10.70,11.07],Feh:[-0.07,-0.05,-0.04,-0.03,-0.01,-0.02,0.],$
             MgFe:[.05,0.08,0.09,0.12,0.15,0.20,0.22]}

   ;Xiangcheng's FIRE data
   readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0pt8.txt',xma_mass08,xma_feh08
   readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0.txt',xma_mass0,xma_feh0
   xma_feh08 = xma_feh08-0.2
   xma_feh0 = xma_feh0-0.2
   testma = 0
   if testma eq 1 then begin ; test for evolution
      good0 = where(xma_mass0 gt 9.3,cgood0)
      good08 = where(xma_mass08 gt 9.3,cgood08)
      ma_linfit_z0 = linfit(xma_mass0(good0)-10.,xma_feh0(good0),sigma=ma_linfit_z0_err)
      ma_linfit_z08 = linfit(xma_mass08(good08)-10.,xma_feh08(good08),sigma=ma_linfit_z08_err)
      ;slopes are not different
      ma_commonslope = wmean([ma_linfit_z0(1),ma_linfit_z08(1)],[ma_linfit_z0_err(1),ma_linfit_z08_err(1)])
      ;linfit with a fixed slope of 0.56684088 dex per log mass
      ma_linfit_fixslope_z0 = [ma_linfit_z0(0),ma_commonslope]
      yfitma0 = curvefit(xma_mass0(good0)-10.,xma_feh0(good0),fltarr(cgood0)+1,ma_linfit_fixslope_z0,ma_linfit_fixslope_z0_err,function_name='linearfit',fita=[1,0])
      ma_linfit_fixslope_z08 = [ma_linfit_z08(0),ma_commonslope]
      yfitma08 = curvefit(xma_mass08(good08)-10.,xma_feh08(good08),fltarr(cgood08)+1,ma_linfit_fixslope_z08,ma_linfit_fixslope_z08_err,function_name='linearfit',fita=[1,0])
      ma_zscore = (ma_linfit_fixslope_z0(0)-ma_linfit_fixslope_z08(0))/sqrt(ma_linfit_fixslope_z0_err(0)^2+ma_linfit_fixslope_z08_err(0)^2)
      print,'ma evolution = ',ma_linfit_fixslope_z0(0)-ma_linfit_fixslope_z08(0),'p/m',sqrt(ma_linfit_fixslope_z0_err(0)^2+ma_linfit_fixslope_z08_err(0)^2)
      print,'z score=',ma_zscore   
   ;  stop 
   endif
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;De Rossi 2017 (EAGLE)
   ;read off from Figure 5
   DeRossi_z0={mass:[9.15,9.48,9.81,10.15,10.5,10.72],feh:[-0.17,-0.08,0.05,0.12,0.17,0.31],feherr:[0.075,0.08,0.07,0.07,0.07,0.03]}
   DeRossi_z1={mass:[9.18,9.47,9.85,10.15,10.55],feh:[-0.39,-0.28,-0.09,0.08,0.17],feherr:[0.06,0.06,0.08,0.07,0.085]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Sybilska et al 2017 (hELENa, IFU from Sauron data)
   aa=read_csv('/scr2/nichal/workspace2/catalogs/sybilska.csv',n_table_header=1,header=header)
   Syb_z0 = {mass:aa.field1,feh:aa.field2}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Lu 2014
   Lu_z0 = {mass:[8.21,8.61,9.0,9.41,9.8,10.2,10.61,11],feh:[-0.53,-0.43,-0.35,-0.26,-0.16,-0.098,-0.075,-0.09],$
            feherr:[0.13,0.12,0.11,0.12,0.11,0.10,0.095,0.095]}
   Lu_z1 = {mass:[8.21,8.61,9,9.4,9.8,10.2,10.6,11],feh:[-0.56,-0.44,-0.35,-0.25,-0.18,-0.13,-0.11,-0.11],$
             feherr:[0.14,0.13,0.11,0.11,0.11,0.10,0.10,0.09]}
   Lu_somerville_z0 = {mass:[8.25,8.65,9.05,9.45,9.85,10.25,10.64],feh:[-1.00,-0.82,-0.64,-0.46,-0.33,-0.15,-0.03],$
                      feherr:[0.06,0.1,0.1,0.1,0.1,0.10,0.095]}
   Lu_lu_z0 = {mass:[8.73,8.95,9.18,9.47,9.76,10.09,10.46,10.88],feh:[-.99,-0.84,-0.70,-0.52,-0.31,-0.06,0.18],$
              feherr:[0.1,0.1,0.1,0.1,0.1,0.1,0.10,0.095]}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   redofit = 0
   if redofit then begin
   ;Test low-mass end slope
   ;1 degree polynomial (linear equation - linfit)
   sdlm = where(sdss.logmstar lt 10.5)
   sclm = where(sci.logmstar lt 10.5)
   lin_mzr_lowmass = linfitmc([sci(sclm).logmstar,sdss(sdlm).logmstar]-10.,[sci(sclm).feh,sdss(sdlm).feh],$
                     [[probsci(sclm).dfeh],[probsdss(sdlm).dfeh]],[[probsci(sclm).probdfeh],[probsdss(sdlm).probdfeh]],$
                     lin_mzr_lowmass_err,chisq=chilin_lowmass)

   lin_mzr_sdss_lowmass = linfitmc(sdss(sdlm).logmstar-10.,sdss(sdlm).feh,probsdss(sdlm).dfeh,probsdss(sdlm).probdfeh,$
                          lin_mzr_sdss_lowmass_err,chisq=chilin_sdss_lowmass,yfit=yfitsdss_lowmass)

   lin_mzr_cl_lowmass = linfitmc(sci(sclm).logmstar-10.,sci(sclm).feh,probsci(sclm).dfeh,probsci(sclm).probdfeh,$
                        lin_mzr_cl_lowmass_err,chisq=chilin_cl_lowmass,yfit=yfitcl_lowmass)

   print,'best linear fit for lowmass end (combined):',sigfig([lin_mzr_lowmass,lin_mzr_lowmass_err],3)
   print,'best linear fit for lowmass end  (sdss)    :',sigfig([lin_mzr_sdss_lowmass,lin_mzr_sdss_lowmass_err],3)
   print,'best linear fit for lowmass end (z=0.4)   :',sigfig([lin_mzr_cl_lowmass,lin_mzr_cl_lowmass_err],3)

   ;Find the best fit MZR parameters
   ;1) functional form from Zahid13
   ;montecarlo
   zahid_mzr = fitzahidmc([sci.logmstar,sdss.logmstar],[sci.feh,sdss.feh],[[probsci.dfeh],[probsdss.dfeh]],$
               [[probsci.probdfeh],[probsdss.probdfeh]],zahid_mzr_sigma,chisq=chizahid)
   zahid_mzr_sdss = fitzahidmc(sdss.logmstar,sdss.feh,probsdss.dfeh,probsdss.probdfeh,zahid_mzr_sdss_sigma,chisq=chizahid_sdss)
   zahid_mzr_cl = fitzahidmc(sci.logmstar,sci.feh,probsci.dfeh,probsci.probdfeh,zahid_mzr_cl_sigma,chisq=chizahid_cl,guess_param=[8.95,10.,0.6])

   ;2) 1 degree polynomial (linear equation - linfit)
   lin_mzr = linfitmc([sci.logmstar,sdss.logmstar]-10.,[sci.feh,sdss.feh],$
               [[probsci.dfeh],[probsdss.dfeh]],[[probsci.probdfeh],[probsdss.probdfeh]],lin_mzr_err,chisq=chilin)

   lin_mzr_sdss = linfitmc(sdss.logmstar-10.,sdss.feh,probsdss.dfeh,probsdss.probdfeh,lin_mzr_sdss_err,chisq=chilin_sdss,yfit=yfitsdss)

   lin_mzr_cl = linfitmc(sci.logmstar-10.,sci.feh,probsci.dfeh,probsci.probdfeh,lin_mzr_cl_err,chisq=chilin_cl,yfit=yfitcl)

   print,'best linear fit (combined):',sigfig([lin_mzr,lin_mzr_err],3)
   print,'best linear fit (sdss)    :',sigfig([lin_mzr_sdss,lin_mzr_sdss_err],3)
   print,'best linear fit (z=0.4)   :',sigfig([lin_mzr_cl,lin_mzr_cl_err],3)

   ;normal z test for slope and intercept
   ;slope
   zslope = (lin_mzr_sdss(1)-lin_mzr_cl(1))/sqrt(lin_mzr_sdss_err(1)^2+lin_mzr_cl_err(1)^2)
   print, 'normal z value for slope of SDSS and Cl0024 =',zslope
   ;The z value is ~-0.5 which is p=0.3. The slopes are not different. So can proceed to do intercept according to
   ;http://www.biostathandbook.com/ancova.html 
   ;find the weighted slope
   commonslope = (lin_mzr_cl(1)/lin_mzr_cl_err(1)^2+lin_mzr_sdss(1)/lin_mzr_sdss_err(1)^2)$
                 /(1./lin_mzr_cl_err(1)^2+1./lin_mzr_sdss_err(1)^2)
   commonslope_err = 1./sqrt(1./lin_mzr_cl_err(1)^2+1./lin_mzr_sdss_err(1)^2)
   commonslope_stdev = sqrt(((lin_mzr_cl(1)-commonslope)^2/lin_mzr_cl_err(1)^2+$
                             (lin_mzr_sdss(1)-commonslope)^2/lin_mzr_sdss_err(1)^2)/$
                            (1./lin_mzr_cl_err(1)^2+1./lin_mzr_sdss_err(1)^2)/2.)

   print, 'common slope',commonslope,commonslope_err,commonslope_stdev
   ;fit with common slope
   sdss_fixslope_pars = [lin_mzr_sdss(0),commonslope]
   yfit_sdss_fixslope = linfit_fixedslope_mc(sdss.logmstar-10.,sdss.feh,probsdss.dfeh,probsdss.probdfeh,$
                        sdss_fixslope_pars,sdss_fixslope_pars_err,chisq=chisq_sdss_fixslope)

   cl_fixslope_pars = [lin_mzr_cl(0),commonslope]
   yfit_cl_fixslope = linfit_fixedslope_mc(sci.logmstar-10.,sci.feh,probsci.dfeh,probsci.probdfeh,$
                        cl_fixslope_pars,cl_fixslope_pars_err,chisq=chisq_cl_fixslope)

   print, 'Fixed slope intercepts:'
   print, 'SDSS:',sdss_fixslope_pars(0),sdss_fixslope_pars_err(0)
   print, 'z=0.4:',cl_fixslope_pars(0),cl_fixslope_pars_err(0)
   ;mean of all mass
   meanmass=mean([sci.logmstar,sdss.logmstar]-10.)
   meanmass = 0. ;set to fit at 10^10 Msun
   print, 'mean mass:',meanmass
   sdss_feh_meanmass = commonslope*meanmass+sdss_fixslope_pars(0)
   sdss_feh_meanmass_err = sqrt(commonslope_err^2*meanmass^2+sdss_fixslope_pars_err(0)^2)
   cl_feh_meanmass = commonslope*meanmass+cl_fixslope_pars(0)
   cl_feh_meanmass_err = sqrt(commonslope_err^2*meanmass^2+cl_fixslope_pars_err(0)^2)
   ;z test comparing two means
   zscore = ((sdss_feh_meanmass-cl_feh_meanmass)-0.)/sqrt(sdss_feh_meanmass_err^2+cl_feh_meanmass_err^2)
   print, 'z score for intercepts:',zscore

   ;3) 2 degree polynomial
   poly_mzr =  polyfitmc([sci.logmstar,sdss.logmstar],[sci.feh,sdss.feh],2,[[probsci.dfeh],[probsdss.dfeh]],$
               [[probsci.probdfeh],[probsdss.probdfeh]],poly_mzr_err,chisq=chipoly)
   poly_mzr_sdss = polyfitmc(sdss.logmstar,sdss.feh,2,probsdss.dfeh,probsdss.probdfeh,poly_mzr_sdss_err,chisq=chipoly_sdss)
   poly_mzr_cl = polyfitmc(sci.logmstar,sci.feh,2,probsci.dfeh,probsci.probdfeh,poly_mzr_cl_err,chisq=chipoly_cl)

   ;4) linear fit with intrinsic scatter
   bestparam_sdss = intrinsic_scatter_linear_nongauss(sdss.logmstar-10.,sdss.feh,$
                    [sdss_fixslope_pars,0.06],probsdss.dfeh,probsdss.probdfeh,fita=[1,1,0])
   print, 'intercept, slope, int scatter(sdss)',sigfig(bestparam_sdss[0:2],3)
   bestparam_cl = intrinsic_scatter_linear_nongauss(sci.logmstar-10.,sci.feh,$
                  [cl_fixslope_pars,0.06],probsci.dfeh,probsci.probdfeh,fita=[1,1,0])
   print, 'intercept, slope, int scatter(z=0.4)',sigfig(bestparam_cl[0:2],3)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Use Akaike's information criterion
   samplesize = n_elements([sci.logmstar,sdss.logmstar])
   k=2
   aic_lin = chilin+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=3
   aic_poly= chipoly+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   aic_zahid = chizahid*(samplesize-3)+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   print, 'AIC linear = ', aic_lin
   print, 'AIC 2degre polynomial = ', aic_poly
   print, 'AIC Zahid =',aic_zahid

   print,'AIC Cl0024:linear(fixed slope),linear, polynomial'
   samplesize=n_Elements(sci.logmstar)
   k=1
   print,chisq_cl_fixslope+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=2
   print,chilin_cl+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=3
   print,chipoly_cl+2.*k+(2.*k*(k+1.))/(samplesize-k-1)

   print,'AIC SDSS:linear(fixed slope),linear, polynomial'
   samplesize=n_Elements(sdss)
   k=1
   print,chisq_sdss_fixslope+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=2
   print,chilin_sdss+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=3
   print,chipoly_sdss+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   save, zahid_mzr, zahid_mzr_cl,zahid_mzr_SDSS,lin_mzr, lin_mzr_err, lin_mzr_sdss,lin_mzr_sdss_err,$
         lin_mzr_cl,lin_mzr_cl_err,commonslope,commonslope_err,commonslope_stdev,$
         sdss_fixslope_pars,sdss_fixslope_pars_err,cl_fixslope_pars,cl_fixslope_pars_err,$
         poly_mzr,poly_mzr_sdss,poly_mzr_cl,poly_mzr_err,poly_mzr_sdss_err,poly_mzr_cl_err,$
         bestparam_sdss,bestparam_cl,filename='mzr_obs_param.sav'

   endif
   restore,'mzr_obs_param.sav' 
   mzr_mass = findgen(104)/40.+8.9
   mzr_zahid= zahid_mzr[0]-alog10(1.+10.^((mzr_mass-zahid_mzr[1])*(-1.*zahid_mzr[2])))-8.9
   mzr_poly = lin_mzr[0]+lin_mzr[1]*(mzr_mass-10.)
   mzr_poly_sdss = sdss_fixslope_pars[0]+sdss_fixslope_pars[1]*(mzr_mass-10.)
   mzr_poly_cl = cl_fixslope_pars[0]+cl_fixslope_pars[1]*(mzr_mass-10.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PLOTTING
   ;misc
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='SDSS_FeH_massv3.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[9,11.5]
      yrange=[-1.,0.35]
      plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      
      ;shade the average region
      polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=fsc_color('pbg2')
     ; polyfill,[bndry_mass_catsdss,reverse(bndry_mass_catsdss)],[hifeh_catsdss,reverse(lofeh_catsdss)],$
     ;          color=fsc_color('pbg1')
      
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 8.9)) = 9.
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      rainbow_colors,n_colors=21
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;draw data points
      cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=40,symsize=1.3
      ;plot representative uncertainties
      lowermass = where(catsdss.logm50 lt 10.,complement=uppermass)
      cgerrplot,[-0.32,-0.32],[11.1+median(catsdss(lowermass).logm16-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm16-catsdss(uppermass).logm50)],$
                              [11.1+median(catsdss(lowermass).logm84-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm84-catsdss(uppermass).logm50)],$
                               /horizontal,color=40,thick=2
      cgerrplot,[11.1,11.3],[-0.32+median(sdss(lowermass).fehlower-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehlower-sdss(uppermass).feh)],$
                            [-0.32+median(sdss(lowermass).fehupper-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehupper-sdss(uppermass).feh)],$
                             color=40,thick=2

      ;Add DeRossi17 
;      oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=fsc_color('springgreen'),$
;      linethick=2,errcolor=fsc_color('springgreen'),psym=1
;      oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),symsize=1.5,color=fsc_color('springgreen')
      
      ;Add Xiangcheng's data
;      oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=fsc_color('forestgreen'),symsize=1.2
      
       ;Add Choi's data
       oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color('red5'),linethick=2,errcolor=10
       oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=fsc_color('red5'),symsize=2

       ;Add Conroy's
       oplot,conroy.mass,conroy.feh,psym=cgsymcat(16),color=fsc_color('ygb6'),symsize=2      
      ;Add Sybilska2017
;      oplot,syb_z0.mass,syb_z0.feh,color=fsc_color('maroon'),thick=3,linestyle=5

      ;Add Gallazzi2014 best linear fit to quiescent galaxies
      massG14 = [10,10.5,11,11.5]
      ZG14 = 0.15*(massG14-11)+0.109
      rmsG14 = 0.14
      oplot,massg14,zg14,color=fsc_color('maroon'),thick=3,linestyle=5
      
      ;Add best fitted line
      oplot,mzr_mass,mzr_poly_sdss,color=0,thick=3 ;linear function
      oplot,mzr_mass,zahid_mzr[0]-alog10(1.+10.^((mzr_mass-zahid_mzr[1])*(-1.*zahid_mzr[2])))-8.9,color=0,linestyle=1,thick=3
;      oplot,mzr_mass,poly_mzr_sdss[0]+poly_mzr_sdss[1]*mzr_mass+poly_mzr_sdss[2]*mzr_mass^2,color=0,linestyle=3,thick=3
       ;Labelling
;      cglegend,title=['Quiescent SDSS','measured in this work','Gallazzi et al. 2005','Choi et al. 2014','Sybilska et al. 2017'],psym=[14,0,15,46,0],location=[10.6,-0.5],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=['navy','navy','pbg1','red5','maroon'],symsize=1.2

;      cglegend,title=['Ma et al. 2016','De Rossi et al. 2017'],psym=[16,24],location=[10.6,-0.87],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=['forestgreen','springgreen'],symsize=1.2
      cglegend,title=['Quiescent SDSS','measured in this work','Gallazzi et al. 2005','Choi et al. 2014','Gallazzi et al. 2014'],psym=[14,0,15,46,0],location=[10.6,-0.5],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=['navy','navy','pbg2','red5','maroon'],symsize=1.2

      oplot,[10.53,10.67],[-0.8,-0.8],linestyle=2,color=fsc_color('maroon'),thick=2
      oplot,[10.6],[-0.645],psym=cgsymcat(15),color=fsc_Color('pbg2'),symsize=1.3
      xyouts,10.1,-0.52,'observations:',charsize=0.8
      xyouts,10.1,-0.88,'simulations:',charsize=0.8
      xyouts,9.1,0.2,'z~0 galaxies',charsize=1.
      ;bracket,10.45,-0.8,10.5,-0.5,/left
      ;bracket,10.45,-0.95,10.5,-0.85,/left
   device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   cleanplot,/silent
   !p.font=0
   psname='Cl0024_FeH_mass2.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[9.5,11.5]
      yrange=[-1.,0.3]
      plot,sci.logmstar,sci.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      ;shade the average region
      ;colorarr=[10,50,203,230,254]
      colorarr= fsc_color(['royalblue','darkorchid','deeppink','maroon','red8'])
      x=[bndry_mass,reverse(bndry_mass)]
      y=[hifeh,reverse(lofeh)]
      polyfill,x,y,color=fsc_color('rose')
      
      ylefthi = interpol(hifeh_sdss,bndry_mass_sdss,[xrange[0]])
      yleftlo = interpol(lofeh_sdss,bndry_mass_sdss,[xrange[0]])
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      goodx = where(x gt xrange[0])
      x= x(goodx)
      y= y(goodx)
      x= [xrange[0],x,xrange[0]]
      y= [ylefthi,y,yleftlo]
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;Add Xiangcheng's data
      ;oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=colorarr(0),symsize=1.2
      ;oplot,xma_mass08,xma_feh08,psym=cgsymcat(16),color=colorarr(4),symsize=1.2
      
      ;Add Munoz15
      
      ;Add DeRossi17 
      ;oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0)
      ;oploterror,derossi_z1.mass,derossi_z1.feh,derossi_z1.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4)      
      ;oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.1
      ;oplot,derossi_z1.mass,derossi_z1.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.1
      
      ;draw data points
      cgerrplot,sci.logmstar,sci.fehlower,sci.fehupper,color='pink',thick=1
      oplot,sci.logmstar,sci.feh,psym=cgsymcat(14),color=colorarr(2),symsize=1.2
      oplot,sci.logmstar,sci.feh,psym=cgsymcat(4),color=fsc_color('maroon'),symsize=1.2
      
      ;Add Choi's data
      vsym,5,/fill,/star
      ;oploterror,z01.mass,z01.feh,z01.feherr,color=colorarr(1),linethick=2,errcolor=colorarr(1)
;      oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color('org7'),linethick=2,errcolor=fsc_color('org7')
      ;oploterror,z06.mass,z06.feh,z06.feherr,color=colorarr(3),linethick=2,errcolor=colorarr(3)
      ;oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=colorarr(1),symsize=2
;      oplot,z03.mass,z03.feh,psym=cgsymcat(46),color=colorarr(2),symsize=1.7
      ;oplot,z06.mass,z06.feh,psym=cgsymcat(46),color=colorarr(3),symsize=2
      ;oplot,z01.mass,z01.feh,psym=cgsymcat(45),symsize=2
;      oplot,z03.mass,z03.feh,psym=cgsymcat(45),symsize=1.7
      ;oplot,z06.mass,z06.feh,psym=cgsymcat(45),symsize=2
      
      ;Add Lu14
      ;oploterror,lu_z0.mass,lu_z0.feh,lu_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0),errthick=2
      ;oploterror,lu_z1.mass,lu_z1.feh,lu_z1.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4)               
      oplot,lu_z0.mass,lu_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.4
      oplot,lu_z0.mass,lu_z0.feh,color=colorarr(0)
      ;oplot,lu_z1.mass,lu_z1.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.2
      oplot,lu_somerville_z0.mass,lu_somerville_z0.feh,psym=cgsymcat(16),color=fsc_color('darkgreen'),symsize=1,linestyle=0
      oplot,lu_somerville_z0.mass,lu_somerville_z0.feh,color=fsc_color('darkgreen')
      oplot,lu_lu_z0.mass,lu_lu_z0.feh,psym=cgsymcat(22),color=fsc_color('darkorchid'),symsize=1,linestyle=0
      oplot,lu_lu_z0.mass,lu_lu_z0.feh,color=fsc_color('darkorchid')
      ;add best fitted MZR curve
      ;oplot,mzr_mass,mzr_zahid,linestyle=2
     ; oplot,mzr_mass,mzr_poly,thick=2
      oplot,mzr_mass,mzr_poly_sdss,thick=3,color=fsc_color('navy'),linestyle=2
      oplot,mzr_mass,mzr_poly_cl,thick=3,color=fsc_color('org7'),linestyle=2
      ;Labelling
      ;al_legend,['z=0','z=[0.1,0.2]','z=[0.3,0.4]','z=[0.6,0.7]','z=[0.8,1]'],psym=15,color=colorarr,box=0,position=[10.6,-0.6]
;      al_legend,['z=0','z=[0.3,0.4]'],psym=15,color=['royalblue','deeppink'],box=0,position=[10.9,-0.75]
      al_legend,['Current work','Croton model','Somerville model','Lu model'],color=['deeppink','royalblue','darkgreen','darkorchid'],psym=[14,24,16,22],position=[10.9,-0.67],box=0,charsize=0.9,symsize=[1.3,1.4,1,1]		
      xyouts,9.6,0.2,'z~0.4 galaxies',charsize=1.
   device,/close
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;TEST KS TEST

   ;ks_test,sdss.logmstar,sdss.feh,sci.logmstar,sci.feh,D_ks,Neff_ks,pvalue   
   ; stop  
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Gas phase metal from Zahid2013 (12+log(O/H))
   gas_oh = [{redshift:0.08,z0:9.121,M0:8.999,gmma:0.85},{redshift:0.29,z0:9.130,M0:9.304,gmma:0.77},$
             {redshift:0.78,z0:9.161,M0:9.661,gmma:0.65},{redshift:1.4,z0:9.06,M0:9.6,gmma:0.7},$
             {redshift:2.26,z0:9.06,M0:9.7,gmma:0.6}]
   ;the equation is 12+log(O/H) = z0-log[1+(M*/M0)^-gmma]
   sun_oh = 8.8
;alpha_Fe = 0.12*mass-1.1 ;linear plot by eye to Choi14 Fig 8, mass is in log scale
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='Formation_Redshift_cl0024.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      	xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      zcat = [1000,2,1,0.7,0.4,0.]
      ;	zcat = [1000.,1.5,0.7,0]
      agecat = galage(zcat,1000.)/1.e9 ;age of universe at that redshif	
      plot,sci.logmstar,sci.feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],$
          xstyle=1,yrange=[-0.8,0.2],/nodata
      
      nage = n_Elements(agecat)-1
      rainbow_colors
      zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
      ;loop over formation times
      for i=0,nage-1 do begin
      	selsdss = where(sdss_ageform gt agecat(i) and sdss_ageform le agecat(i+1),cselsdss)
      	if cselsdss gt 0 then begin
      		cgplot,sdss(selsdss).logmstar,sdss(selsdss).feh,psym=16,/overplot,color=zcolor(i),symsize=0.7
      	endif
      endfor
      for i=0,nage-1 do begin
              sel = where(ageform gt agecat(i) and ageform le agecat(i+1), csel)
              if csel gt 0 then begin
                cgerrplot,sci(sel).logmstar,sci(sel).fehlower,sci(sel).fehupper,color=zcolor(i)-10,thick=0.5
      		cgplot,sci(sel).logmstar,sci(sel).feh,psym=14,/overplot,color=zcolor(i),symsize=1.3
      		cgplot,sci(sel).logmstar,sci(sel).feh,psym=4,/overplot,color='darkgray',symsize=1.3
              endif
      endfor
      ;add gas phase MZR
;      marr = findgen(101)/40.+9. ;logmass from 9 to 11.5
;      a_fe = 0.12*marr-1.1
;      for i=0,n_Elements(gas_oh)-1 do begin
;      	o_h = gas_oh(i).z0-alog10(1.+10.^((marr-gas_oh(i).m0)*(-1.)*gas_oh(i).gmma))-sun_oh;-a_fe
;      	colornowi = value_locate(zcat,gas_oh(i).redshift)
;     	colornow = zcolor(colornowi)
;      	oplot,marr,o_h,color=colornow,thick=0.5
;      endfor
      oplot,mzr_mass,mzr_poly_sdss,thick=2,linestyle=2
      oplot,mzr_mass,mzr_poly_cl,thick=2,linestyle=2
      xyouts,[11.2,11.2],[0.,0.13],['z~0.4','z~0'],charsize=1
      ;Labelling
      zarr_str = strarr(n_elements(zcat)-1)
      for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+$
                                    '<z$\tex_{form}$<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
      zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(F3.1)'),2)
      al_Legend,zarr_str,psym=15,color=zcolor,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
      al_Legend,['SDSS subsample','Cl0024 z~0.4'],psym=[16,14],symsize=[0.5,1.3],color=0,box=0,thick=2,$
                 charsize=1,position=[10.6,-0.3],font=0
      xyouts,9.1,0.1,'(a)',charsize=1.2
   device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='FEH_deviation_fromz0line_obs_probdist.eps'
   if skip then sdss_fixslope_pars=[-0.0528328,0.158635]
   fehideal = sdss_fixslope_pars[0]+(sci.logmstar-10.)*sdss_fixslope_pars[1]
   sdss_fehideal = sdss_fixslope_pars[0]+(sdss.logmstar-10.)*sdss_fixslope_pars[1]
   x=[sdss_ageform,ageform]
   y=[sdss.feh-sdss_fehideal,sci.feh-fehideal]
   dx = 0.5*[sdss.ageupper-sdss.agelower,sci.ageupper-sci.agelower]
   dy = 0.5*[sdss.fehupper-sdss.fehlower,sci.fehupper-sci.fehlower]
   dxarr = -1.*[[probsdss.dage],[probsci.dage]] ;since x axis is ageuni-agegal
   dyarr = [[probsdss.dfeh],[probsci.dfeh]]
   probdxdyarr = transpose([[[probsdss.probdfehdage]],[[probsci.probdfehdage]]],[1,0,2]) ;age feh nobjs
   ;;method 1 linfitexymc
   ;feh_dev_par = linfitexymc(x,y,dx,dy,dxarr,dyarr,probdxdyarr,sigma_a_b,paramarr,/fitwith_linfit)
   ;feh_dev_par_sdss = linfitexymc(sdss_ageform,sdss.feh-sdss_fehideal,0.5*[sdss.ageupper-sdss.agelower],0.5*[sdss.fehupper-sdss.fehlower],-1.*probsdss.dage,probsdss.dfeh,transpose(probsdss.probdfehdage,[1,0,2]),sigma_a_b_sdss,paramarr_sdss,/fitwith_linfit)
   ;feh_dev_par_cl = linfitexymc(ageform,sci.feh-fehideal,0.5*[sci.ageupper-sci.agelower],0.5*[sci.fehupper-sci.fehlower],-1.*probsci.dage,probsci.dfeh,transpose(probsci.probdfehdage,[1,0,2]),sigma_a_b_cl,paramarr_cl,/fitwith_linfit,miny=-0.9,/plot)
   ;;method 2
   ;fitexy,ageform,sci.feh-fehideal,Acl,Bcl,x_sig=0.5*[sci.ageupper-sci.agelower],y_sig=0.5*[sci.fehupper-sci.fehlower],sigma_a_b_cl
   ;feh_dev_par_cl = [Acl,Bcl]
   ;;method 3  Evan's ab96 code /raid/idl/enk/ab96.pro from Akritas and Bershady 1996
   ;feh_dev_par = ab96(x,y,dx,dy,3,sigma=sigma_a_b)
   ;xcl = ageform
   ;ycl = sci.feh-fehideal
   ;dxcl = 0.5*[sci.ageupper-sci.agelower]
   ;dycl = 0.5*[sci.fehupper-sci.fehlower]
   ;feh_dev_par_cl = ab96(xcl,ycl,dxcl,dycl,3,sigma=sigma_a_b_cl) 
   ;feh_dev_par_sdss = ab96(sdss_ageform,sdss.feh-sdss_fehideal,0.5*[sdss.ageupper-sdss.agelower],0.5*[sdss.fehupper-sdss.fehlower],3,sigma=sigma_a_b_sdss)  
   ;;method 4 Bayesian. Summing prob from each cells
   redolinfit = 1
   ;stop
   if redolinfit eq 1 then begin
      ;feh_dev_par_sdss = linfitprobgrid(sdss_ageform,sdss.feh-sdss_fehideal,-1.*probsdss.dage,probsdss.dfeh,transpose(probsdss.probdfehdage,[1,0,2]),sigma_a_b_sdss,outfitfile='sdss_prob_deltafeh_age.fits',outepsfile='sdss_plot_deltafeh_age.eps',fitrange=[0,11.5],/plot,badid=[0,127,47,59,79,21,123,63],stepsize=[0.001,0.001],nstep=5000,initial_guess=[-0.35,0.055],/gridsearch,/skipmapmaking,checkfolder='checksdss')
      ; stop
      feh_dev_par_cl = linfitprobgrid(ageform,sci.feh-fehideal,-1.*probsci.dage,probsci.dfeh,transpose(probsci.probdfehdage,[1,0,2]),sigma_a_b_cl,outfitfile='cl0024_prob_deltafeh_age.fits',outepsfile='cl0024_plot_deltafeh_age.eps',fitrange=[0.,8.],/plot,badid=[22,53,42,46],stepsize=[0.001,0.001],nstep=3000,initial_guess=[-0.35,0.055],/gridsearch,checkfolder='checkcl0024',/skipmapmaking)
      stop
     feh_Dev_par = linfitprobgrid(x,y,dxarr,dyarr,probdxdyarr,sigma_a_b,outfitfile='all_prob_deltafeh_age.fits',outepsfile='all_plot_deltafeh_age.eps',fitrange=[0,11.5],/plot,badid=[59,79,173,201],stepsize=[0.0005,0.0005],nstep=5000)
      save, sigma_a_b,sigma_a_b_sdss,sigma_a_b_cl,filename='linfitprobgrid_param_linfitprobgrid2.sav'
     stop
   endif else begin
      restore,'180130linfitprobgrid_param_linfitprobgrid2.sav'
      feh_dev_par_sdss = sigma_a_b_sdss[1,*]
      feh_dev_par_cl = sigma_a_b_cl[1,*]
      feh_dev_par = sigma_a_b[1,*]
   endelse
   print,'feh deviation (intercept, slope)'
   print,'combined',sigma_a_b
   print,'sdss:',sigma_a_b_sdss
   print,'cl:',sigma_a_b_cl

   xfitrange = [0,13]
   yfitrange = [-0.4,0.4]

   reduceimage,'all_prob_deltafeh_age.fits',xfitrange,yfitrange,img,ximg,yimg
   reduceimage,'cl0024_prob_deltafeh_age.fits',xfitrange,yfitrange,img_cl,ximg_cl,yimg_cl
   reduceimage,'sdss_prob_deltafeh_age.fits',xfitrange,yfitrange,img_sdss,ximg_sdss,yimg_sdss

   badimgcl = where(img_cl lt -0.3)
   img_cl(badimgcl) = -0.4
   badimgsdss = where(img_sdss lt -0.)
   img_sdss(badimgsdss) = -0.1


   plotposition=[0.15,0.15,0.95,0.9]

   set_plot,'x'
   loadct,60
   cgdisplay,/pixmap
   cgimage,img_sdss,position=plotposition
   bgimage = cgsnapshot()
   wdelete

   set_plot,'ps'
   device, filename = psname,xsize = 15,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;
      loadct,54
      outval = 99.
      good = where(finite(img),cgood,complement=bad)
      img(bad) = min(img(good))
      loadct,55
      cgimage,img_cl,transparent=50,alphabackgroundimage=bgimage,alphafgposition=plotposition
      plot,sdss_ageform,sdss.feh-sdss_fehideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle=delta+'[Fe/H]',xrange=xfitrange,xstyle=9,yrange=yfitrange,/nodata,position=plotposition,/noerase

      axis,xaxis=1,xticks=4,xtickv=[1.513,3.223, 5.747,8.422 ,12.161],xtickn=['4','2','1','0.5','0.1'],xtitle='z'

      cgplot,sdss_ageform,sdss.feh-sdss_fehideal,psym=16,/overplot,symsize=0.5,color='ygb5'
      cgplot,ageform,sci.feh-fehideal,psym=14,/overplot,color='org4'
      oplot,[0,14],feh_dev_par(0)+feh_dev_par(1)*[0,14],thick=2
      ;oplot,[0,14],feh_dev_par_sdss(0)+feh_dev_par_sdss(1)*[0,14],thick=2,linestyle=2,color=fsc_color('blu5')
      ;oplot,[0,14],feh_dev_par_cl(0)+feh_dev_par_cl(1)*[0,14],thick=2,linestyle=2,color=fsc_color('org6')

      xyouts,0.5,0.32,'(c) real observations',charsize=1.2
      xyouts,1,-0.38,'older galaxies',charsize=0.8
      xyouts,9.4,-0.38,'younger galaxies',charsize=0.8
      arrow,0.9,-0.37,0.3,-0.37,/data
      arrow,12.1,-0.37,12.7,-0.37,/data
   device,/close

   mkplottest,sci,sdss,ageform,sdss_ageform,fehideal,sdss_fehideal

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;Write table for latex
   make_catalog,sci.ra,sci.dec,sci.logmstar,sci.feh,sci.fehupper-sci.feh,sci.fehlower-sci.feh,$
                sci.age,sci.ageupper-sci.age,sci.agelower-sci.age,sci.snfit,sci.objname
   stop	
end

