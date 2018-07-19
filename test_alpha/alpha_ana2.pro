pro alpha_ana2,redo=redo
   ;this code is to specifically track why there is an offset in alpha measurements
   ;that does not depend on signal to noise

   common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize, rest
   common get_sps_alpha, element

   ;read data
   alphafit = mrdfits('/scr2/nichal/workspace4/test_alpha/datasps/mockdata_alpha/sps_fit01.fits.gz',1)
   ;below is the first returned parameters when nloop = 0 (fit contdiv spec)
   alphafit[15].feh = -0.17818448
   alphafit[15].age = 2.8420257
   alphafit[15].vdisp = 103.12794
   alphafit[15].zfit = 0.55002154
   alphafit[15].alphafe = -0.19633903

   nwalkers = 5

   str = alphafit
   name='alphafit_5para'
   element= ['Mg','O','Si','Ca','Ti']

   ;get correct/input values
   nstr = n_elements(str)
   izmet = fltarr(nstr) ;i is for initial
   iage = fltarr(nstr)
   isn = fltarr(nstr)
   ialpha = fltarr(nstr)
   for i=0,nstr-1 do begin
       exts = strsplit(str[i].objname,'_',/extract)
       mop = (strmid(exts[0],4,1) eq 'm') ? -1. : 1.
       izmet[i] = mop*(float(fix(strmid(exts[0],5,1)))+0.01d*float(fix(strmid(exts[0],7,2))))
       iage[i] = float(fix(strmid(exts[1],0,1)))
       mop = (strmid(exts[2],5,1) eq 'm') ? -1. : 1
       ialpha[i] = mop*(float(fix(strmid(exts[2],6,1)))+0.01d*float(fix(strmid(exts[2],8,2))))
       isn[i] = float(fix(strmid(exts(3),2)))
   endfor

   sel = where(izmet eq -0.2 and isn eq 30,csel) 
   str = str(Sel)
   izmet = izmet(sel)
   isn = isn(sel)
   iage = iage(sel)
   ialpha = ialpha(sel)
   if keyword_Set(redo) then begin
      chibestval = fltarr(csel)
      chirealval = fltarr(csel) 
      chibestval_contdiv = fltarr(csel)
      chirealval_contdiv = fltarr(csel)

      ;compare chisq of the best fit and of the truth
      for i=0,csel-1 do begin
         science = str[i]
         znow = science.zspec
         reallambda = science.lambda
         won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
         dlam = science.dlam/2.35
         xmp = reallambda[won]
         ymp = science.contdiv[won]/science.spscont[won]
         ymp2 = science.contdiv[won]
         wonfit = won
         dataivar = science.contdivivar
         datalam = science.lambda
         contmask = science.contmask
         rest =0
         normalize = 0
       
         spsrealval = get_sps_alpha_obs(xmp,[izmet[i],iage[i],300./2.35,0.55,ialpha[i]])
         spsbestval = get_sps_alpha_obs(xmp, [science.feh,science.age,science.vdisp,science.zfit,science.alphafe])
       
         chibestval[i] = total((spsbestval-ymp)^2*science.contdivivar[won]*(science.spscont[won])^2)
         chirealval[i] = total((spsrealval-ymp)^2*science.contdivivar[won]*(science.spscont[won])^2)

         normalize = 1
         spsrealval = get_sps_alpha_obs(xmp,[izmet[i],iage[i],300./2.35,0.55,ialpha[i]])
         spsbestval = get_sps_alpha_obs(xmp, [science.feh,science.age,science.vdisp,science.zfit,science.alphafe])
         chibestval_contdiv[i] = total((spsbestval-ymp2)^2*science.contdivivar[won])
         chirealval_contdiv[i] = total((spsrealval-ymp2)^2*science.contdivivar[won])
      endfor
       save,chibestval,chirealval,chibestval_contdiv,chirealval_contdiv,filename='savedpara_alphaana2.sav'
   endif else restore,'savedpara_alphaana2.sav'
   ;calculate differences
   dage = alog10(str.age/iage)
   dageu = alog10(str.ageupper/iage)
   dagel = alog10(str.agelower/iage)
   dz  = str.feh-izmet 
   dzu = str.fehupper-izmet
   dzl = str.fehlower-izmet
   dalpha = str.alphafe-ialpha
   dalphau = str.alphafeupper-ialpha
   dalphal = str.alphafelower-ialpha
   dchi = chibestval-chirealval
   dchi_contdiv = chibestval_contdiv-chirealval_contdiv
   ;symbols
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
   chisq = '!9'+string("143B)+'!x!E2!N'
   Delta = '!9'+string("104B)+'!x'

   ;plotting
   set_plot,'ps'
   !p.font=0
   psname='plots/realval_bestfitval_compare.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,1,1]
      plot,ialpha,dchi,psym=cgsymcat(15),xtitle='['+alpha+'/zmet]',ytitle=delta+chisq+'(best-real)'
      oplot,[ialpha[0]],[dchi[0]],psym=cgsymcat(15),color=fsc_color('indianred')
   device,/close

   psname='plots/realval_bestfitval_contdiv_compare.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi = [0,1,1]
      plot,ialpha,dchi_contdiv,psym=cgsymcat(15),xtitle='['+alpha+'/zmet]',ytitle=delta+chisq+'(best-real)'
      oplot,[ialpha[0]],[dchi_contdiv[0]],psym=cgsymcat(15),color=fsc_color('indianred')

   device,/close
stop
end 
