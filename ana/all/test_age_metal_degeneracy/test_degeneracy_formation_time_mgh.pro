pro test_degeneracy_formation_time_mgh
common mzr_lin_param, marparam
common runname, name
;this program is to test whether age-degeneracy can cause the scatter in the observed Mass-Mg Relation
;and create a separation into formation time.
;It first assumes that galaxies follow a tight correlation of MZR from the best-linear fit
;Then galaxies adopt random values within chisq=chisq+min(chisq)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   name='mgh'
   marparam = [0.09,0.10]
   marr = findgen(101)/40.+9. ;logmass from 9 to 11.5
;open or make the mock set of galaxies
   if file_test(name+'_gal_degeneracy.fits') eq 0 then make_mock_gal_degeneracy_mgh,/sdss
   sdss = mrdfits(name+'_gal_degeneracy.fits',1)
   if file_test(name+'_cl0024_degeneracy.fits') eq 0 then make_mock_gal_degeneracy_mgh,/cl0024
   sci = mrdfits(name+'_cl0024_degeneracy.fits',1)
   ;= {sn:0.,mass:0.,fehideal:0.,ageideal:0.,feh:0.,dfehupper,dfehlower,dageupper,dagelower.,age:0.,ageerr:0.}

;   redouncertainty=0
;   if redouncertainty then begin
;      namesci = strtrim(string(indgen(n_Elements(sci))),2)
;      namesdss = strtrim(string(indgen(n_Elements(sdss))),2)
;      get_additional_uncertainties_mockcl0024,sci.fehideal,sci.ageideal,sci.sn,namesci,outname='chisqcube_'+name+'cl0024_sci.fits'
;      get_prob_dist,filein='chisqcube_'+name+'cl0024_sci.fits',fileout=name+'cl0024_feh_age_probdist.fits'
;      get_additional_uncertainties_mocksdss,sdss.fehideal,sdss.ageideal,sdss.sn,namesdss,outname='chisqcube_'+name+'sdss_sci.fits'
;      get_prob_dist,filein='chisqcube_'+name+'sdss_sci.fits',fileout=name+'sdss_feh_age_probdist.fits'
;   endif
   probsci = mrdfits(name+'cl0024_feh_age_probdist.fits',1)
   probsdss = mrdfits(name+'sdss_feh_age_probdist.fits',1)
 
   sdss_ageform = (galage(0.05,1000.)/1.e9-sdss.ageideal)>0.
   ageform = (galage(sci.z,1000.)/1.e9-sci.ageideal)>0.

   redolinfit=0
   if redolinfit then begin
 ;     feh_dev_par_sdss1 = linfitprobgridmock(sdss_ageform,sdss.fehideal-sdss.fehideal,-1.*probsdss.dage,probsdss.dfeh,$
 ;                        transpose(probsdss.probdfehdage,[1,0,2]),sigma_a_b_sdss1,$
 ;                        outfitfile=name+'sdss_prob_deltafeh_age.fits',$
 ;                        outepsfile=name+'sdss_plot_deltafeh_age.eps',fitrange=[0,11.5])
      feh_dev_par_sdss = linfitprobgridmock2(sdss_ageform,sdss.fehideal-sdss.fehideal,-1.*probsdss.dage,probsdss.dfeh,$
                         transpose(probsdss.probdfehdage,[1,0,2]),sigma_a_b_sdss,$
                         outfitfile=name+'sdss_prob_deltafeh_age_linfitprobgridmock2.fits',$
                         outepsfile=name+'sdss_plot_deltafeh_age_linfitprobgridmock2.eps',fitrange=[0,11.5],/plot)

      feh_dev_par_cl= linfitprobgridmock2(ageform,sci.fehideal-sci.fehideal,-1.*probsci.dage,probsci.dfeh,$
                         transpose(probsci.probdfehdage,[1,0,2]),sigma_a_b_cl,$
                         outfitfile=name+'cl0024_prob_deltafeh_age_linfitprobgridmock2.fits',$
                         outepsfile=name+'cl0024_plot_deltafeh_age_linfitprobgridmock2.eps',fitrange=[0,8],/plot)
  ;    feh_dev_par_cl1= linfitprobgridmock(ageform,sci.fehideal-sci.fehideal,-1.*probsci.dage,probsci.dfeh,$
  ;                       transpose(probsci.probdfehdage,[1,0,2]),sigma_a_b_cl1,$
  ;                       outfitfile=name+'cl0024_prob_deltafeh_age.fits',$
  ;                       outepsfile=name+'cl0024_plot_deltafeh_age.eps',fitrange=[0,8])

      x=[sdss_ageform,ageform]
      y=[sdss.fehideal-sdss.fehideal,sci.fehideal-sci.fehideal]
      dxarr = -1.*[[reform(probsdss.dage)],[reform(probsci.dage)]] ;since x axis is ageuni-agegal
      dyarr = [[probsdss.dfeh],[probsci.dfeh]]
      probdxdyarr = transpose([[[probsdss.probdfehdage]],[[probsci.probdfehdage]]],[1,0,2])

      feh_Dev_par = linfitprobgridmock2(x,y,dxarr,dyarr,probdxdyarr,sigma_a_b,$
                    outfitfile=name+'all_prob_deltafeh_age_linfitprobgridmock2.fits',$
                    outepsfile=name+'all_plot_deltafeh_age_linfitprobgridmock2.eps',$
                    fitrange=[0,11.5],/plot)
   ;   feh_Dev_par1 = linfitprobgridmock(x,y,dxarr,dyarr,probdxdyarr,sigma_a_b1,$
   ;                 outfitfile=name+'all_prob_deltafeh_age.fits',$
   ;                 outepsfile=name+'all_plot_deltafeh_age.eps',$
   ;                 fitrange=[0,11.5])

    ;  save, sigma_a_b1,sigma_a_b_sdss1,sigma_a_b_cl1,filename='linfitprobgrid_'+name+'param_coadd.sav'
      save, sigma_a_b,sigma_a_b_sdss,sigma_a_b_cl,filename='linfitprobgrid_'+name+'.sav'
   endif 
      restore,'linfitprobgrid_'+name+'.sav'
      feh_dev_par_sdss = sigma_a_b_sdss[1,*]
      feh_dev_par_cl = sigma_a_b_cl[1,*]
      feh_dev_par = feh_dev_par_cl
  

   sdss_ageform = (galage(0.05,1000.)/1.e9-sdss.age)>0.
   ageform = (galage(sci.z,1000.)/1.e9-sci.age)>0.
	
   set_plot,'ps'
   !p.font =  0 
   sunsym= sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plotting
   psname='Formation_Redshift_'+name+'.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      zcat = [1000,2,1,0.7,0.4,0.]
      agecat = galage(zcat,1000.)/1.e9 ;age of universe at that redshif	
      plot,sdss.mass,sdss.feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],xstyle=1,yrange=[-0.8,0.2],/nodata
      
      nage = n_Elements(agecat)-1
      rainbow_colors
      zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
      ;loop over formation times
      for i=0,nage-1 do begin
         selsdss = where(sdss_ageform gt agecat(i) and sdss_ageform le agecat(i+1),cselsdss)
         if cselsdss gt 0 then begin
         	cgplot,sdss(selsdss).mass,sdss(selsdss).feh,psym=16,/overplot,color=zcolor(i),symsize=0.7
         endif
      endfor
      
      for i=0,nage-1 do begin
         sel = where(ageform gt agecat(i) and ageform le agecat(i+1), csel)
         if csel gt 0 then begin
            cgerrplot,sci(sel).mass,sci(sel).dfehlower+sci(sel).feh,sci(sel).dfehupper+sci(sel).feh,color=zcolor(i)-10,thick=0.5
            cgplot,sci(sel).mass,sci(sel).feh,psym=14,/overplot,color=zcolor(i),symsize=1.3
            cgplot,sci(sel).mass,sci(sel).feh,psym=4,/overplot,color='darkgray',symsize=1.3
         endif
      endfor
      ;ADD THE BEST FITTED LINEAR FUNCTION
      oplot,marr,(marr-10.)*mzrparam_cl[1]+mzrparam_cl[0],linestyle=2
      oplot,marr,(marr-10.)*mzrparam_sdss[1]+mzrparam_sdss[0],linestyle=2
      
      ;Labelling
      zarr_str = strarr(n_elements(zcat)-1)
      for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z$\tex_{form}$<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
      zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(F3.1)'),2)
      al_Legend,zarr_str,psym=15,color=zcolor,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
      al_Legend,['mock SDSS','mock Cl0024 z~0.4'],psym=[16,14],symsize=[0.5,1.3],color=0,box=0,thick=2,charsize=1,position=[10.6,-0.3],font=0
      xyouts,9.1,0.1,'(b)',charsize=1.2
   device,/close

   psname='FEH_deviation_'+name+'.eps'
   Delta = '!9'+string("104B)+'!x'
   Delta = '!9'+string("104B)+'!x'
   device, filename = psname,xsize = 15,ysize = 10, $
   		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   	plot,sdss_ageform,sdss.feh-sdss.fehideal,psym=1,xtitle='Age of Universe at Formation',ytitle=delta+'[Fe/H]',xrange=[0,14],xstyle=9,yrange=[-0.4,0.4],/nodata,position=[0.15,0.15,0.95,0.9]
   	axis,xaxis=0,xstyle=1,xrange=[0,14]
   	axis,xaxis=1,xticks=4,xtickv=[1.558,3.316, 5.903, 8.628,13.712],xtickn=['4','2','1','0.5','0'],xtitle='z'
   	cgplot,sdss_ageform,sdss.feh-sdss.fehideal,psym=16,/overplot,symsize=0.5,color='ygb5'
   	cgplot,ageform,sci.feh-sci.fehideal,psym=14,/overplot,color='org4'
   	x = [sdss_ageform,ageform]
   	y = [sdss.feh-sdss.fehideal,sci.feh-sci.fehideal]
   	dx = 0.5*[sdss.dageupper-sdss.dagelower,sci.dageupper-sci.dagelower]
   	dy = 0.5*[sdss.dfehupper-sdss.dfehlower,sci.dfehupper-sci.dfehlower]
   	fitexy,x,y,A,B,x_sig=dx,y_sig=dy,sigma_a_b,chi_sq,q	
   	print,A,B,sigma_a_b ;   -0.0781566   0.00968366     0.014947358    0.0022390772
   	oplot,[0,14],A+B*[0,14],color=cgcolor('black'),thick=2
   	cgerrplot,[12.],-0.25+[median([sdss.dfehlower,sci.dfehlower])],$
   		-0.25+[median([sdss.dfehupper,sci.dfehupper])]
   	cgerrplot,-0.25,12.+[median([sdss.dagelower,sci.dagelower])],$
   		12.+[median([sdss.dageupper,sci.dageupper])],/horizontal
   	
           xyouts,0.5,0.32,'(d) mock observations without evolution',charsize=1.2
   	xyouts,1,-0.38,'older galaxies',charsize=0.8
   	xyouts,10,-0.38,'younger galaxies',charsize=0.8
   	arrow,0.9,-0.37,0.3,-0.37,/data
   	arrow,13.1,-0.37,13.7,-0.37,/data
   
   device,/close		

   xfitrange = [0,13]
   yfitrange = [-0.4,0.4]
   reduceimage,name+'all_prob_deltafeh_age.fits',xfitrange,yfitrange,img,ximg,yimg
   reduceimage,name+'cl0024_prob_deltafeh_age.fits',xfitrange,yfitrange,img_cl,ximg_cl,yimg_cl
   reduceimage,name+'sdss_prob_deltafeh_age.fits',xfitrange,yfitrange,img_sdss,ximg_sdss,yimg_sdss

   imsize = size(img,/dimension)
   badimgcl = where(img_cl lt -0.1)
   img_cl(badimgcl) = -0.2
   badimgsdss = where(img_sdss lt 0.)
   img_sdss(badimgsdss) = -0.1

   plotposition=[0.15,0.15,0.95,0.9]
   set_plot,'x'
   loadct,60
   cgdisplay,/pixmap
   cgimage,img_sdss,position=plotposition
   bgimage = cgsnapshot()
   wdelete

   set_plot,'ps'

   psname='FEH_deviation_fromz0line_'+name+'_probdist.eps'
   device, filename = psname,xsize = 15,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;
      loadct,54
      outval = 99.
      good = where(finite(img_cl),cgood,complement=bad,ncomplement=cbad)
      if cbad gt 0 then img_cl(bad) = min(img_cl(good))
      loadct,55
      cgimage,img_cl,transparent=50,alphabackgroundimage=bgimage,alphafgposition=plotposition,missing_value=-99.
      plot,sdss_ageform,sdss.feh-sdss.fehideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle=delta+'[Fe/H]',xrange=xfitrange,xstyle=9,yrange=yfitrange,/nodata,position=plotposition,/noerase

      axis,xaxis=1,xticks=4,xtickv=[1.513,3.223, 5.747,8.422 ,12.161],xtickn=['4','2','1','0.5','0.1'],xtitle='z'

      cgplot,sdss_ageform,sdss.feh-sdss.fehideal,psym=16,/overplot,symsize=0.5,color='ygb5'
      cgplot,ageform,sci.feh-sci.fehideal,psym=14,/overplot,color='org4'
      oplot,[0,14],feh_dev_par(0)+feh_dev_par(1)*[0,14],thick=2
    ;  oplot,[0,14],feh_dev_par_sdss(0)+feh_dev_par_sdss(1)*[0,14],thick=2,linestyle=2,color=fsc_color('blu5')
    ;  oplot,[0,14],feh_dev_par_cl(0)+feh_dev_par_cl(1)*[0,14],thick=2,linestyle=2,color=fsc_color('org6')
 
;      oplot,[0,14],A+B*[0,14],color=cgcolor('black'),thick=2

      xyouts,0.5,0.32,'(d) mock observations without evolution',charsize=1.2
      xyouts,1,-0.38,'older galaxies',charsize=0.8
      xyouts,9.4,-0.38,'younger galaxies',charsize=0.8
      arrow,0.9,-0.37,0.3,-0.37,/data
      arrow,12.1,-0.37,12.7,-0.37,/data
   device,/close
   stop
end 



