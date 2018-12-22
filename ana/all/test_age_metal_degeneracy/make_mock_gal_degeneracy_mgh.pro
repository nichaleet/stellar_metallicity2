pro make_mock_gal_degeneracy_mgh,sdss=sdss,cl0024=cl0024,ms0451=ms0451
;This program is a nest inside test_degeneracy_formation_time2.pro
;This is the part where it makes a fits files of mock galaxies whose the real MZR follows
;a tight MZR but scattered by age-metal degeneracy
;The SN,age,mass are taken directly from SDSS data
common mar_lin_param, marparam
common runname, name

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;setting
	if keyword_set(sdss) then begin 
		sci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_sdss_ana_afterev_final.fits',1)
		marparam = marparam_sdss
		outfile = name+'_gal_degeneracy.fits'
                proboutfile = name+'sdss_feh_age_probdist.fits'
		type = 'SDSS'
	endif
	if keyword_Set(cl0024) then begin
		sci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_cl0024_ana_afterrev_final.fits',1)
		marparam = marparam_cl
		outfile = name+'_cl0024_degeneracy.fits'
                proboutfile = name+'cl0024_feh_age_probdist.fits'
		type = 'Cl0024'
	endif

	ngals = n_elements(sci)

	nfeh = 120
	nage = 126
	grid_Feh = rebin(findgen(nfeh)/120.-0.8,nfeh,nage) ;range from [-0.8,0.192] dex ;
	grid_Age = transpose(rebin(findgen(126)/10.+0.5 ,nage,nfeh)) ;range from [0.5,13] Gyr
	deltafeh = grid_feh[1,0]-grid_feh[0,0]
        deltaage = grid_age[0,1]-grid_age[0,0]
	feharr = grid_feh[*,0]
	agearr = grid_age[0,*]

	;output structures
	indi = {sn:0.,mass:0.,fehideal:0.,ageideal:0.,feh:0.,dfehupper:0.,dfehlower:0.,age:0.,dageupper:0.,dagelower:0.,z:0.,vdisp:0.}
	str  = replicate(indi,ngals)

        probfeh_age_str = {objname:'name',dfeh:feharr,probdfeh:dblarr(nfeh),dage:agearr,probdage:dblarr(nage),probdfehdage:dblarr(nfeh,nage),dageupper:0.,dagelower:0.,dfehupper:0.,dfehlower:0.}
        probfeh_age = replicate(probfeh_age_str,ngals)
	
	;assign mass as the observed mass
	str.mass=sci.logmstar
	;assign ideal metallicities
	str.fehideal = marparam[0]+(str.mass-10.)*marparam[1]
	;get ideal age as the observed age
	str.ageideal = sci.age
	str.z = sci.zfit
	str.sn = sci.sn
	str.vdisp = sci.vdisp

        ;assign things in probfeh_age
        for i=0,ngals-1 do begin
           probfeh_age[i].objname = strtrim(string(i),2)
           probfeh_age[i].dfeh = feharr-str[i].fehideal
           probfeh_age[i].dage = agearr-str[i].ageideal
        end

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;run a loop for each 'galaxy' to get the 'observed feh and observed age'
      ;  if type eq 'SDSS' then begin
      ;      istart = 80
      ;      str = mrdfits('temp_degeneracy6.fits',9)
      ;  endif else istart = 0
        istart = 0
	for i=istart,ngals-1 do begin
		print,type+' now doing galaxy number',i
		;create the synthesis spec with template from sci
		if keyword_set(sdss) then spec_arr = make_spec_sdss_specify_template(grid_feh,grid_age,sci(i))		
		if keyword_set(cl0024) then spec_arr = make_spec_specify_template(grid_feh,grid_age,sci(i))

		;find the matched spectrum
		dist = abs(spec_arr.feh-str(i).fehideal)+abs(spec_arr.age-str(i).ageideal)
		loc = where(dist eq min(dist))
		print,strtrim(string(fix(i)),2)+' FeH=',str(i).fehideal,' use ',spec_arr(loc).feh
		print,strtrim(string(fix(i)),2)+' AGE=',str(i).ageideal,' use ',spec_arr(loc).age

		;create noised spectrum
		;use the noise as a fn of wl from observed spectra
		refspecstr = spec_arr(loc)
		won = where(sci(i).fitmask eq 1 and finite(sci(i).contdiv) and finite(sci(i).contdivivar) $
			and sci(i).contdivivar gt 0 and sci(i).lambda/(1.+sci(i).zfit) gt 3500. $
			and sci(i).lambda/(1.+sci(i).zfit) lt 7400., con)
		if n_Elements(refspecstr.spec) ne con then stop,'oh oh'
		refspec_err = abs(refspecstr.spec/(sci(i).contdiv(won)*Sqrt(sci(i).contdivivar(won))))
		npix = con
		refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

		;make chisq arr
		chisqarr = fltarr(nfeh,nage)           
		for iz=0,nfeh-1 do for ia=0,nage-1 do begin
			if grid_feh[iz,ia] ne spec_arr[iz,ia].feh or grid_age[iz,ia] ne spec_arr[iz,ia].age $
				then stop,'emm sth is wront'
			chisqarr[iz,ia] = total((refspec-spec_arr[iz,ia].spec)^2/refspec_err)
		endfor

		Lgrid = -0.5*(chisqarr-min(chisqarr))
		probgrid = exp(double(Lgrid))
		volume = total(probgrid)*deltafeh*deltaage
		probgrid = probgrid/volume

                ;SAVING THE IDEAL PROB DISTRIBUTION
                probfeh_age(i).probdfehdage = probgrid

                probfeh = fltarr(nfeh)
                probage = fltarr(nage)
                for ii=0,nfeh-1 do probfeh(ii) = int_tabulated(agearr,probgrid[ii,*])
                for jj=0,nage-1 do probage(jj) = int_tabulated(feharr,probgrid[*,jj])
                probfeh = probfeh/int_tabulated(feharr,probfeh) ;make sure total prob is 1
                probage = probage/int_tabulated(agearr,probage) ;make sure total prob is 1
                probfeh_age(i).probdfeh = probfeh
                probfeh_age(i).probdage = probage

                cumprobfeh = fltarr(nfeh)
                for jj=1,nfeh-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
                cumprobage = fltarr(nage)
                for jj=1,nage-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])

                probfeh_age(i).dageupper = interpol(agearr,cumprobage,0.84)-str(i).ageideal
                probfeh_age(i).dagelower = interpol(agearr,cumprobage,0.16)-str(i).ageideal
                probfeh_age(i).dfehupper = interpol(feharr,cumprobfeh,0.84)-str(i).fehideal
                probfeh_age(i).dfehlower = interpol(feharr,cumprobfeh,0.16)-str(i).fehideal
                ;;;;;;;;;;;;FINISH SAVING IDEAL PROB DISTRIBUTION

		;select the points according to prob
		maxprob = max(probgrid)
		wgoodchi = where(probgrid gt 0. and probgrid gt maxprob*0.05,cgoodchi) ;should be all non zero prob and within 2.5 sigma
		prob = probgrid(wgoodchi)
		cumuprob =total(prob,/cumulative)*deltafeh*deltaage  ;should be zero to close to 1
		randomnum = randomu(seed)*cumuprob(cgoodchi-1)
		chosen = value_locate(cumuprob,randomnum)+1
		str(i).feh = grid_feh(wgoodchi(chosen))
		str(i).age = grid_age(wgoodchi(chosen))

		;now get the uncertainties of these detected point
		;gotta make a new noised spec with these 'observed' values
                dist = abs(spec_arr.feh-str(i).feh)+abs(spec_arr.age-str(i).age)
                loc = where(dist eq min(dist))
		print,'observed values'
                print,strtrim(string(fix(i)),2)+' FeH=',str(i).feh,' use ',spec_arr(loc).feh
                print,strtrim(string(fix(i)),2)+' AGE=',str(i).age,' use ',spec_arr(loc).age

                ;create noised spectrum
                ;use the noise as a fn of wl from observed spectra
                refspecstr = spec_arr(loc)
                if n_Elements(refspecstr.spec) ne con then stop,'oh oh'
                refspec_err = abs(refspecstr.spec/(sci(i).contdiv(won)*Sqrt(sci(i).contdivivar(won))))
                refspec = refspecstr.spec+randomn(seed,npix)*refspec_err
 		;make chisq arr
                chisqarr = fltarr(nfeh,nage)
                for iz=0,nfeh-1 do for ia=0,nage-1 do begin
                        if grid_feh[iz,ia] ne spec_arr[iz,ia].feh or grid_age[iz,ia] ne spec_arr[iz,ia].age $
                                then stop,'emm sth is wront'
                        chisqarr[iz,ia] = total((refspec-spec_arr[iz,ia].spec)^2/refspec_err)
                endfor
                Lgrid = -0.5*(chisqarr-min(chisqarr))
                probgrid = exp(double(Lgrid))
                volume = total(probgrid)*deltafeh*deltaage
                probgrid = probgrid/volume

		probfeh = fltarr(nfeh)
		probage = fltarr(nage)
		for ii=0,nfeh-1 do probfeh(ii) = int_tabulated(agearr,probgrid[ii,*])
		for jj=0,nage-1 do probage(jj) = int_tabulated(feharr,probgrid[*,jj])
		probfeh = probfeh/int_tabulated(feharr,probfeh)
		probage = probage/int_tabulated(agearr,probage)

		cumprobfeh = fltarr(nfeh)
		for jj=1,nfeh-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
		cumprobage = fltarr(nage)
		for jj=1,nage-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])

		midagevalue = interpol(agearr,cumprobage,0.5)
		midfehvalue = interpol(feharr,cumprobfeh,0.5)

		str(i).dageupper = interpol(agearr,cumprobage,0.84)-midagevalue
		str(i).dagelower = interpol(agearr,cumprobage,0.16)-midagevalue
		str(i).dfehupper = interpol(feharr,cumprobfeh,0.84)-midfehvalue
		str(i).dfehlower = interpol(feharr,cumprobfeh,0.16)-midfehvalue
		if i mod 10 eq 0 then mwrfits,str,'temp_degeneracy6.fits',/silent,/create
	endfor

	mwrfits,str,outfile,/silent,/create
        mwrfits,probfeh_age,proboutfile,/silent,/create       
end


