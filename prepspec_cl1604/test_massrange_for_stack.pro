pro test_massrange_for_stack,massrange=massrange

       if ~keyword_Set(massrange) then massrange=[10.,10.4,10.6,10.75,11.05,11.5]
       inputfile = '/scr2/nichal/workspace4/sps_fit/data/spline_cl1604/sps_fit01.fits.gz'
       data = mrdfits(inputfile,1)
       ;calculate sn per angstrom
       snperpix = data.sn
       snperang = snperpix
       for i=0,n_elements(snperpix)-1 do begin
          dlam = median(-1.*ts_diff(data[i].lambda,1)/(1.+data[i].zcat)) ;rest lambda
          snperang[i] = snperpix[i]/sqrt(dlam)
       ;   print, dlam
       endfor
       ;get data for each mass range
       print, 'mass range,  ngal, est. final SN'
       for i=0,n_elements(massrange)-2 do begin
          sel = where(data.logmstar_sed_bl ge massrange[i] and data.logmstar_sed_bl lt massrange[i+1] $
                      and snperang gt 5 and snperang le 15 and data.oiiew gt -5.,ndup)
          if ndup lt 3 then stop,'Warning: there are less than 3 galaxies in this mass range'
          print, [massrange[i:i+1],ndup, sqrt(total((snperang(sel))^2))], format='(f6.2,":",f6.2,3x,I3,3x,f6.2)'
       endfor
     ; stop
end
