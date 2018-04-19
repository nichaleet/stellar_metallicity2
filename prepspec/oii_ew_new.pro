function oii_ew_new,contdiv,lambda,contdivivar,z,widtherr=widtherr,znew=znew
   fixz=0
   znew=0.
   lambdain = lambda
   zinput = z

   linewaves = [2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85]
   linenames = ['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', '[NeIII]', 'HeI', 'H8', 'CaH', '[NeIII]', 'CaK', 'He', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', 'Mgb', 'Mgb', 'Mgb', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]']
   linecolors = ['blue', 'white', 'blue', 'white', 'white', 'white', 'white', 'blue', 'blue', 'white', 'red', 'blue', 'red', 'white', 'white', 'red', 'white', 'white', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'white', 'blue', 'blue', 'blue', 'blue']

   while fixz eq 0 do begin
      lamrange=[3720.,3735.]
      midlam = mean(lamrange)
      lambda  =lambdain/(z+1.)
      dlambda = lambda-shift(lambda,1)
      inrange = where(lambda gt lamrange[0] and lambda le lamrange[1] and finite(contdiv) and finite(contdivivar),cinrange)
  
      width = total(dlambda(inrange)*(1.-contdiv(inrange)))
      widtherr = sqrt(total((dlambda(inrange))^2/contdiv(inrange)))
      ;plotting
      window, 0, xsize=1200,ysize=800
      !p.multi = [0,1,2]
      ;first plot
      rangeplot=[3700,3750]
      plot,lambda,contdiv,xrange=rangeplot,/nodata
      if width gt 0 then factor = 0 else factor = 1
      polyfill,[midlam-width/2.,midlam+width/2.,midlam+width/2.,midlam-width/2.],[0,0,1,1]+factor,color=fsc_color('yellow')
      oplot,lambda,contdiv
      oplot,[3726,3726],!y.crange,color=fsc_Color('red'),linestyle=1
      oplot,[3729,3729],!y.crange,color=fsc_Color('red'),linestyle=1
      oplot,[lamrange[0],lamrange[0]],!y.crange,color=fsc_color('green'),linestyle=2
      oplot,[lamrange[1],lamrange[1]],!y.crange,color=fsc_color('green'),linestyle=2
      xyouts, 0.2,0.3,'EW= '+sigfig(width,3)+'+/-'+sigfig(widtherr,3),/normal
      ;second plot
      rangeplot=[3500,5300]
      inrange = where(lambda gt rangeplot[0] and lambda lt rangeplot[1])
      yrange = [min(contdiv(inrange),/nan)>(-2),max(contdiv(inrange),/nan)<(10)]
      plot,lambda,contdiv,xrange=rangeplot,yrange=yrange
      n = n_elements(linewaves)
      for i=0,n-1 do begin
         if (linewaves)[i] le !X.CRANGE[0] or (linewaves)[i] ge !X.CRANGE[1] then continue
         oplot, [(linewaves)[i], (linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((linecolors)[i])
         xyouts, (linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (linenames)[i], orientation=90, alignment=1, color=fsc_color((linecolors)[i])
      endfor
      ;read terminal option
      strfixz=''
      read,strfixz,prompt='fix redshift? (y/n)'
      if strfixz eq 'y' or strfixz eq 'o' then fixz=0 else fixz=1
      if strfixz eq 'y' then begin
      	read,newoiiwl,prompt='type in the location of oii that you see: '
      	z = (newoiiwl*(z+1)/3728.)-1.
      endif
      if strfixz eq 'o' then begin
        z= zinput
      endif
      znew = z
   endwhile
			
   return,width
end
