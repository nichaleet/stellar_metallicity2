pro plot_sdss_spec

;http://skyserver.sdss.org/dr12/en/tools/quicklook/summary.aspx?id=1237648720693821757
;http://skyserver.sdss.org/dr12/en/tools/quicklook/summary.aspx?id=1237648720693821797

   sf = READ_CSV('sfgal_sdss.csv', HEADER=SfHeader, $
      N_TABLE_HEADER=1, TABLE_HEADER=sftableheader)
      ;field 1 = wl, field2 = flux, field3 = bestfit, field 4 = skyflux
   passive = READ_CSV('passivegal_sdss.csv', HEADER=passiveHeader, $
      N_TABLE_HEADER=1, TABLE_HEADER=passivetableheader)
   zsf = 0.1514758
   zpassive = 0.1747471

   linewaves = [2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85]
   linenames = ['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', ' ', ' ', 'H8', 'CaK', ' ', 'CaH', ' ', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', ' ', 'Mgb', ' ', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]']
   linecolors = ['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue']

   set_plot,'ps'
   !p.font = 0
   !p.charsize = 1.5
   psname = 'example_spec.eps'
   device, filename = psname,xsize = 25,ysize = 10,decomposed=1,color=1, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated;

      plot,sf.field1/(1.+zsf),sf.field3,xrange=[3700,7000],/nodata,xstyle=1,$
        ytitle='f!D!9'+String("154B)+'!X!N (10!U-17!N erg/s/cm!U2!N/Ang)',$
        xtitle = 'wavelength(Ang)'
      aa = string("154B)  
      oplot, passive.field1/(1.+zpassive),passive.field3,color=fsc_color('red')
      oplot,sf.field1/(1.+zsf),sf.field3,color=fsc_Color('blue')

      n = n_elements(linewaves)
      for j=0,n-1 do begin
      if (linewaves)[j] le !X.CRANGE[0] or (linewaves)[j] ge !X.CRANGE[1] then continue
      oplot, [(linewaves)[j], (linewaves)[j]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((linecolors)[j])
      xyouts, (linewaves)[j]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (linenames)[j], orientation=90, alignment=1, color=fsc_color((linecolors)[j]),charsize=0.8
   endfor
   device,/close
end
