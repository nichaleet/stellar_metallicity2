pro reduced_cl1604ana
   z = 0.96
   file = file_search('/scr2/nichal/WMKO/LRIS/2017may_cl1604/reduce/cl1604_1.r.*.msdc.fits',count=cfile)
   for i=0,cfile-1 do begin
      spec = readfits(file[i],hdr)
      npix = n_elements(spec)
      lambda = findgen(npix)*sxpar(hdr,'CDELT1')+sxpar(hdr,'CRVAL1')
      plot,lambda/(1.+z),spec,yrange=[0,3*median(spec)]
      
      ;signal to noise
      contpara = poly_fit(lambda,spec,6,yfit=cont)
      oplot,lambda/(1.+z),cont,color=fsc_color('red')

      ;calculate signal to noise
      dev = abs((spec-cont)/cont)
      avgdev = mean(dev)
      w = where(dev lt 3.0*avgdev, c)
      if c gt 0 then sn_perpix = 1.0/median(dev[w])
      print,file[i]
      print,'SN:',sn_perpix

      stop
   endfor
end
