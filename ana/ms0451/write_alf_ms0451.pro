pro write_alf_ms0451_original

   sci = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
   hisn = max(sci.snfit,sel) 
   sci = sci(sel)
   lambda = sci.lambda/(1.+sci.zfit)
   flux = sci.contdiv
   error = 1./sqrt(sci.contdivivar)
   weight = sci.fitmask
   res = sci.dlam/sci.lambda*3.e5 ;km/s
 
   dir = '/scr2/nichal/alf/indata/'
   outname = dir+sci.objname+'_input.dat'
   print, outname+'   written'
   writecol,outname,lambda,flux,error,weight,res
    
end

pro write_alf_ms0451
   sciall = mrdfits('/scr2/nichal/workspace4/sps_fit/data/spline_ms0451/sps_fit03.fits.gz',1)
   goodsn = where(sciall.snfit gt 26 and sciall.snfit lt 30.,csel)
   for i=0,csel-1 do begin
      sel = goodsn(i)
      sci = sciall(sel)

      lambda = sci.lambda/(1.+sci.zfit)
      flux = sci.contdiv
      error = 1./sqrt(sci.contdivivar)
      weight = sci.fitmask
      res = sci.dlam/sci.lambda*3.e5 ;km/s
      goodwl = where(lambda gt 3700. and lambda lt 6000.)
      dir = '/scr2/nichal/alf/indata/'
      outname = dir+sci.objname+'_input.dat'
      print, outname+'   written'
      writecol,outname,lambda(goodwl),flux(goodwl),error(goodwl),weight(goodwl),res(goodwl)
  endfor
end

