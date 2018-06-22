pro responsefn_ana
   respn = mrdfits('atlas_ssp_abund_krp.fits',1)
   wsel = where(respn.agegyr eq 3. and respn.zmet eq 0.,csel)
   str = respn(wsel)

   lambda = str.lam
   nl = n_elements(lambda)

   ;do mg/fe
   pspec = str.mgp/str.solar
   mspec = str.mgm/str.solar

   set_plot,'ps'
   !p.font=0
   epsname = 'responsefn_mg.eps'
   device, filename = epsname,xsize = 15,ysize = 30, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      nframe = 5
      pixfr = fix(nl/nframe)
      heightfr = (0.95-0.05)/nframe
      for i=0,nframe-1 do begin
          ipix = i*pixfr
          fpix = (i+1)*pixfr < nl
          loc = [0.1,heightfr*(nframe-i-1)+0.05,0.95,heightfr*(nframe-i)]
          yrange = minmax([mspec[ipix:fpix],pspec[ipix:fpix]])+[-0.1,0.1]
          plot,lambda[ipix:fpix],mspec[ipix:fpix],/nodata,position=loc,/noerase,$
               xrange=[lambda(ipix),lambda(fpix)],xstyle=1,yrange=yrange
          oplot,lambda[ipix:fpix],pspec[ipix:fpix],color=fsc_color('seagreen')
          oplot,lambda[ipix:fpix],mspec[ipix:fpix],color=fsc_color('lightseagreen')
      endfor
   device,/close

   ;alpha/Fe vs ca/fe vs mg/fe (plus only)
   mgspec = str.mgp/str.solar
   aspec = str.aFep/str.solar
   caspec = str.cap/str.solar
   epsname = 'responsefn_mg_ca_alpha.eps' 
   goodpix = where(lambda gt 3500. and lambda lt 6000.)
   device,filename = epsname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
        plot,lambda(goodpix),mgspec(goodpix),/nodata,yrange=[0.8,1.1],xtitle='lambda'
        oplot,lambda(goodpix),mgspec(goodpix),color=fsc_color('seagreen')
        oplot,lambda(goodpix),caspec(goodpix),color=fsc_color('darkgray')
   device,/close

   stop
end
