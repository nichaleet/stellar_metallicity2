function linfitprobgrid_general,x,y,dxarr,dyarr,probdxdyarr,sigmaout,outfitfile=outfitfile,outepsfile=outepsfile,$
         xrange=xrange,yrange=yrange,nxout=nxout,nyout=nyout,fitrange=fitrange,plot=plot,badid=badid,$
         stepsize=stepsize,nstep=nstep,initial_guess=initial_guess,gridsearch=gridsearch,$
         skipmapmaking=skipmapmaking,checkfolder=checkfolder,p0range=p0range,p1range=p1range
;Likelihood maximization (multiply prob together and fit linear function)
;INPUT::
;       x     = array of x values of size nobjs
;       y     = array of y values of size nobjs
;       dxarr = array of dx (so x could be x+dx) of size ndx X nobjs
;       dyarr = array of dy (so y could be y+dy) of size ndy X nobjs
;       probdxdyarr = array of probabilities of each associated x+dxarr,y+dxarr. The size is ndx X ndy X nobjs
;       outfitfile = string of file name to create a fits file (image) of the final total probability (sum over all obj)
;       fitrange = range in x to fit linear line
;       FOR MCMC:
;       stepsize = step size in intercept, slope
;       nstep = one number, default = 3000 
;       initial_guess = first guess in intercept, slope
;OPTIONAL INPUT::
;       xrange,yrange = Each is a 2-element array. So,the input x and y don't have to be the same of each object. 
;                       The program will interpolate them to the same grid for all object. 
;                       xrange and yrange are this common grid range. Default is minmax of all objects. 
;       nxout,nyout = number of elements in each x and y dimension of the common grid above. Default is 50.
;       /gridsearch = set if you want to do grid search instead of mcmc. if so, you gotta specify p0range and 
;                     p1range for intercept and slope range accordingly. The resolution of grid is 25x25 
;       outepsfile = output best fit param plots
;OUTPUT:: sigmaout = [[intercept16,intercept50,intercept84],[slope16,slope50,slope84]] 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   dimen = size(probdxdyarr,/dimensions)
   ngridx = dimen(0)
   ngridy = dimen(1)
   nobjs = dimen(2)
 
   xarr = transpose(rebin(x,nobjs,ngridx))+dxarr ;dimension of ndx X nobjs
   yarr = transpose(rebin(y,nobjs,ngridy))+dyarr ;dimension of ndy X nobjs

   ;assign value of where prob is 0 = minprob/1.e10
   zeroval = 0. ;min(probdxdyarr(where(probdxdyarr ne 0)))/1.d20
   ;assigning grid out as the min and max in each dimension
   if ~keyword_set(nxout) then nxout = 50
   if ~keyword_set(nyout) then nyout = 50
   if ~keyword_set(xrange) then xrange = minmax(xarr)
   if ~keyword_set(yrange) then yrange = minmax(yarr)
   xdelt = (xrange[1]-xrange[0])/nxout
   ydelt = (yrange[1]-yrange[0])/nyout
   xout = xrange[0]+findgen(nxout)*xdelt ;dimension of nxout
   yout = yrange[0]+findgen(nyout)*ydelt ;dimension of nyout
   xgridout = rebin(xout,nxout,nyout) ;dimension of nxout X nyout
   ygridout = transpose(rebin(yout,nyout,nxout)) ; dimension of nxout X nyout

   ;interpolate probability into the new grid and make sure that the cumulative prob is still 1.
   probxygrid = dblarr(nxout,nyout,nobjs)
   volumearr = dblarr(nobjs)
   for i=0,nobjs-1 do begin
         probin = double(probdxdyarr[*,*,i])
         probgridnow = interp2d(probin,xarr[*,i],yarr[*,i],xout,yout,/grid) 
         wzero = where(probgridnow le 0.,czero,complement=wnotzero,ncomplement=cnotzero) ;the probability must not be smaller or equal to zero
         ;if czero gt 0 then print, 'obj',i,' fraction of zero prob:',float(czero)/(czero+cnotzero)
         if czero gt 0 then probgridnow(wzero) = 0 ;min(probgridnow(wnotzero))
         volume = int_tabulated_2d(xgridout,ygridout,probgridnow) ;this should be close to 1
         volumearr(i) = volume
         probgridnow = probgridnow/volume
         probxygrid[*,*,i] = probgridnow
         
         writecheckfits = 0
         if writecheckfits then begin 
            mkhdr,header,probgridnow,/image
            paraname = ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2']
            paraval = [xout[0],yout[0],0,0,xdelt,ydelt]
            for ii=0,n_elements(paraname)-1 do sxaddpar,header,paraname(ii),paraval(ii)
            sxaddpar,header,'CTYPE1','pixel'
            sxaddpar,header,'CYPTE2','pixel'
            if i eq 0 then mwrfits,probgridnow,'linfitprobgrid_checkinterpolategrid.fits',header,/create,/silent else $
                           mwrfits,probgridnow,'linfitprobgrid_checkinterpolategrid.fits',header,/silent
         endif
   endfor
   ;getting cumulative probability of all objects
    totalprob = total(probxygrid,3) ;summing over all objects
    sumlogprob = alog10(totalprob)
    wnan = where(~finite(sumlogprob),cwnan,complement=wgood)
    minlog = min(sumlogprob(wgood))
    if cwnan gt 0 then sumlogprob(wnan) = minlog;-2

   if keyword_set(outfitfile) then begin
      mkhdr,header,sumlogprob,/image
      paraname = ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2']
      paraval = [xout[0],yout[0],0,0,xdelt,ydelt]
      for ii=0,n_elements(paraname)-1 do sxaddpar,header,paraname(ii),paraval(ii)
      writefits,outfitfile,sumlogprob,header
   endif

   if keyword_set(outfitfile) then print,'now doing ',outfitfile
   xfit = xout(where(xout gt fitrange[0] and xout lt fitrange[1]))

   if keyword_set(gridsearch) then begin
      np0 = 25. 
      np1 = 25.
      p0arr = (p0range[1]-p0range[0])/(np0)*findgen(np0+1)+p0range[0]
      p1arr = (p1range[1]-p1range[0])/(np1)*findgen(np1+1)+p1range[0]
      if ~file_test('./linfitprobgrid_gridsearch',/directory) then file_mkdir,'linfitprobgrid_gridsearch'
      gridsearch_linfitprobgrid,xfit,xout,yout,probxygrid,p0arr,p1arr,'linfitprobgrid_gridsearch/gridsearchsdss_chisq.fits','linfitprobgrid_gridsearch',returnvalues,badid=badid,skipmapmaking=skipmapmaking
      stop
   endif else begin
   ;mcmc stuff
     ; if ~keyword_set(stepsize) then stepsize = [0.001,0.001]
      if ~keyword_set(nstep) then nstep = 3000
     ; if ~keyword_set(initial_guess) then initial_guess = [-0.4,0.06]
      mcmc_linfitprobgrid2,xfit,xout,yout,probxygrid,2,nstep,initial_guess,stepsize,['intercept','slope'],bestfitparam,plot=plot,badid=badid
   endelse

   bestfitparam2=bestfitparam
   bestintercepts = bestfitparam[*,0]
   bestslopes = bestfitparam[*,1]  


   returnval = bestfitparam2[1,*] ;median value
   sigmaout = [[bestintercepts],[bestslopes]]
   ;plotting
   if keyword_Set(outepsfile) then  plotprobgrid,sumlogprob,xout,yout,outepsfile,'age of univ','delta feh','sum log prob',bestintercepts,bestslopes
return,returnval

end
