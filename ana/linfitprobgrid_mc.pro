function linfitprobgrid_mc,x,y,dxarr,dyarr,probdxdyarr,sigmaout
;monte carlo to the probability distribution
;INPUT::
;       x     = array of x values of size nobjs
;       y     = array of y values of size nobjs
;       dxarr = array of dx (so x could be x+dx) of size ndx X nobjs
;       dyarr = array of dy (so y could be y+dy) of size ndy X nobjs
;       probdxdyarr = array of probabilities of each associated x+dxarr,y+dxarr. The size is ndx X ndy X nobjs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   dimen = size(probdxdyarr,/dimensions)
   ngridx = dimen(0)
   ngridy = dimen(1)
   nobjs = dimen(2)
 
   xarr = transpose(rebin(x,nobjs,ngridx))+dxarr ;dimension of ndx X nobjs
   yarr = transpose(rebin(y,nobjs,ngridy))+dyarr ;dimension of ndy X nobjs

   ;turn the values of each object into 1 dimension array and 
   ;getting cumulative prob for each obj

   cumprobdxdyarr = fltarr(ngridx*ngridy,nobjs)
   x1demarr = cumprobdxdyarr 
   y1demarr = cumprobdxdyarr
   for nn=0,nobjs-1 do begin
     curprob = probdxdyarr[*,*,nn]
     npix = n_elements(curprob)
     cumprob = 0
     for ii=0,npix-1 do begin
        if ii eq 0 then cumprob = curprob[ii]  else cumprob = cumprob+curprob[ii]
        ind = array_indices(curprob,ii)
        cumprobdxdyarr[ii,nn]=cumprob
        x1demarr[ii,nn] = xarr(ind(0),nn)
        y1demarr[ii,nn] = yarr(ind(1),nn)
     endfor 
     cumprobdxdyarr[*,nn]=cumprobdxdyarr[*,nn]/cumprob
   endfor

   ;mc stuff
   nmc = 1000
   pararr = fltarr(2,nmc)
   sigmapararr = fltarr(2,nmc)
   for i=0,nmc-1 do begin
     randomval = randomu(seed,nobjs)
     xnow = fltarr(nobjs)
     ynow = fltarr(nobjs)
     for nn=0,nobjs-1 do begin
        loc = value_locate(cumprobdxdyarr[*,nn],randomval(nn))+1
        xnow[nn] = x1demarr[loc,nn]
        ynow[nn] = y1demarr[loc,nn]
     endfor
     pararr[*,i]=linfit(xnow,ynow,sigma=sigma)
     sigmapararr[*,i] = sigma
   endfor
   set_plot,'x'
   !p.multi = [0,2,1]
   plothist,pararr[0,*]
   plothist,pararr[1,*]
   meanerr,pararr[0,*],sigmapararr[0,*],xmeanp0,sigmamp0,sigmadp0,sigmasp0
   meanerr,pararr[1,*],sigmapararr[1,*],xmeanp1,sigmamp1,sigmadp1,sigmasp1
   returnval = [xmeanp0,xmeanp1]
   sigmaout = [sigmasp0,sigmasp1]
;   stop
   return,returnval

end
