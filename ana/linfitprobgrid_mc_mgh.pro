function linfitprobgrid_mc_mgh,x,y1,y2,dxarr,dy1arr,dy2arr,probcube,sigmaout
;monte carlo to the probability distribution. Linfit to x and y1+y2
;INPUT::
;       x     = array of x values of size nobjs
;       y1     = array of y1 values of size nobjs
;       y2     = array of y2 values
;       dxarr = array of dx (so x could be x+dx) of size ndx X nobjs
;       dyarr = array of dy (so y could be y+dy) of size ndy X nobjs
;       probcubearr = array of probabilities of each associated x+dxarr,y+dxarr. The size is ndx X ndy X nobjs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   dimen = size(probcube,/dimensions)
   ngridx = dimen(0)
   ngridy1 = dimen(1)
   ngridy2 = dimen(2)
   nobjs = dimen(3)
 
   xarr = transpose(rebin(x,nobjs,ngridx))+dxarr ;dimension of ndx X nobjs
   y1arr = transpose(rebin(y1,nobjs,ngridy1))+dy1arr ;dimension of ndy X nobjs
   y2arr = transpose(rebin(y2,nobjs,ngridy2))+dy2arr ;dimension of ndy X nobjs

   ;turn the values of each object into 1 dimension array and 
   ;getting cumulative prob for each obj

   cumprobcubearr = fltarr(ngridx*ngridy1*ngridy2,nobjs)
   x_1dimarr = cumprobcubearr 
   y1_1dimarr = cumprobcubearr
   y2_1dimarr = cumprobcubearr

   for nn=0,nobjs-1 do begin
     curprob = probcube[*,*,*,nn]
     npix = n_elements(curprob)
     cumprob = 0
     for ii=0,npix-1 do begin
        if ii eq 0 then cumprob = curprob[ii]  else cumprob = cumprob+curprob[ii]
        ind = array_indices(curprob,ii)
        cumprobcubearr[ii,nn]=cumprob
        x_1dimarr[ii,nn] = xarr[ind(0),nn]
        y1_1dimarr[ii,nn] = y1arr[ind(1),nn]
        y2_1dimarr[ii,nn] = y2arr[ind(2),nn]
     endfor 
     cumprobcubearr[*,nn]=cumprobcubearr[*,nn]/cumprob
   endfor

   ;mc stuff
   nmc = 5000
   pararr = fltarr(2,nmc)
   sigmapararr = fltarr(2,nmc)
   for i=0,nmc-1 do begin
     randomval = randomu(seed,nobjs)
     xnow = fltarr(nobjs)
     ynow = fltarr(nobjs)
     for nn=0,nobjs-1 do begin
        loc = value_locate(cumprobcubearr[*,nn],randomval(nn))+1
        xnow[nn] = x_1dimarr[loc,nn]
        ynow[nn] = y1_1dimarr[loc,nn]+y2_1dimarr[loc,nn]
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
