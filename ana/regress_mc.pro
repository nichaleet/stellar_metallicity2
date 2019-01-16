function regress_mc, x,y,dyarr,probdyarr,const=const,sigmaconst=sigmaconst,$
               sigmacoeff=sigmacoeff,rsq=rsq,stderror=stderror,adj_rsq = adj_rsq
;PROPOSE
;Perform regression with monte carlo when probability distribution of y is not gaussian
;INPUT 
;x: x is array with [nx,nobjs] in size. Each x will be in the regression model
;      y = const+a0x0+a1x1+...anxn   
;y: y is array of what you fit (arrays of nobjs size)
;dyarr and probdyarr: they are [ngrid,nobjs] array 
;                     (each dyarr[*,i] and probdyarr[*,i] is dy to be added to y
;                      and the probability for ith object)
;NL Jan2019
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;setting
   nmc = 1000. ;do 1000 times
   
   dummy = size(x,/dimensions)
   nparams = dummy[0]
   nobjs = dummy[1]

   coeffarr = fltarr(nparams,nmc)
   sigmaarr = fltarr(nparams,nmc)
   constarr = fltarr(nmc)

   dummy = size(dyarr,/dimensions)
   ngrid = dummy[0]
   if nobjs ne dummy[1] then stop,'Stopped:REGRESS_MC-->wrong input dimensions'

   ;getting cumulative probability
   cumprob = fltarr(ngrid,nobjs)
   dy = fltarr(nobjs)
   for i=0,nobjs-1 do begin
      curyarr = dyarr[*,i]
      curprob = probdyarr[*,i]
      for j=1,ngrid-1 do cumprob[j,i] = int_tabulated(curyarr[0:j],curprob[0:j])
      dy(i) = 0.5*(interpol(dyarr[*,i],cumprob[*,i],0.84)-interpol(dyarr[*,i],cumprob[*,i],0.16))
   endfor

   ;run over nmc loops to get parameters
   for i=0,nmc-1 do begin
      ysamp = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         ysamp(j) = y(j)+interpol(dyarr[*,j],cumprob[*,j],randomarr(j))
      endfor
      coeffarr[*,i] = regress(x,ysamp,const=const,sigma=sigma)
      sigmaarr[*,i] = sigma
      constarr[i] = const
   endfor

   ;calculate the output parameters
   pout = mean(coeffarr,dimension=2)
   sigmacoeff = stddev(coeffarr,dimension=2)
   const = mean(constarr)
   sigmaconst = stddev(constarr)

   if arg_present(rsq) or arg_present(stderror) or arg_present(adj_rsq) then begin
      yfit = reform(const+x##pout)
      sse = total((yfit-y)^2)
      ybar = mean(y)
      sst = total((y-ybar)^2)
      stderror = sqrt(sse/nobjs)
      rsq = 1.-(sse/sst)
      adj_rsq = 1.-(1.-rsq)*(nobjs-1.)/(nobjs-nparams-1.)
   endif
return,pout
end
