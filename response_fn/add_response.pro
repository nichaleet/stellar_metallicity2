function add_response, lambda, spec, zage, element, abund,silent=silent
;Function to interpolate the modelled response fn from alf to the input
;lambda(xin) and add it to the input spectra(yin)
;PARAMETER
;lambda  : lambda in rest wl
;spec    : spectrum to be added with response fn
;zage    : array of [zmet,age (in Gyr)]
;element : name array of elements to be added. Choices are Fe,Mg,O,C,N,Si,Ca
;          (see response_interp.pro for updated available choices)
;abundace: float array of each element abundance
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
   zh = zage[0]
   age = zage[1]

   rstr = response_interp(element,abund,zh,age,silent=silent)
   ;interpolate to input lambda
   rout = interpol(rstr.rspec,rstr.lambda,lambda)
   specout = rout*spec
   return,specout 
end
