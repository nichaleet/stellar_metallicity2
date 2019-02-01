pro chrisholm2017_scaling
  ;mass loading factor from UV absorption line of nearby galaxies
  ;the paper provide the scaling for all mass. we calculate here 
  ;only for those gt 10^9 solar mass 
  mass = [9.3,9.5,10.1,10.3,10.7]
  nta = [0.91,2,1.2,0.23,0.54]
  ntaerr = [0.33,0.57,0.51,0.08,0.18]
  lognta = alog10(nta)
  logntaerr = abs(ntaerr/nta/alog(10.))
  par = linfit(mass,lognta,measure_errors=logntaerr,sigma=sigma)
  print, 'intercept    slope'
  print, par
  print, sigma
end
 
