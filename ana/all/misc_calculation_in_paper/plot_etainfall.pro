pro plot_etainfall
  ;kappa is stellar metallicity in unit of yield
  ;alpha is infall metallicity in unit of solar metallicity

  alpha = [0.01,0.05,0.1]
  logkappa = findgen(100)/100*(-0.2+1.5)-1.5
  kappa = 10.^logkappa
  color = ['ygb5','seagreen','firebrick']

  ;metallicity assume yield=3*solar 
  metal = alog10(kappa*3.)
  set_plot,'ps'
  !p.font = 0
  !p.charsize = 1.2
  ntaletter = "150B
  proptoletter = "265B
  sunsym = sunsymbol()
  xrange = minmax(logkappa)
  yrange = [0.5,1.1]
  psname='nta_overestimate.eps'
  device, filename = psname,xsize = 15,ysize = 10, $
       xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  for i=0,n_elements(alpha)-1 do begin
    etaovers = (1.-kappa)/((1./(1.-alpha[i]))-kappa)
    if i eq 0 then begin
        plot,logkappa,etaovers,xtitle='log(Z!L*!N/yield)',ytitle='Fraction !9'+string(ntaletter)+'!x under-estimated',xrange=xrange,yrange=yrange,xstyle=9,ystyle=1,/nodata,position=[0.15,0.15,0.95,0.85],psym=1
    endif
    oplot,logkappa,etaovers,color=fsc_color(color[i]),linestyle=i,thick=2
  endfor
  axis,xaxis=1,xtitle='[Z!L*!N/Z'+sunsym+']',xrange=minmax(metal),xstyle=1
  xyouts,-1.45,0.7,'Inflow metallicity'
  al_legend,['0.01 Z!L*!N','0.05 Z!L*!N','  0.1 Z!L*!N'],colors=color,linestyle=[0,1,2],thick=2,box=0,/bottom,/left,linsize=0.5
  device,/close
  stop
end
