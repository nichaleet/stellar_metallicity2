pro make_targetlist
  ;read in table5 spectral detail
  readcol,'gdds_sa22_table5.txt',gdds_t5,agn,oii,oiii,hiBal,loBal,Fe2375,Fe2600,Mg2800,Mg2852,$
          HandK,Balmer,D4000,Temp,Class,conf,z_t5,$
          format='(A7,I1,I1,I1,I1,I1,I1,I1,I1,I1,I1,I1,I1,I1,I3,I2,F6.3)',comment='#'
  ;read in table4 
  readcol,'gdds_sa22_table4.txt',gdds_t4,z_t4,sp,ihr,imin,xsec,ideg,imn,xsc,conf,ovlap,$
          vmag,imag,zmag,kmag,format='(A7,f7.4,A2,I2,I2,f5.2,I3,I2,f5.2,I2,I1,F5.2,F5.2,F5.2,F5.2)',$
          comment='#'
  ;first select with redshift and spectral class
  zmin = 0.85 ;actually should be 0.8793 to get Mg5170 in Yband, but relaxed
  zmax = 1.176  ;to get Mg5170 in Yband
  good = where(z_t5 ge zmin and z_t5 le zmax and class lt 100 and agn eq 0, cgood)
  
  match,gdds_t5(good),gdds_t4,suba,subb,count=cmatch
  finalist = {gdds:gdds_t4(subb),agn:agn(good(suba)),oii:oii(good(suba)),oiii:oiii(good(suba)),$
              hibal:hiBal(good(suba)),lobal:loBal(good(suba)),fe2375:Fe2375(good(suba)),$
              fe2600:Fe2600(good(suba)),mg2800:Mg2800(good(suba)),mg2852:Mg2852(good(suba)),$
              handk:HandK(good(suba)),balmer:Balmer(good(suba)),d4000:D4000(good(suba)),$
              temp:Temp(good(suba)),class:Class(good(suba)),conf:conf(good(suba)),z:z_t5(good(suba)),$
              sp:sp(subb),ihr:ihr(subb),imin:imin(subb),xsec:xsec(subb),ideg:ideg(subb),$
              imn:imn(subb),xsc:xsc(subb),vmag:vmag(subb),imag:imag(subb),zmag:zmag(subb),$
              kmag:kmag(subb)}
  ;high to low priority: class 1,11,10
  priority = intarr(cmatch)
  for i=0,cmatch-1 do begin
      if finalist.class(i) eq 1 then basep = 200
      if finalist.class(i) eq 11 then basep = 100
      if finalist.class(i) eq 10 then basep = 50
      priority(i) = basep+fix((24-finalist.zmag(i))*10)
  endfor
  write_ds9_regionfile,(finalist.ihr+finalist.imin/60.+finalist.xsec/3600.)*15.,$
                        finalist.ideg+finalist.imn/60.+finalist.xsc/3600.,$
                        comment=finalist.gdds, filename='gdds_sa22_targets.reg',color='yellow'
 
  writecol,'gdds_sa22_targetlist.txt','sa'+finalist.gdds,priority,finalist.zmag,$
           finalist.ihr,finalist.imin,finalist.xsec,finalist.ideg,finalist.imn,finalist.xsc,$
           fltarr(cmatch)+2000.,fltarr(cmatch)+2000.,fltarr(cmatch),fltarr(cmatch),$
           fmt='(A,1x,I5,2x,F4.1,2x,I2,2x,I2,2x,F5.2,2x,I2,2x,I2,2x,F5.2,2x,F6.1,2x,F6.1,2x,F3.1,2x,F3.1)'
end
