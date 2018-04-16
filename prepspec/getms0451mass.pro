function getms0451mass, cat
;read combined spectra, which has moran catalog in one of the extensions (should be 6)
;with the info of colors
;calculate stellar mass using kcorrect then add the mass to one of the field in the structure
;INPUT: cat should be a structure containing B(VRIJK)_AUTO, F814W_AUTO,...

   cgals = n_elements(cat)
   color = dblarr(7,cgals)
   ;magnitudes are in AB mag
   color[0,*] = cat.NUV_AUTO
   color[1,*] = cat.B_AUTO
   color[2,*] = cat.V_AUTO
   color[3,*] = cat.R_AUTO
   color[4,*] = cat.I_AUTO
   color[5,*] = cat.K_AUTO
   color[6,*] = cat.F814W_AUTO

   badcolor = where(color eq 99. or color eq 0., cbadcolor)
   if cbadcolor gt 0 then color(badcolor) = 1./0.
   badkcolor = where(cat.k_auto gt 24.,cbadkcolor)
   if cbadkcolor gt 0 then color[5,badkcolor] = 1./0.

   colorerr = dblarr(7,cgals)      ;b,v,r,i,j,k,f814
   colorerr[0,*] = cat.NUV_AUTO_ERR
   colorerr[1,*] = cat.B_AUTO_ERR
   colorerr[2,*] = cat.V_AUTO_ERR
   colorerr[3,*] = cat.R_AUTO_ERR
   colorerr[4,*] = cat.I_AUTO_ERR
   colorerr[5,*] = cat.K_AUTO_ERR
   colorerr[6,*] = cat.F814W_AUTO_ERR
   baderr = where(finite(color) and colorerr eq 99,cbaderr)
   if cbaderr gt 0 then colorerr(baderr) = 0.02

   filterlist=['galex_NUV.par','cfh12k_B.par','cfh12k_V.par','cfh12k_R.par',$
               'cfh12k_I.par','lco_wirc_Ks.par','wfpc2_f814w.par']

   zlist = cat.z

   mass = dblarr(cgals)
   gotmass = bytarr(cgals)
   count = 0
   for i0=0,1 do for i1=0,1 do for i2=0,1 do for i3=0,1 do for i4=0,1 do $
   for i5=0,1 do for i6=0,1 do begin
      good = where(finite(color[0,*]) eq i0 and finite(color[1,*]) eq i1 and $
                   finite(color[2,*]) eq i2 and finite(color[3,*]) eq i3 and $
                   finite(color[4,*]) eq i4 and finite(color[5,*]) eq i5 and $
                   finite(color[6,*]) eq i6,cgood)
      filteron = where([i0,i1,i2,i3,i4,i5,i6],cfilteron)
      count += cgood 
      if cgood gt 0 and cfilteron gt 3 then begin
         filternow = filterlist(filteron)
         print, strjoin(filternow,' '),cgood, format = '("doing:",A," for ",I3," galaxies.")'
         colornow = color[*,good]
         colorerrnow = colorerr[*,good]
         colornow = colornow[filteron,*]
         colorerrnow = colorerr[filteron,*]
         zlistnow = zlist(good)
         kcorrect,colornow,colorerrnow,zlistnow,kcorrect1,chi2=chi2,filterlist=filternow,$
                 /magnitude,mass=massnow,b300=b300
         mass(good) = massnow
         gotmass(good) = 1
      endif      
   endfor
   if count ne cgals then stop

   mass = mass/0.52 ;fix for the hubble constant h=0.7
   mass = alog10(mass)
   return, mass   
end
