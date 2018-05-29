pro makems0451regionfile,list,cat
   nobj = n_elements(list)
   ra = [] 
   dec = []
   comment = []
   for i=0,nobj-1 do begin
      nspec = list[i].count
      if nspec gt 1 then namelist = strtrim(list[i].objname[0:nspec-1],2)+'_'+strtrim(string(i),2)+'_'+strtrim(string(indgen(nspec),format='(I)'),2) $
                    else namelist = strtrim(list[i].objname[0:nspec-1],2)+'_'+strtrim(string(i),2)     
      ra = [ra,list[i].ra[0:nspec-1]]
      dec= [dec,list[i].dec[0:nspec-1]]
      comment = [comment,namelist]
      print, namelist
   endfor
   write_ds9_regionfile, ra, dec,comment=comment, filename='ms0451_observed_spec.reg', symbol='circle', color='green'
   write_ds9_regionfile, cat.ra,cat.dec,comment=strtrim(cat.specname,2),filename='ms0451_cat.reg',symbol='circle',color='red'
end
