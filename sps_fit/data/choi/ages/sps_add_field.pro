pro sps_add_field,filein,tag,data,itag=itag,fileout=fileout
;add field to array of structure
;  itag, where itag=0 indicates that the new field should be created at
;  the beginning of the structure, itag=1 indicates that the new field
;  should be the second field in the structure, etc.
   if ~keyword_set(fileout) then fileout = filein
   str = mrdfits(filein,1)
   nstr = n_Elements(str)
   temp = str[0]
   if n_Elements(data) eq 1 then data = rebin(data,nstr)
   if ~keyword_set(itag) then itag = n_tags(temp)
   struct_add_field,temp,tag,data[0],itag=itag
   strout = replicate({temp}, nspec)
   for i=0,nstr-1 do begin
      temp = strout[i]
      struct_assign,str[i],temp,/nozero
      strout[i] = temp
   endfor
   strout.(itag) = data

end
