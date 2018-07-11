pro combinesps,files,nameout,dir=dir
;e.g. combinesps,['sps_fit01.fits.gz','sps_fit02.fits.gz','sps_fit03.fits.gz','sps_fit04.fits.gz'],'sps_fitall.fits',dir='/scr2/nichal/workspace4/test_shortwl/datasps/test_shortwl_cl1604/'
;combinesps,['sps_fit01.fits.gz','sps_fit02.fits.gz','sps_fit03.fits.gz'],'sps_fitall.fits',dir='/scr2/nichal/workspace4/test_shortwl/datasps/test_fullwl_cl1604/'
   ;use the first file as the main
   if ~keyword_set(dir) then dir = './'
   nfiles = n_Elements(files)
   if strmid(dir,strlen(dir)-1,1) ne '/' then dir=dir+'/'
   nameout = dir+nameout
   str = mrdfits(dir+files[0],1,/silent)
   if nfiles gt 1 then for i=1,nfiles-1 do begin
      if file_test(dir+files[i]) eq 0 then stop,'cannot find file '+dir+files[i]
      newstr = mrdfits(dir+files[i],1,/silent)
      wreplace = where(newstr.chisq lt str.chisq,nreplace)
      if nreplace gt 0 then begin
         str(wreplace) = newstr(wreplace) 
      endif
   endfor
   mwrfits,str,nameout,/create,/silent
   spawn, 'gzip -f '+nameout
   print, 'done writing '+ nameout
   ;stop
end
