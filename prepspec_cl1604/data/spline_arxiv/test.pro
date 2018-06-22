pro test
nmax = 0
file = file_search('spline*.fits.gz',count=cfile)
for i=0,cfile-1 do begin
   spec = mrdfits(file[i],1)
   nlambda = n_elements(spec.lambda)
   if nlambda gt nmax then nmax=nlambda
endfor
print, nmax

stop
end
