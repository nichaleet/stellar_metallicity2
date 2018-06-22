pro make_alignstar_belli14
 readcol,'/scr2/nichal/workspace4/imagefile/Belli14/alignstars.tsv',_RAJ2000,_DEJ2000,RAJ2000,DEJ2000,twoMASS,Jmag,e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,Qflg,Rflg,Bflg,Cflg,Xflg,Aflg,format='F,F,F,F,A,F,F,F,F,F,F,A,F,F,F,F,F'
        radec,raj2000,dej2000,ihr,imin,xsec,ideg,imn,xsc
        nstars=n_elements(twomass)
        starnames='star'+strtrim(string(indgen(nstars)+1),2)+'b'
        star_pri = intarr(nstars)-1
        starlist = {name:starnames,prior:star_pri,mag:jmag,rahr:ihr,ramin:imin,rasec:xsec,decdeg:ideg,decmin:imn,decsec:xsc}
        final_list = [soa2aos(starlist)]

        ;create a mosfire list
        nobjs = n_elements(final_list)
        twokarr = fltarr(nobjs)+2000.
        zeroarr = fltarr(nobjs)
        writecol,'belli14_alignstarlist.txt',final_list.name,final_list.prior,final_list.mag,final_list.rahr,final_list.ramin,final_list.rasec,final_list.decdeg,final_list.decmin,final_list.decsec,twokarr,twokarr,zeroarr,zeroarr,fmt='(A,1x,I4,2x,F4.1,2x,I2,2x,I2,2x,F5.2,2x,I2,2x,I2,2x,F5.2,2x,F6.1,2x,F6.1,2x,F3.1,2x,F3.1)'
end
