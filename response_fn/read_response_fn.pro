pro read_response_fn
   ;change the alf text file to structure format
   krpa_file = file_search('/scr2/nichal/workspace4/alf/infiles/atlas_ssp_t*.krpa.s100',count=cfile)
   response_arr=[]
   for i=0,cfile-1 do begin
      readcol,krpa_file(i),lam, Solar, Nap, Nam, Cap, Cam, Fep, Fem, Cp, Cm, aFep,$
         Np, Nm, asFep, Tip, Tim, Mgp, Mgm, Sip, Sim,Tp, Tm, Crp, Mnp, Bap, Bam,$
         Nip, Cop, Eup, Srp, Kp, Vp, Cup, Nap6, Nap9,comment='#' 
      extens = strsplit(krpa_file(i),'_',/extract)
      agegyr = float(strmid(extens[2],1))
      zmetsign = (strmid(extens[3],1,1) eq 'p') ?  1.: -1.
      zmet = zmetsign*float(strmid(extens[3],2,4))
      response = {agegyr:agegyr,zmet:zmet,lam:lam, Solar:solar, Nap:Nap, Nam:Nam,$
                  Cap:Cap, Cam:Cam, Fep:Fep,$
                  Fem:Fem, Cp:Cp, Cm:Cm, aFep:aFep,Np:Np, Nm:Nm, asFep:asFep, Tip:Tip,$
                  Tim:Tim, Mgp:Mgp, Mgm:Mgm, Sip:Sip, Sim:Sim,Tp:Tp, Tm:Tm, Crp:Crp,$
                  Mnp:Mnp, Bap:Bap, Bam:Bam,Nip:Nip, Cop:Cop, Eup:Eup, Srp:Srp, Kp:Kp,$
                  Vp:Vp, Cup:Cup, Nap6:Nap6, Nap9:Nap9} 
      response_arr=[response_arr,response]
   endfor
   mwrfits,response_arr,'atlas_ssp_abund_krp.fits',/create,/silent
end
