; another program to correct for the detunes intensities
; for HMI_laser2.pro
; I wrote another code to make sure the first one is working properly...

OPENR,1,'CPT/T_OBS_710420_2.txt' ; contains the T_OBS values for FSNs 710420 to 714327

TOBS=FLTARR(3908)
READF,1,TOBS ; T_OBS
CLOSE,1
FSN=FINDGEN(3908)+710420


EXPOSURE=FLTARR(3908)
OPENR,1,'CPT/EXP_710420_2.txt' ; contains exposure times
READF,1,EXPOSURE ; T_OBS
CLOSE,1

OPENR,1,'CPT/20071015_LASER_NEWPORT_corrected.dat'
d=FLTARR(2,49599)
READF,1,d ; time and laser intensities
CLOSE,1

detunes_FSN=[[710420,710479], $
[710480,710539], $ 
[710540,710599], $ 
[710600,710659], $ 
[710660,710719], $ 
[710720,710779], $ 
[710780,710839], $ 
[710840,710899], $ 
[710900,710959], $ 
[710960,711019], $ 
[711020,711079], $ 
[711080,711139], $ 
[711140,711199], $ 
[711200,711259], $ 
[711260,711319], $ 
[711320,711379], $ 
[711380,711439], $ 
[711440,711499], $ 
[711500,711559], $ 
[711560,711619], $ 
[711620,711679], $ 
[711680,711739], $ 
[711740,711799], $ 
[711800,711859], $ 
[711860,711919], $ 
[711920,711979], $ 
[711980,712039], $ 
[712040,712099], $ 
[712100,712159], $ 
[712160,712219], $ 
[712220,712279], $ 
[712280,712339], $ 
[712408,712467], $ 
[712468,712527], $ 
[712528,712587], $ 
[712588,712647], $ 
[712648,712707], $
[712708,712767], $ 
[712768,712827], $ 
[712828,712887], $ 
[712888,712947], $ 
[712948,713007], $ 
[713008,713067], $ 
[713068,713127], $ 
[713128,713187], $ 
[713188,713247], $ 
[713248,713307], $ 
[713308,713367], $ 
[713368,713427], $ 
[713428,713487], $ 
[713488,713547], $ 
[713548,713607], $ 
[713608,713667], $ 
[713668,713727], $ 
[713728,713787], $ 
[713788,713847], $ 
[713848,713907], $ 
[713908,713967], $ 
[713968,714027], $ 
[714028,714087], $ 
[714088,714147], $ 
[714148,714207], $ 
[714208,714267], $ 
[714268,714327]] 

;minten=FLTARR(64)
;sinten=minten
;FOR i=0,63 DO BEGIN
;    a=WHERE(FSN GE detunes_FSN[0,i] AND FSN LE detunes_FSN[1,i],na)
;    b=WHERE(d[0,*] GE TOBS[a[0]] AND d[0,*] LE TOBS[a[na-1]],nb)
;    minten[i]=MEAN(d[1,b])
;    sinten[i]=SIGMA(d[1,b])
;    PLOT,d[0,b],(d[1,b]-minten[i])/minten[i],xst=1,yst=1,psym=3
;    PRINT,FSN[a[0]],FSN[a[na-1]],minten[i],sinten[i],(MAX(d[1,b])-MIN(d[1,b]))/MEAN(d[1,b])
;read,pause
;ENDFOR

minten=FLTARR(3908)
FOR i=0,3907 DO BEGIN
    a=WHERE(d[0,*] GE TOBS[i] AND d[0,*] LE TOBS[i]+EXPOSURE[i]/1000.)
    minten[i]=MEAN(d[1,a])
ENDFOR

read,pasue

END
