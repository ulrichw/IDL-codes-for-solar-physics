list1 ='listSun060207_230329' ; in Calmode  VIGNETTING
list2 ='listSun060208_180841' ; in Calmode  VIGNETTINGlist
list3 ='listSun060224_225806' ; in Calmode  NO APPARENT PROBLEM

list4 ='listSun060616_205620' ; in Calmode
list5 ='listSun060616_211450' ; in Calmode
list6 ='listSun060616_224250' ; in Calmode

nx =128

RESTORE,'SEQUENCE_'+STRTRIM(list1,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
jeez1=fltarr(27)
for i=0,26 do jeez1[i]=total(total(imx[*,*,i],1),1)

RESTORE,'SEQUENCE_'+STRTRIM(list2,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
jeez2=fltarr(27)
for i=0,26 do jeez2[i]=total(total(imx[*,*,i],1),1)

RESTORE,'SEQUENCE_'+STRTRIM(list3,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
jeez3=fltarr(27)
for i=0,26 do jeez3[i]=total(total(imx[*,*,i],1),1)

RESTORE,'SEQUENCE_'+STRTRIM(list4,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
jeez4=fltarr(27)
for i=0,26 do jeez4[i]=total(total(imx[*,*,i],1),1)

RESTORE,'SEQUENCE_'+STRTRIM(list5,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
jeez5=fltarr(27)
for i=0,26 do jeez5[i]=total(total(imx[*,*,i],1),1)

RESTORE,'SEQUENCE_'+STRTRIM(list6,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
jeez6=fltarr(27)
for i=0,26 do jeez6[i]=total(total(imx[*,*,i],1),1)

PRINT,sigma(jeez1),MEAN(jeez1)
PRINT,sigma(jeez2),MEAN(jeez2)
PRINT,sigma(jeez3),MEAN(jeez3)
PRINT,sigma(jeez4),MEAN(jeez4)
PRINT,sigma(jeez5),MEAN(jeez5)
PRINT,sigma(jeez6),MEAN(jeez6)



jeez1=jeez1/MEAN(jeez1)
jeez2=jeez2/MEAN(jeez2)
jeez3=jeez3/MEAN(jeez3)
jeez4=jeez4/MEAN(jeez4)
jeez5=jeez5/MEAN(jeez5)
jeez6=jeez6/MEAN(jeez6)

WINDOW,0,retain=2
set_plot,'ps'
device,file='yo.ps'
plot,jeez1,yrange=[0.5,1.5],yst=1;,psym=10
oplot,jeez2;,psym=10
oplot,jeez3;,psym=10

oplot,jeez4,thick=3;,psym=10
oplot,jeez5,thick=3;,psym=10
oplot,jeez6,thick=3;,psym=10
device,/close



READ,pause

END



