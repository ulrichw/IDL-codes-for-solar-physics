filename=STRARR(12)
filename[0]="/SUM2/D22154010/S00000/file.fits"
filename[1]="/SUM4/D22154012/S00000/file.fits"
filename[2]="/SUM6/D22154014/S00000/file.fits"
filename[3]="/SUM8/D22154016/S00000/file.fits"
filename[4]="/SUM0/D22154018/S00000/file.fits"
filename[5]="/SUM10/D22154020/S00000/file.fits"
filename[6]="/SUM12/D22154022/S00000/file.fits"
filename[7]="/SUM3/D22154024/S00000/file.fits"
filename[8]="/SUM5/D22154026/S00000/file.fits"
filename[9]="/SUM7/D22154028/S00000/file.fits"
filename[10]="/SUM9/D22154030/S00000/file.fits"
filename[11]="/SUM1/D22154032/S00000/file.fits"

LCP=fltarr(6,4096,4096)
RCP=LCP

FOR i=0,5 DO BEGIN
LCP[i,*,*]=readfits(filename[i*2])
RCP[i,*,*]=readfits(filename[i*2+1])
ENDFOR

a=WHERE(LCP EQ MIN(LCP))
LCP[a]=0.0
a=WHERE(RCP EQ MIN(RCP))
RCP[a]=0.0

LCP1=FLTARR(6)
RCP1=LCP1
FOR i=0,5 DO BEGIN
LCP1[i]=LCP[i,2048,2048]
RCP1[i]=RCP[i,2048,2048]
ENDFOR

LCP0=REFORM(rebin(LCP,6,1,1))
RCP0=REFORM(rebin(RCP,6,1,1))

plot,LCP1,psym=4
oplot,RCP1,col=180,psym=4

END
