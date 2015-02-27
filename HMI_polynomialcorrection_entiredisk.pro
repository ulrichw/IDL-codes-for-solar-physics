; idl code that calculates the DATAMEDN-OBS_VR VALUE AT 5 DIFFERENT LOCATIONS
; SOLAR CENTER, AND 4 AREAS NEAR LIMB

PRO HMI_polynomialcorrection_entiredisk,draw

IF draw EQ 1 THEN GOTO,draw

; TEMP must include the filepaths
; TEMP2 must include OBS_VR,CRPIX1,CRPIX2,RSUN_OBS,CDELT1
nfile=19200 ; 10 days worth of data

OPENR,1,'temp'
filename=STRARR(nfile)
READF,1,filename
CLOSE,1
OPENR,1,'temp2'
data=FLTARR(5,nfile)
READF,1,data
CLOSE,1

DATAMEDN=FLTARR(5,nfile)
rotation0=FLTARR(nfile)
rotation1=FLTARR(nfile)

issue=0l

FOR i=0,nfile-1 DO BEGIN
    PRINT,i
    IF(FINITE(data[0,i]) NE 0) THEN BEGIN
        d=fitsio_read_image(STRTRIM(filename[i],1)+STRTRIM("/Dopplergram.fits",1))
        
        x0=[data[1,i]-1760,data[1,i]+1760,data[1,i],data[1,i],data[1,i]]
        y0=[data[2,i],data[2,i],data[2,i]-1760,data[2,i]+1760,data[2,i]]
        rotation0[i]=2000.0*ABS(X0[0]-data[1,i])/(data[3,i]/data[4,i])
        rotation1[i]=2000.0*ABS(X0[1]-data[1,i])/(data[3,i]/data[4,i])
        DATAMEDN[0,i]=MEDIAN(d[x0[0]-50:x0[0]+50,y0[0]-50:y0[0]+50]) ; left limb (west)
        DATAMEDN[1,i]=MEDIAN(d[x0[1]-50:x0[1]+50,y0[1]-50:y0[1]+50]) ; right limb (east)
        DATAMEDN[2,i]=MEDIAN(d[x0[2]-50:x0[2]+50,y0[2]-50:y0[2]+50]) ; down (north)
        DATAMEDN[3,i]=MEDIAN(d[x0[3]-50:x0[3]+50,y0[3]-50:y0[3]+50]) ; up (south)
        DATAMEDN[4,i]=MEDIAN(d[x0[4]-50:x0[4]+50,y0[4]-50:y0[4]+50]) ; disk center
    ENDIF ELSE BEGIN
        data[*,i]=data[*,i+1]
        filename[i]=filename[i+1]
        i=i-1
        issue=issue+1
    ENDELSE

END

PRINT,'ISSUES:',issue
SAVE,data,datamedn,rotation0,rotation1,FILE='polynomialcorrection_20101010-21.bin'

draw:

SET_PLOT,'PS'
!P.MULTI=0
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=22,ysize=20,/color,bits=24
LOADCT,3

restore,'polynomialcorrection_20110301-11.bin'
;restore,'polynomialcorrection_20101010-21.bin'
a=WHERE(datamedn[0,*] ne 0.0)
b=SORT(data[0,a])
plot ,data[0,a[b]]+rotation0[a[b]],SMOOTH(datamedn[0,a[b]]-data[0,a[b]]-rotation0[a[b]],25,/edge_truncate),xrange=[-5500,5500],yrange=[-250,250],yst=1,charsize=1.0,tit='!17 2011/03/01-2011/03/11',xtit='OBSVR (m/s)',ytit='DATAMEDN-OBSVR (m/s)',xst=1
oplot,data[0,a[b]]-rotation1[a[b]],SMOOTH(datamedn[1,a[b]]-data[0,a[b]]+rotation1[a[b]],25,/edge_truncate)
oplot,data[0,a[b]],SMOOTH(datamedn[2,a[b]]-data[0,a[b]],25,/edge_truncate)
oplot,data[0,a[b]],SMOOTH(datamedn[3,a[b]]-data[0,a[b]],25,/edge_truncate)
oplot,data[0,a[b]],SMOOTH(datamedn[4,a[b]]-data[0,a[b]],25,/edge_truncate),col=180 ; disk center

RESTORE,'polynomialcorrection_20101010.bin'
data0=data
rotation00=rotation0
rotation01=rotation1
datamedn0=datamedn
RESTORE,'polynomialcorrection_20110101.bin'
data1=data
rotation10=rotation0
rotation11=rotation1
datamedn1=datamedn
RESTORE,'polynomialcorrection_20110301.bin'
data2=data
rotation20=rotation0
rotation21=rotation1
datamedn2=datamedn
RESTORE,'polynomialcorrection_20110701.bin'
data3=data
rotation30=rotation0
rotation31=rotation1
datamedn3=datamedn

!P.MULTI=[0,2,2]

plot ,data0[0,*]+rotation00[*],SMOOTH(datamedn0[0,*]-data0[0,*]-rotation00[*],15,/edge_truncate),xrange=[-5500,5500],yrange=[-250,250],yst=1,charsize=1.0,tit='!17 2010/10/10',xtit='OBSVR (m/s)',ytit='DATAMEDN-OBSVR (m/s)',xst=1
oplot,data0[0,*]-rotation01[*],SMOOTH(datamedn0[1,*]-data0[0,*]+rotation01[*],15,/edge_truncate)
oplot,data0[0,*],SMOOTH(datamedn0[2,*]-data0[0,*],15,/edge_truncate)
oplot,data0[0,*],SMOOTH(datamedn0[3,*]-data0[0,*],15,/edge_truncate)
oplot,data0[0,*],SMOOTH(datamedn0[4,*]-data0[0,*],15,/edge_truncate)

plot ,data1[0,*]+rotation10[*],SMOOTH(datamedn1[0,*]-data1[0,*]-rotation10[*],15,/edge_truncate),xrange=[-5500,5500],yrange=[-250,250],yst=1,charsize=1.0,tit='!17 2011/01/01',xtit='OBSVR (m/s)',ytit='DATAMEDN-OBSVR (m/s)',xst=1
oplot,data1[0,*]-rotation11[*],SMOOTH(datamedn1[1,*]-data1[0,*]+rotation11[*],15,/edge_truncate)
oplot,data1[0,*],SMOOTH(datamedn1[2,*]-data1[0,*],15,/edge_truncate)
oplot,data1[0,*],SMOOTH(datamedn1[3,*]-data1[0,*],15,/edge_truncate)
oplot,data1[0,*],SMOOTH(datamedn1[4,*]-data1[0,*],15,/edge_truncate)

plot ,data2[0,*]+rotation20[*],SMOOTH(datamedn2[0,*]-data2[0,*]-rotation20[*],15,/edge_truncate),xrange=[-5500,5500],yrange=[-250,250],yst=1,charsize=1.0,tit='!17 2011/03/01',xtit='OBSVR (m/s)',ytit='DATAMEDN-OBSVR (m/s)',xst=1
oplot,data2[0,*]-rotation21[*],SMOOTH(datamedn2[1,*]-data2[0,*]+rotation21[*],15,/edge_truncate)
oplot,data2[0,*],SMOOTH(datamedn2[2,*]-data2[0,*],15,/edge_truncate)
oplot,data2[0,*],SMOOTH(datamedn2[3,*]-data2[0,*],15,/edge_truncate)
oplot,data2[0,*],SMOOTH(datamedn2[4,*]-data2[0,*],15,/edge_truncate)

plot ,data3[0,*]+rotation30[*],SMOOTH(datamedn3[0,*]-data3[0,*]-rotation30[*],15,/edge_truncate),xrange=[-5500,5500],yrange=[-250,250],yst=1,charsize=1.0,tit='!17 2011/07/01',xtit='OBSVR (m/s)',ytit='DATAMEDN-OBSVR (m/s)',xst=1
oplot,data3[0,*]-rotation31[*],SMOOTH(datamedn3[1,*]-data3[0,*]+rotation31[*],15,/edge_truncate)
oplot,data3[0,*],SMOOTH(datamedn3[2,*]-data3[0,*],15,/edge_truncate)
oplot,data3[0,*],SMOOTH(datamedn3[3,*]-data3[0,*],15,/edge_truncate)
oplot,data3[0,*],SMOOTH(datamedn3[4,*]-data3[0,*],15,/edge_truncate)

DEVICE,/CLOSE
SET_PLOT,'X'

READ,pause

END
