; tracked cube of quiet Sun at disk center: 512x512x400 (IDENT = QUIET
; SUN, T_START: 2010.05.29_4:15_TAI); 
LCP2=fitsio_read_image("/SUM12/D156481812/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM0/D156483395/S00000/trackedLCP3.fits")
OPENR,1,'velocitiesQUIETSUN' ; file containing OBS_VR
vel=FLTARR(400)
READF,1,vel
CLOSE,1
L3=TOTAL(TOTAL(LCP3[246:266,246:266,*],1),1)/21./21.
L2=TOTAL(TOTAL(LCP2[246:266,246:266,*],1),1)/21./21.
LCP3=0
LCP2=0
LCP1=fitsio_read_image("/SUM17/D156480000/S00000/trackedLCP1.fits")
LCP4=fitsio_read_image("/SUM7/D156484565/S00000/trackedLCP4.fits")
L1=TOTAL(TOTAL(LCP1[246:266,246:266,*],1),1)/21./21.
L4=TOTAL(TOTAL(LCP4[246:266,246:266,*],1),1)/21./21.
LCP1=0
LCP4=0
LCP0=fitsio_read_image("/SUM9/D156478242/S00000/trackedLCP0.fits")
LCP5=fitsio_read_image("/SUM0/D156485787/S00000/trackedLCP5.fits")
L0=TOTAL(TOTAL(LCP0[246:266,246:266,*],1),1)/21./21.
L5=TOTAL(TOTAL(LCP5[246:266,246:266,*],1),1)/21./21.
LCP0=0
LCP5=0

dlamdv=6173.3433/3.d8
lam=-vel*dlamdv
plot ,lam+0.172,L0,xrange=[-0.3,0.3],xst=1
oplot,lam+0.1032,L1
oplot,lam+0.0344,L2
oplot,lam-0.0344,L3
oplot,lam-0.1032,L4
oplot,lam-0.172,L5


; tracked cube of quiet Sun at disk center: 128x128x961 (IDENT = QUIET
; SUN 2, T_START: 2010.05.29_0:0_TAI); 
LCP2=fitsio_read_image("/SUM5/D156800752/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM19/D156803850/S00000/trackedLCP3.fits")
LCP1=fitsio_read_image("/SUM4/D156798543/S00000/trackedLCP1.fits")
LCP4=fitsio_read_image("/SUM13/D156806880/S00000/trackedLCP4.fits")
LCP0=fitsio_read_image("/SUM9/D156794960/S00000/trackedLCP0.fits")
LCP5=fitsio_read_image("/SUM17/D156810196/S00000/trackedLCP5.fits")
L0=TOTAL(TOTAL(LCP0[54:74,54:74,*],1),1)/21./21.
L1=TOTAL(TOTAL(LCP1[54:74,54:74,*],1),1)/21./21.
L2=TOTAL(TOTAL(LCP2[54:74,54:74,*],1),1)/21./21.
L3=TOTAL(TOTAL(LCP3[54:74,54:74,*],1),1)/21./21.
L4=TOTAL(TOTAL(LCP4[54:74,54:74,*],1),1)/21./21.
L5=TOTAL(TOTAL(LCP5[54:74,54:74,*],1),1)/21./21.

OPENR,1,'velocitiesQUIETSUN2' ; file containing OBS_VR
vel=FLTARR(961)
READF,1,vel
CLOSE,1
dlamdv=6173.3433/3.d8
lam   =-vel*dlamdv
plot ,lam+0.172,L0,xrange=[-0.3,0.3],xst=1
oplot,lam+0.1032,L1
oplot,lam+0.0344,L2
oplot,lam-0.0344,L3
oplot,lam-0.1032,L4
oplot,lam-0.172,L5

;UNTRACKED cube of quiet Sun at disk center and 2 other locations: 4096x4096x961 (IDENT = QUIET
; SUN 2, T_START: 2010.05.29_0:0_TAI);
list=STRARR(961)
OPENR,1,'listoflev1pmay29'
READF,1,list
CLOSE,1

L0b=FLTARR(961)
L1b=L0b
L2b=L0b
L3b=L0b
L4b=L0b
L5b=L0b
L0c=FLTARR(961)
L1c=L0c
L2c=L0c
L3c=L0c
L4c=L0c
L5c=L0c
L0d=FLTARR(961)
L1d=L0d
L2d=L0d
L3d=L0d
L4d=L0d
L5d=L0d
L0e=FLTARR(961)
L1e=L0e
L2e=L0e
L3e=L0e
L4e=L0e
L5e=L0e
L0f=FLTARR(961)
L1f=L0f
L2f=L0f
L3f=L0f
L4f=L0f
L5f=L0f

CRPIX1=2042.1-1.
CRPIX2=2047.1-1.
CRPIX3=1042.1-1.
CRPIX4=2047.1-1.
CRPIX5=542.1-1.
CRPIX6=2047.1-1.
CRPIX7=1542.1-1.
CRPIX8=2047.1-1.
CRPIX9=190.1-1.
CRPIX10=2047.1-1.

FOR i=0,960 DO BEGIN 
    PRINT,i 
    d=fitsio_read_image(list[i]+"/LCP0.fits") 
    L0b[i]=TOTAL(TOTAL(d[CRPIX1-10:CRPIX1+10,CRPIX2-10:CRPIX2+10],1),1)/21./21. 
    L0c[i]=TOTAL(TOTAL(d[CRPIX3-10:CRPIX3+10,CRPIX4-10:CRPIX4+10],1),1)/21./21. 
    L0d[i]=TOTAL(TOTAL(d[CRPIX5-10:CRPIX5+10,CRPIX6-10:CRPIX6+10],1),1)/21./21.
    L0e[i]=TOTAL(TOTAL(d[CRPIX7-10:CRPIX7+10,CRPIX8-10:CRPIX8+10],1),1)/21./21.
    L0f[i]=TOTAL(TOTAL(d[CRPIX9-10:CRPIX9+10,CRPIX10-10:CRPIX10+10],1),1)/21./21.

     d=fitsio_read_image(list[i]+"/LCP1.fits") 
    L1b[i]=TOTAL(TOTAL(d[CRPIX1-10:CRPIX1+10,CRPIX2-10:CRPIX2+10],1),1)/21./21.
    L1c[i]=TOTAL(TOTAL(d[CRPIX3-10:CRPIX3+10,CRPIX4-10:CRPIX4+10],1),1)/21./21. 
    L1d[i]=TOTAL(TOTAL(d[CRPIX5-10:CRPIX5+10,CRPIX6-10:CRPIX6+10],1),1)/21./21.
    L1e[i]=TOTAL(TOTAL(d[CRPIX7-10:CRPIX7+10,CRPIX8-10:CRPIX8+10],1),1)/21./21.
    L1f[i]=TOTAL(TOTAL(d[CRPIX9-10:CRPIX9+10,CRPIX10-10:CRPIX10+10],1),1)/21./21.

    d=fitsio_read_image(list[i]+"/LCP2.fits") 
    L2b[i]=TOTAL(TOTAL(d[CRPIX1-10:CRPIX1+10,CRPIX2-10:CRPIX2+10],1),1)/21./21. 
    L2c[i]=TOTAL(TOTAL(d[CRPIX3-10:CRPIX3+10,CRPIX4-10:CRPIX4+10],1),1)/21./21. 
    L2d[i]=TOTAL(TOTAL(d[CRPIX5-10:CRPIX5+10,CRPIX6-10:CRPIX6+10],1),1)/21./21. 
    L2e[i]=TOTAL(TOTAL(d[CRPIX7-10:CRPIX7+10,CRPIX8-10:CRPIX8+10],1),1)/21./21.
    L2f[i]=TOTAL(TOTAL(d[CRPIX9-10:CRPIX9+10,CRPIX10-10:CRPIX10+10],1),1)/21./21.

    d=fitsio_read_image(list[i]+"/LCP3.fits") 
    L3b[i]=TOTAL(TOTAL(d[CRPIX1-10:CRPIX1+10,CRPIX2-10:CRPIX2+10],1),1)/21./21. 
    L3c[i]=TOTAL(TOTAL(d[CRPIX3-10:CRPIX3+10,CRPIX4-10:CRPIX4+10],1),1)/21./21. 
    L3d[i]=TOTAL(TOTAL(d[CRPIX5-10:CRPIX5+10,CRPIX6-10:CRPIX6+10],1),1)/21./21. 
    L3e[i]=TOTAL(TOTAL(d[CRPIX7-10:CRPIX7+10,CRPIX8-10:CRPIX8+10],1),1)/21./21.
    L3f[i]=TOTAL(TOTAL(d[CRPIX9-10:CRPIX9+10,CRPIX10-10:CRPIX10+10],1),1)/21./21.

    d=fitsio_read_image(list[i]+"/LCP4.fits") 
    L4b[i]=TOTAL(TOTAL(d[CRPIX1-10:CRPIX1+10,CRPIX2-10:CRPIX2+10],1),1)/21./21. 
    L4c[i]=TOTAL(TOTAL(d[CRPIX3-10:CRPIX3+10,CRPIX4-10:CRPIX4+10],1),1)/21./21. 
    L4d[i]=TOTAL(TOTAL(d[CRPIX5-10:CRPIX5+10,CRPIX6-10:CRPIX6+10],1),1)/21./21. 
    L4e[i]=TOTAL(TOTAL(d[CRPIX7-10:CRPIX7+10,CRPIX8-10:CRPIX8+10],1),1)/21./21.
    L4f[i]=TOTAL(TOTAL(d[CRPIX9-10:CRPIX9+10,CRPIX10-10:CRPIX10+10],1),1)/21./21.

    d=fitsio_read_image(list[i]+"/LCP5.fits") 
    L5b[i]=TOTAL(TOTAL(d[CRPIX1-10:CRPIX1+10,CRPIX2-10:CRPIX2+10],1),1)/21./21. 
    L5c[i]=TOTAL(TOTAL(d[CRPIX3-10:CRPIX3+10,CRPIX4-10:CRPIX4+10],1),1)/21./21. 
    L5d[i]=TOTAL(TOTAL(d[CRPIX5-10:CRPIX5+10,CRPIX6-10:CRPIX6+10],1),1)/21./21. 
    L5e[i]=TOTAL(TOTAL(d[CRPIX7-10:CRPIX7+10,CRPIX8-10:CRPIX8+10],1),1)/21./21.
    L5f[i]=TOTAL(TOTAL(d[CRPIX9-10:CRPIX9+10,CRPIX10-10:CRPIX10+10],1),1)/21./21.

ENDFOR


OPENR,1,'velocitiesQUIETSUN2' ; file containing OBS_VR
vel=FLTARR(961)
READF,1,vel
CLOSE,1
dlamdv=6173.3433/3.d8
lam   =-vel*dlamdv

RESTORE,'improvedline.bin' ; from HMI_filter_lookup.pro
!p.multi=[0,2,3]
plot ,lam+0.172 ,L0b/MAX(L0b),xrange=[-0.22,0.22],xst=1,yrange=[0.53,1],yst=1
oplot,lam+0.1032,L1b/MAX(L0b)
oplot,lam+0.0344,L2b/MAX(L0b)
oplot,lam-0.0344,L3b/MAX(L0b)
oplot,lam-0.1032,L4b/MAX(L0b)
oplot,lam-0.172 ,L5b/MAX(L0b)
oplot,wave+0.1720,line[0:300]    -0.01,col=180,thick=2
oplot,wave+0.1032,line[301:601]  -0.01,col=180,thick=2
oplot,wave+0.0344,line[602:902]  -0.01,col=180,thick=2
oplot,wave-0.0344,line[903:1203] -0.01,col=180,thick=2
oplot,wave-0.1032,line[1204:1504]-0.01,col=180,thick=2
oplot,wave-0.1720,line[1505:1805]-0.01,col=180,thick=2

plot ,lam+0.172,L0e/MAX(L0e),xrange=[-0.22,0.22],xst=1,yrange=[0.53,1],yst=1
oplot,lam+0.1032,L1e/MAX(L0e)
oplot,lam+0.0344,L2e/MAX(L0e)
oplot,lam-0.0344,L3e/MAX(L0e)
oplot,lam-0.1032,L4e/MAX(L0e)
oplot,lam-0.172 ,L5e/MAX(L0e)
oplot,wave+0.1720+0.0125,line[0:300]    -0.01,col=180,thick=2
oplot,wave+0.1032+0.0125,line[301:601]  -0.01,col=180,thick=2
oplot,wave+0.0344+0.0125,line[602:902]  -0.01,col=180,thick=2
oplot,wave-0.0344+0.0125,line[903:1203] -0.01,col=180,thick=2
oplot,wave-0.1032+0.0125,line[1204:1504]-0.01,col=180,thick=2
oplot,wave-0.1720+0.0125,line[1505:1805]-0.01,col=180,thick=2

plot ,lam+0.172,L0c/MAX(L0c),xrange=[-0.22,0.22],xst=1,yrange=[0.53,1],yst=1
oplot,lam+0.1032,L1c/MAX(L0c)
oplot,lam+0.0344,L2c/MAX(L0c)
oplot,lam-0.0344,L3c/MAX(L0c)
oplot,lam-0.1032,L4c/MAX(L0c)
oplot,lam-0.172 ,L5c/MAX(L0c)
oplot,wave+0.1720+0.025,line[0:300]    -0.01,col=180,thick=2
oplot,wave+0.1032+0.025,line[301:601]  -0.01,col=180,thick=2
oplot,wave+0.0344+0.025,line[602:902]  -0.01,col=180,thick=2
oplot,wave-0.0344+0.025,line[903:1203] -0.01,col=180,thick=2
oplot,wave-0.1032+0.025,line[1204:1504]-0.01,col=180,thick=2
oplot,wave-0.1720+0.025,line[1505:1805]-0.01,col=180,thick=2

plot ,lam+0.172,L0d/MAX(L0d),xrange=[-0.22,0.22],xst=1,yrange=[0.53,1],yst=1
oplot,lam+0.1032,L1d/MAX(L0d)
oplot,lam+0.0344,L2d/MAX(L0d)
oplot,lam-0.0344,L3d/MAX(L0d)
oplot,lam-0.1032,L4d/MAX(L0d)
oplot,lam-0.172 ,L5d/MAX(L0d)
oplot,wave+0.1720+0.037,line[0:300]    -0.01,col=180,thick=2
oplot,wave+0.1032+0.037,line[301:601]  -0.01,col=180,thick=2
oplot,wave+0.0344+0.037,line[602:902]  -0.01,col=180,thick=2
oplot,wave-0.0344+0.037,line[903:1203] -0.01,col=180,thick=2
oplot,wave-0.1032+0.037,line[1204:1504]-0.01,col=180,thick=2
oplot,wave-0.1720+0.037,line[1505:1805]-0.01,col=180,thick=2

plot ,lam+0.172,L0f/MAX(L0f),xrange=[-0.22,0.22],xst=1,yrange=[0.53,1],yst=1
oplot,lam+0.1032,L1f/MAX(L0f)
oplot,lam+0.0344,L2f/MAX(L0f)
oplot,lam-0.0344,L3f/MAX(L0f)
oplot,lam-0.1032,L4f/MAX(L0f)
oplot,lam-0.172 ,L5f/MAX(L0f)
oplot,wave+0.1720+0.05,line[0:300]    -0.005,col=180,thick=2
oplot,wave+0.1032+0.05,line[301:601]  -0.005,col=180,thick=2
oplot,wave+0.0344+0.05,line[602:902]  -0.005,col=180,thick=2
oplot,wave-0.0344+0.05,line[903:1203] -0.005,col=180,thick=2
oplot,wave-0.1032+0.05,line[1204:1504]-0.005,col=180,thick=2
oplot,wave-0.1720+0.05,line[1505:1805]-0.005,col=180,thick=2


oplot,wave+0.1720,line[0:300]    -0.01,col=180,thick=2
oplot,wave+0.1032,line[301:601]  -0.01,col=180,thick=2
oplot,wave+0.0344,line[602:902]  -0.01,col=180,thick=2
oplot,wave-0.0344,line[903:1203] -0.01,col=180,thick=2
oplot,wave-0.1032,line[1204:1504]-0.01,col=180,thick=2
oplot,wave-0.1720,line[1505:1805]-0.01,col=180,thick=2


SAVE,L0b,L1b,L2b,L3b,L4b,L5b,L0c,L1c,L2c,L3c,L4c,L5c,L0d,L1d,L2d,L3d,L4d,L5d,L0e,L1e,L2e,L3e,L4e,L5e,L0f,L1f,L2f,L3f,L4f,L5f,file='lineatdiskcenter.bin'


END
