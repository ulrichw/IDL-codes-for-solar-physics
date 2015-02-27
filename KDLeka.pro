; program to analyze the results of the DRMS module HMIfilters
; using the files by K.D. Leka

PRO KDLeka

nx=512

d=FLTARR(24,nx,nx)

d[0,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/I0.fits"),nx,nx)
d[1,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/I1.fits"),nx,nx)
d[2,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/I2.fits"),nx,nx)
d[3,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/I3.fits"),nx,nx)
d[4,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/I4.fits"),nx,nx)
d[5,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/I5.fits"),nx,nx)
d[6,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/Q0.fits"),nx,nx)
d[7,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/Q1.fits"),nx,nx)
d[8,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/Q2.fits"),nx,nx)
d[9,*,*] =REBIN(readfits("/SUM0/D21866037/S00000/Q3.fits"),nx,nx)
d[10,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/Q4.fits"),nx,nx)
d[11,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/Q5.fits"),nx,nx)
d[12,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/U0.fits"),nx,nx)
d[13,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/U1.fits"),nx,nx)
d[14,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/U2.fits"),nx,nx)
d[15,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/U3.fits"),nx,nx)
d[16,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/U4.fits"),nx,nx)
d[17,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/U5.fits"),nx,nx)
d[18,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/V0.fits"),nx,nx)
d[19,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/V1.fits"),nx,nx)
d[20,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/V2.fits"),nx,nx)
d[21,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/V3.fits"),nx,nx)
d[22,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/V4.fits"),nx,nx)
d[23,*,*]=REBIN(readfits("/SUM0/D21866037/S00000/V5.fits"),nx,nx)

!p.multi=[0,2,3]
loadct,3
set_plot,'ps'
device,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=24,/color,bits=24
for i=0,23 do tvim,d[i,*,*],/scale,barwidth=0.5
device,/close

REAd,pause


END
