; PROGRAM TO DO THE CONE-AVERAGING OF THE MICHELSON PHASE MAPS

FUNCTION phase,ii

IF(ii EQ 0) THEN BEGIN
    OPENR,1,'Michelson1a.txt'   ; phase map for the WB Michelson
    x=FLTARR(30)                ; coordinates of the center of the tile
    y=x                         ; that is why the diameter of the image
    READF,1,x                   ; is larger than x[0]+x[N_ELEMENTS(x)-1]
    phase=FLTARR(31,30)
    READF,1,phase
    y=phase[0,*]
    phase=phase[1:30,*]
    CLOSE,1
    a=WHERE(phase EQ 0.0,COMPLEMENT=b)
    min1=MIN(phase[b])
    max1=MAX(phase[b])
    av1=MEAN(phase[b])
    phase[b]=phase[b]-av1
    min1=min1-av1
    max1=max1-av1
    DEFSYSV,'!dx1',x[1]-x[0]
    DEFSYSV,'!x1',x[0]
ENDIF ELSE BEGIN
    OPENR,1,'Michelson2a.txt'   ; phase map for the NB Michelson
    x=FLTARR(32)
    y=x
    READF,1,x
    phase=FLTARR(33,32)
    READF,1,phase
    y=phase[0,*]
    phase=phase[1:32,*]
    CLOSE,1
    a=WHERE(phase EQ 0.0,COMPLEMENT=b)
    min2=MIN(phase[b])
    max2=MAX(phase[b])
    av2=MEAN(phase[b])
    phase[b]=phase[b]-av2
    min2=min2-av2
    max2=max2-av2
    DEFSYSV,'!dx1',x[1]-x[0]
    DEFSYSV,'!x1',x[0]
ENDELSE

RETURN,phase
END


FUNCTION contrast,ii

IF(ii EQ 0) THEN BEGIN
    OPENR,1,'Michelson1b.txt'   ; contrast map for the WB Michelson
    x=FLTARR(30)
    y=x
    READF,1,x
    contrast=FLTARR(31,30)
    READF,1,contrast
    y=contrast[0,*]
    contrast=contrast[1:30,*]
    CLOSE,1
    DEFSYSV,'!dx2',x[1]-x[0]
    DEFSYSV,'!x2',x[0]
ENDIF ELSE BEGIN
    OPENR,1,'Michelson2b.txt'   ; contrast map for the NB Michelson
    x=FLTARR(31)
    y=x
    READF,1,x
    contrast=FLTARR(32,31)
    READF,1,contrast
    y=contrast[0,*]
    contrast=contrast[1:31,*]
    CLOSE,1
    DEFSYSV,'!dx2',x[1]-x[0]
    DEFSYSV,'!x2',x[0]
ENDELSE

RETURN,contrast
END


;----------------------------------------------------------------------
;
; MAIN PROGRAM
; this code computes the phase and contrast maps
; for the Michelsons simulating a cone averaging
;
; For solar radius of 976 arcseconds
; the image sizes are:
; 9.833 mm (radius) for the Broad-band Michelson
; 5.806 mm for the Narrow-Band Michelson
; the sizes of the light cones are:
; 4.63 mm (radius) for the Broad-band Michelson
; 8.69 mm for the Narrow-band Michelson
; DIAMETER OF DICK SHINE'S IMAGES = 33.7 mm (SO XSIZE=39.466 MM) (Dick
; Shine, e-mail from March 10, 2006
;----------------------------------------------------------------------


PRO Michelsoncalib,draw

;set_plot,'x'
;!p.multi=[0,1,3];4]
;window,0,retain=2,ysize=800

largeur   = FLTARR(2); provided by Jesper Schou
largeur[0]= 4.63    ; radius (in mm) of the beam on the Wide Band Michelson aperture
largeur[1]= 8.69    ; radius (in mm) of the beam on the Narrow Band Michelson aperture
nx1       = 325;215 ;860.
length    = 39.4685;34.9 ; image size in mm
dx        = length/nx1
phase     = FLTARR(2,nx1,nx1)

SET_PLOT,'ps'

IF draw EQ 1 THEN GOTO,draw

!p.multi=[0,1,3]                ;4]
device,file='yo.ps',bits=24,xsize=20,xoffset=0,ysize=26,yoffset=0.5,/color
loadct,3

FOR ii=0,1 DO BEGIN  ; ii=0 FOR WB Michelson, ii=1 FOR NB Michelson

   ;phase1    = PHASE(ii)   *!pi/180.0
   ;contrast1 = CONTRAST(ii)


    ;nx1=N_ELEMENTS(contrast1[*,0]) ; for ii=1 the diameters of the phase and contrast maps are not the same
    ;IF(ii EQ 1) THEN BEGIN
    ;    IX=FINDGEN(nx1)*!dx2+!x2-!x1
    ;    IY=FINDGEN(nx1)*!dx2+!x2-!x1
    ;    phase2=BILINEAR(phase1,IX,IY)
    ;ENDIF


    IF(ii EQ 1) THEN BEGIN
       ;phase1=READFITS('michN00_18jul2005.fits') 
        phase1=READFITS('michNB_28feb2006.fits')
       ;phase1=READFITS('michNold_9mar2006.fits')
       ;phase1 = REVERSE(phase1,1)
       ;phase1 = ROT(phase1,90)
    ENDIF ELSE BEGIN 
       ;phase1=READFITS('michBE00_20jul2005.fits')
        phase1=READFITS('michWB_24feb2006.fits')
       ;phase1=READFITS('michWold_9mar2006.fits')
       ;phase1 = ROT(phase1,180)
    ENDELSE
    phase1 = phase1 * !dpi/180.d0
    phase1 = REBIN(phase1,nx1,nx1)

   
    IF( ii EQ 0) THEN ny1=ROUND(largeur[0]*2./dx) ELSE ny1=ROUND(largeur[1]*2/dx);300
    disc1=SHIFT(DIST(ny1,ny1)*dx,ny1/2,ny1/2)
    ;ROUND(largeur[ii]/!dx1*2.0)
    ;disc1=SHIFT(DIST(ny1,ny1)*!dx1,ny1/2,ny1/2)
    a=WHERE(disc1 le largeur[ii],COMPLEMENT=b)
    disc1[a]=1.0
    na=FLOAT(N_ELEMENTS(a))
    disc1[b]=0.0

    resultat = FLTARR(2,nx1,nx1) ; average phase

; the transmission profile is (1+contrast cos(---+2 phase))/2
    FOR i=ny1/2,nx1-1-ny1/2 DO BEGIN
        PRINT,i
        FOR j=ny1/2,nx1-1-ny1/2 DO BEGIN
          
            IF(ny1 MOD 2 EQ 0) THEN BEGIN
                discp=disc1*phase1[i-ny1/2:i+ny1/2-1,j-ny1/2:j+ny1/2-1] ;phase
                ;discc=disc1*contrast1[i-ny1/2:i+ny1/2-1,j-ny1/2:j+ny1/2-1] ;contrast
            ENDIF ELSE BEGIN
                discp=disc1*phase1[i-ny1/2:i+ny1/2,j-ny1/2:j+ny1/2] ;phase
                ;discc=disc1*contrast1[i-ny1/2:i+ny1/2,j-ny1/2:j+ny1/2] ;contrast
            ENDELSE               


           ;youp=TOTAL(discc[a]*exp(COMPLEX(0,1)*discp[a]))
            youp=TOTAL(exp(COMPLEX(0,1)*discp[a]))
           ;resultat[0,i,j]=1./na*SQRT(youp*TOTAL(discc[a]*exp(-COMPLEX(0,1)*discp[a]))) ; contrast
           resultat[0,i,j]=1./na*SQRT(youp*TOTAL(exp(-COMPLEX(0,1)*discp[a]))) ; contrast
           IF(resultat[0,i,j] NE 0.0) THEN resultat[1,i,j]=ALOG(youp/resultat[0,i,j]/na)/COMPLEX(0,1)  ELSE resultat[1,i,j]=0.0 ; phase
        ENDFOR
       ;tvim,phase1,/scale,tit='!17',xtit='!17x',ytit='!17y',stit='!17Phase (in degrees)'
       ;tvim,contrast1,/scale,tit='!17',xtit='!17x',ytit='!17y',stit='!17Contrast'
       ;tvim,resultat[1,*,*]*180.0/!pi,/scale
       ;tvim,resultat[0,*,*],/scale
    ENDFOR
    
    disc=dist(nx1*2,nx1*2)*dx/2.;!dx1/2.
    disc=shift(disc,nx1,nx1)
    a=WHERE(FINDGEN(nx1*2) MOD 2 EQ 1)
    disc2=FLTARR(nx1,nx1)
    FOR iii=0,N_ELEMENTS(a)-1 DO disc2[iii,*]=disc[a[iii],a] ; because there is no (x,y)=(0,0)

    IF(ii EQ 0) THEN a=where(disc2 LE 9.833,COMPLEMENT=b) ; 10.4 mm is the actual image radius for the WB Michelson
                                                         ; I use 10.5 because it is a value of x
    IF(ii EQ 1) THEN a=where(disc2 LE 5.806,COMPLEMENT=b)  ;  6.1 mm is the image radius for the NB Michelson
    disc2[a]=1.0
    disc2[b]=0.0
    resultat[0,*,*]=resultat[0,*,*]*disc2
    resultat[1,*,*]=resultat[1,*,*]*disc2
    SAVE,resultat,file='M_Shine_'+STRTRIM(STRING(ii),1)+'.bin'
    
    PRINT,'minimum/max phase',MIN(resultat[1,*,*])*180./!pi,MAX(resultat[1,*,*])*180./!pi
    temp=REFORM(resultat[0,*,*])
    PRINT,'minimum contrast',MIN(temp[a])
    temp[WHERE(temp EQ 0.0)]=2.0
    temp2=REFORM(resultat[1,*,*])
    temp2[WHERE(temp2 EQ 0.0)]=MAX(resultat[1,*,*])+5.0*!pi/180.

    ;a=WHERE(contrast1 NE 0.0,COMPLEMENT=b)
    ;contrast1[b]=1.0
    

    tvim,phase1*180./!pi,/scale,tit='!17Element M'+STRTRIM(STRING(ii+1),1)+': phase map',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length];,xrange=[!x1-!dx1/2.,-!x1+!dx1/2.],yrange=[!x1-!dx1/2.,-!x1+!dx1/2.]
    ;tvim,contrast1,/scale,tit='!17Element M'+STRTRIM(STRING(ii+1),1)+': contrast map',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Contrast',range=[MIN(contrast1[a]),MAX(contrast1[a])],pcharsize=1.5,xrange=[!x1-!dx1/2.,-!x1+!dx1/2.],yrange=[!x1-!dx1/2.,-!x1+!dx1/2.]
    tvim,temp2*180./!pi,/scale,tit='!17Transmission I(!7k!17)=(1+Contrast cos(2!7pk!17/FSR+ Phase))/2',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length];,xrange=[!x1-!dx1/2.,-!x1+!dx1/2.],yrange=[!x1-!dx1/2.,-!x1+!dx1/2.]
    tvim,temp,/scale,range=[MIN(temp),MAX(resultat[0,*,*])],tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Contrast',pcharsize=1.5,xrange=[0,length],yrange=[0,length];,xrange=[!x1-!dx1/2.,-!x1+!dx1/2.],yrange=[!x1-!dx1/2.,-!x1+!dx1/2.]
    
    phase[ii,*,*] = phase1*180./!pi

ENDFOR    

read,pause
SAVE,phase,FILE='Michelsons_phase.bin'

DEVICE,/close


;--------------------------------------------------------------------------------------------------------
; FOR ANGULAR DEPENDENCE TEST
;--------------------------------------------------------------------------------------------------------

draw:

!p.multi=[0,1,3]
DEVICE,file='yo2.ps',bits=24,xsize=26,xoffset=-2,ysize=27,yoffset=0,/color

; BROAD-BAND MICHELSON

phase1=READFITS('michWB_24feb2006.fits')
;phase1=READFITS('michBE00_20jul2005.fits')
phase1 = REBIN(phase1,nx1,nx1)
;phase1 = ROT(phase1,180)

xshift0 = 3.63786/dx ; 3.63786 = 55.*4.63/70. ; 70=calmode image radius, 55 mm is where the holes are located in radial distance, 4.63 is the light cone radius
yshift0 = 0
cercle0 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift0,nx1/2+yshift0)
a = WHERE(ABS(cercle0-9.833)/9.833 LT 0.01,COMPLEMENT=b)
cercle0[a]=2000*(phase1[a])
cercle0[b]=1.0

xshift1 = -3.63786/dx
yshift1 = 0
cercle1 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift1,nx1/2+yshift1)
a = WHERE(ABS(cercle1-9.833)/9.833 LT 0.01,COMPLEMENT=b)
cercle1[a]=2000*(phase1[a])
cercle1[b]=1.0

tvim,phase1[*,*]*cercle0*cercle1,/scale,tit='!17Broad-Band Michelson: phase map',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length],range=[-24,24]


yshift2 = -3.63786/dx
xshift2 = 0
cercle2 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift2,nx1/2+yshift2)
a = WHERE(ABS(cercle2-9.833)/9.833 LT 0.01,COMPLEMENT=b)
cercle2[a]=2000*(phase1[a])
cercle2[b]=1.0

yshift3 = 3.63786/dx
xshift3 = 0
cercle3 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift3,nx1/2+yshift3)
a = WHERE(ABS(cercle3-9.833)/9.833 LT 0.01,COMPLEMENT=b)
cercle3[a]=2000*(phase1[a])
cercle3[b]=1.0

tvim,phase1*cercle2*cercle3,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length],range=[-24,24]

yshift4 = 0
xshift4 = 0
cercle4 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift4,nx1/2+yshift4)
a = WHERE(ABS(cercle4-9.833)/9.833 LT 0.01,COMPLEMENT=b)
cercle4[a]=2000*(phase1[a])
cercle4[b]=1.0

tvim,phase1*cercle4,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length],range=[-24,24]


; NARROW BAND MICHELSON

;phase1=READFITS('michN00_18jul2005.fits')
phase1=READFITS('michNB_28feb2006.fits')
phase1 = REBIN(phase1,nx1,nx1)
;phase1 = REVERSE(phase1,1)
;phase1 = ROT(phase1,90)

xshift0 = 6.82786/dx ; 8.69*55./70.
yshift0 = 0
cercle0 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift0,nx1/2+yshift0)
a = WHERE(ABS(cercle0-5.806)/5.806 LT 0.01,COMPLEMENT=b)
cercle0[a]=2000*(phase1[a])
cercle0[b]=1.0

xshift1 = -6.82786/dx
yshift1 = 0
cercle1 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift1,nx1/2+yshift1)
a = WHERE(ABS(cercle1-5.806)/5.806 LT 0.01,COMPLEMENT=b)
cercle1[a]=2000*(phase1[a])
cercle1[b]=1.0

tvim,phase1[*,*]*cercle0*cercle1,/scale,tit='!17Narrow-Band Michelson: phase map',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length],range=[-42,49]

yshift2 = -6.82786/dx
xshift2 = 0
cercle2 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift2,nx1/2+yshift2)
a = WHERE(ABS(cercle2-5.806)/5.806 LT 0.01,COMPLEMENT=b)
cercle2[a]=2000*(phase1[a])
cercle2[b]=1.0

yshift3 = 6.82786/dx
xshift3 = 0
cercle3 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift3,nx1/2+yshift3)
a = WHERE(ABS(cercle3-5.806)/5.806 LT 0.01,COMPLEMENT=b)
cercle3[a]=2000*(phase1[a])
cercle3[b]=1.0

tvim,phase1*cercle2*cercle3,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length],range=[-42,49]

yshift4 = 0
xshift4 = 0
cercle4 = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift4,nx1/2+yshift4)
a = WHERE(ABS(cercle4-5.806)/5.806 LT 0.01,COMPLEMENT=b)
cercle4[a]=2000*(phase1[a])
cercle4[b]=1.0

tvim,phase1*cercle4,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length],range=[-42,49]


DEVICE,/close

phase1=READFITS('michN00_18jul2005.fits')
phase1 = REBIN(phase1,215,215)

phase1 = REVERSE(phase1,1)
;phase1 = REVERSE(phase1,2)
phase1 = ROT(phase1,90);100
;phase1 = REVERSE(phase1,1)
;phase1 = REVERSE(phase1,2)

yshift = [-6.82786,0,6.82786,0,0]/dx
xshift = [0,0,0,-6.82786,6.82786]/dx

DEVICE,file='yo3.ps',bits=24,xsize=26,xoffset=-2,ysize=27,yoffset=0,/color
!p.multi=[0,2,3]

tvim,phase1,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length]
FOR i=0,4 DO BEGIN
cercle = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift[i],nx1/2+yshift[i])
a = WHERE(cercle LE 5.806,COMPLEMENT=b)
cercle[a]=1.0
cercle[b]=0.0

tvim,phase1*cercle,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length]
ENDFOR

yshift = [-3.63786,0,3.63786,0,0]/dx
xshift = [0,0,0,-3.63786,3.63786]/dx

phase1=READFITS('michBE00_20jul2005.fits')
phase1 = REBIN(phase1,215,215)
phase1 = ROT(phase1,180);145)

tvim,phase1,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length]
FOR i=0,4 DO BEGIN
cercle = SHIFT(DIST(nx1,nx1)*dx,nx1/2+xshift[i],nx1/2+yshift[i])
a = WHERE(cercle LE 9.833,COMPLEMENT=b)
cercle[a]=1.0
cercle[b]=0.0

tvim,phase1*cercle,/scale,tit='!17',xtit='!17x (mm)',ytit='!17y (mm)',stit='!17Phase (in degrees)',pcharsize=1.5,xrange=[0,length],yrange=[0,length]
ENDFOR

DEVICE,/CLOSE

set_plot,'x'


END
