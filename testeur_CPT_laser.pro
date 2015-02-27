; PROGRAM TO PRODUCE ARTIFICIAL DATA TO TEST CPT_laser.pro
; AND SEE WHAT EFFECT WRONG FSRs HAVE

PRO testeur_CPT_laser

lamref      = 6173.3433d0
nx          = 256              ; number of rows
ny          = 256              ; number of columns

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"
nseq        = 27
nlam        = 40000            ; number of wavelengths
dlam        = 1.d0/2.d3       ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
dpi         = 2.d0*!dpi
seed        = 1l
Inten       = DBLARR(nx,ny,nseq)

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172d0;-0.001057600d0 ; for the narrow-band Michelson
FSR[1]      = 0.344d0;-0.002076830d0 ; for the broad-band  Michelson
FSR[2]      = 0.690d0;+0.000483467d0 ; for E1
FSR[3]      = 1.407d0          ; for E2
FSR[4]      = 2.779d0          ; for E3
FSR[5]      = 5.682d0          ; for E4
FSR[6]      = 11.354d0         ; for E5

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' 
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q           = READFITS('blocker11.fits')      ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lamref,lam);,/LSQUADRATIC) ; I center the profile


; TUNING POSITIONS
tuning      = DBLARR(3,nseq)
tuning[*,0] = [         0.d0,          0.d0,          0.d0]
tuning[*,1] = [        80.d0,          0.d0,          0.d0]
tuning[*,2] = [       160.d0,          0.d0,          0.d0]
tuning[*,3] = [         0.d0,         80.d0,          0.d0]
tuning[*,4] = [        80.d0,         80.d0,          0.d0]
tuning[*,5] = [       160.d0,         80.d0,          0.d0]
tuning[*,6] = [         0.d0,        160.d0,          0.d0]
tuning[*,7] = [        80.d0,        160.d0,          0.d0]
tuning[*,8] = [       160.d0,        160.d0,          0.d0]
tuning[*,9] = [         0.d0,          0.d0,         80.d0]
tuning[*,10]= [        80.d0,          0.d0,         80.d0]
tuning[*,11]= [       160.d0,          0.d0,         80.d0]
tuning[*,12]= [         0.d0,         80.d0,         80.d0] 
tuning[*,13]= [        80.d0,         80.d0,         80.d0]
tuning[*,14]= [       160.d0,         80.d0,         80.d0]
tuning[*,15]= [         0.d0,        160.d0,         80.d0]
tuning[*,16]= [        80.d0,        160.d0,         80.d0]
tuning[*,17]= [       160.d0,        160.d0,         80.d0]
tuning[*,18]= [         0.d0,          0.d0,        160.d0]
tuning[*,19]= [        80.d0,          0.d0,        160.d0]
tuning[*,20]= [       160.d0,          0.d0,        160.d0]
tuning[*,21]= [         0.d0,         80.d0,        160.d0]
tuning[*,22]= [        80.d0,         80.d0,        160.d0]
tuning[*,23]= [       160.d0,         80.d0,        160.d0]
tuning[*,24]= [         0.d0,        160.d0,        160.d0]
tuning[*,25]= [        80.d0,        160.d0,        160.d0]
tuning[*,26]= [       160.d0,        160.d0,        160.d0]


FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [6.d0,6.d0,-6.d0]*!dpi/180.d0


; DATA

Phig    = FLTARR(nx,ny,7)
FOR i=0,6 DO phig[*,*,i] = SHIFT(DIST(nx,nx),nx/2,nx/2)/FLOAT(nx)*!dpi/2.d0+i*10.d0/180.d0*!dpi;(RANDOMU(seed,nx,ny,7)*160.-80.)*!dpi/180.d0
Bg      = FLTARR(nx,ny,7);RANDOMU(seed,nx,ny,7)*0.1+0.9
FOR i=0,6 DO Bg[*,*,i] = SHIFT(DIST(nx,nx),nx/2,nx/2)/3.5/FLOAT(nx)+0.8d0
I0g     = FLTARR(nx,ny)+9500.;RANDOMU(seed,nx,ny)*500.+9500.

; ARTIFICIAL INTENSITIES

lam0    = 6173.3232d0-lamref;6175.23d0-lamref

FOR iii=0,nx-1 DO BEGIN
    FOR jjj=0,ny-1 DO BEGIN
        FOR j=0,nseq-1 DO BEGIN
            profileg                 = INTERPOL(blocker0,lam,lam0)
            profileg                 = profileg * I0g[iii,jjj]
            FOR i=0,2 DO profileg    = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0+Phig[iii,jjj,i]+tuning[i,j]))/2.d0
            FOR i=3,6 DO profileg    = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0+Phig[iii,jjj,i]))/2.d0
            Inten[iii,jjj,j] = profileg
        ENDFOR
    ENDFOR        
ENDFOR


SAVE,I0g,Bg,Phig,lam0,FILE='testeur1.bin'
SAVE,Inten,FILE='testeur2.bin'

xcenter  = nx/2             ; center of the solar disk, depends on the image !
ycenter  = nx/2
distance = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center
anglim   = 960.d0

; WE APPLY A MASK TO KEEP ONLY THE PHASES AND CONTRASTS ON THE ELEMENT
; APERTURE
Intenmap    = TOTAL(Inten[*,*,*],3)
Intenmapmax = MAX(Intenmap)

a = WHERE(REFORM(Bg[*,*,0]) EQ 0.0 OR distance[*,*] GT anglim)

FOR i=0,2 DO BEGIN
    temp=REFORM(Phig[*,*,i])
    temp[a]=-10000.d0
    phig[*,*,i]=temp[*,*]
    temp=REFORM(Bg[*,*,i])
    temp[a]=-10000.d0
    Bg[*,*,i]=temp[*,*]
ENDFOR

a = WHERE(FINITE(Phig) EQ 0)
IF(a[0] NE -1) THEN Phig[a]= -10000.0
a = WHERE(FINITE(Bg) EQ 0)
IF(a[0] NE -1) THEN Bg[a]  = -10000.0

a = WHERE(REFORM(Bg[*,*,0]) EQ -10000.d0,COMPLEMENT=ba)

; ADD 180 DEGREES TO PHASES < -180
FOR i=0,2 DO BEGIN
temp = REFORM(Phig[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi
        Phig[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+dpi
        Phig[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy   = DBLARR(3)
moyb  = moy
mini  = DBLARR(3)
maxi  = DBLARR(3)
mini2 = DBLARR(3)
maxi2 = DBLARR(3)
FOR i=0,2 DO BEGIN
    temp    = REFORM(Phig[*,*,i])
    a       = WHERE(temp NE -10000.d0,COMPLEMENT=b)
    moy[i]  = MEAN(temp[a])
    mini[i] = MIN(temp[a])*180./!dpi
    maxi[i] = MAX(temp[a])*180./!dpi
    temp[b] = moy[i]
    Phig[*,*,i]=temp
    temp    = REFORM(Bg[*,*,i])
    a       = WHERE(temp NE -10000.d0,COMPLEMENT=b)
    moyb[i] = MEAN(temp[a])
    temp[b] = moyb[i]
    mini2[i]= MIN(temp[a])
    maxi2[i]= MAX(temp[a])
    Bg[*,*,i]=temp
ENDFOR

; WE PLOT THE RESULT
SET_PLOT,'ps'
!p.multi=[0,2,3]
 DEVICE,file='temp.ps',bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27,/color

loadct,4

TVIM,phig[*,*,0]*180.d0/!dpi,/scale,tit='!17NB Michelson, !7l!17='+STRTRIM(STRING(moy[0]*180./!dpi),1),xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[0],maxi[0]],pcharsize=1.5,barwidth=0.5

TVIM,phig[*,*,1]*180.d0/!dpi,/scale,tit='!17WB Michelson, !7l!17='+STRTRIM(STRING(moy[1]*180./!dpi),1),xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[1],maxi[1]],pcharsize=1.5,barwidth=0.5

TVIM,phig[*,*,2]*180.d0/!dpi,/scale,tit='!17E1, !7l!17='+STRTRIM(STRING(moy[2]*180./!dpi),1),xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[2],maxi[2]],pcharsize=1.5,barwidth=0.5

TVIM,Bg[*,*,0],/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',pcharsize=1.5,barwidth=0.5,range=[0.85,1.05];,range=[mini2[0],maxi2[0]]
temp = Bg[*,*,0]
hist = histogram(temp[ba],binsize=0.001,min=0.7,max=1.05)
PLOT,FINDGEN(N_ELEMENTS(hist))*0.001+0.7,hist/FLOAT(N_ELEMENTS(ba)),psym=10,tit='!7r!17='+STRING(SIGMA(temp[ba]))+'/!7l!17='+STRING(MEAN(temp[ba])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5,xst=1

TVIM,Bg[*,*,1],/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',pcharsize=1.5,barwidth=0.5,range=[0.85,1.05];,range=[mini2[1],maxi2[1]]
temp = Bg[*,*,1]
hist = histogram(temp[ba],binsize=0.001,min=0.7,max=1.05)
PLOT,FINDGEN(N_ELEMENTS(hist))*0.001+0.7,hist/FLOAT(N_ELEMENTS(ba)),psym=10,tit='!7r!17='+STRING(SIGMA(temp[ba]))+'/!7l!17='+STRING(MEAN(temp[ba])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5,xst=1

TVIM,Bg[*,*,2],/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',pcharsize=1.5,barwidth=0.5,range=[0.85,1.05];,range=[mini2[2],maxi2[2]]
temp = Bg[*,*,2]
hist = histogram(temp[ba],binsize=0.001,min=0.7,max=1.05)
PLOT,FINDGEN(N_ELEMENTS(hist))*0.001+0.7,hist/FLOAT(N_ELEMENTS(ba)),psym=10,tit='!7r!17='+STRING(SIGMA(temp[ba]))+'/!7l!17='+STRING(MEAN(temp[ba])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5,xst=1


TVIM,Intenmap/Intenmapmax,/scale,tit='!17Laser Illumination Estimate',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,pcharsize=1.5,stit='Relative Intensity'
TVIM,I0g/MAX(I0g),/scale,tit='!17Actual Laser Illumination',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,pcharsize=1.5,stit='Relative Intensity'

DEVICE,/close
PRINT,'LASER WAVELENGTH:',lam0+lamref
PRINT,'AVERAGE PHASES (FOR NARROW-BAND MICHELSON, BROAD-BAND MICHELSON, AND E1)'
FOR i=0,2 DO PRINT,moy[i]*180.d0/!dpi
PRINT,'AVERAGE CONTRASTS'
FOR i=0,2 DO PRINT,moyb[i]

READ,pause

END
