; THIS PROGRAM ANALYZE THE DOPPLERGRAMS
; TAKEN BY THE LINDOP SEQUENCE

PRO HMI_lindop,draw


!p.multi = 0

dpi         = 2.d0*!dpi
lamref      = 6173.3433        ; Fe I 6173 line central wavelength in air
nseq        = 5                ; number of positions in the co-tune sequence
nlam        = 30000            ; number of wavelengths
dlam        = 1.d0/2.d3       ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
nx          = 256
ndop        = 17
lam0g       = FLTARR(nx,nx,ndop)

base='/tmp20/schou/hmi060301b/'
list='listDop_060301_170309'

IF draw EQ 1 THEN GOTO,draw

; FIRST WE HAVE A LOOK AT WHAT THE FILTERS LOOK LIKE
;-------------------------------------------------------

; CO-TUNE SEQUENCE FOR THE SUNLIGHT TESTS
tuningJS     = DBLARR(nseq,3)
tuningJS[0,*]= [  129 , 124 ,   55]
tuningJS[1,*]= [  135 , 112 ,   91]
tuningJS[2,*]= [   81 , 100 ,   67]
tuningJS[3,*]= [   87 ,  88 ,   43]
tuningJS[4,*]= [   93,   76 ,   79]
tuningJS     = REVERSE(tuningJS,2)
FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]

; PHASE MAPS
phase    = DBLARR(3)
phase[0] = -55.75    ; NB Michelson
phase[1] = 130.80    ; WB Michelson
phase[2] = 106.14    ; E1 
phase    = phase*!dpi/180.d0

; COMPUTATION OF THE FILTER PROFILES

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom
FSR[0]      = 0.172457         ; for the narrow-band Michelson
FSR[1]      = 0.344242         ; for the broad-band  Michelson
FSR[2]      = 0.702            ; for E1
FSR[3]      = 1.405            ; for E2
FSR[4]      = 2.779            ; for E3
FSR[5]      = 5.682            ; for E4
FSR[6]      = 11.354           ; for E5

RESTORE,'frontwindow.bin' 
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam)
q           = READFITS('blocker11.fits')
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam)

profilef = DBLARR(nlam)+1.d0
profilef  = blocker0
FOR i=0,3 DO profilef  = profilef * 0.5d0 * (1.d0+COS(dpi/FSR[i+3]*lam))
profile  = DBLARR(nseq,nlam)

FOR j=0,nseq-1 DO BEGIN
    profile[j,*] = profilef
    FOR i=0,2 DO profile[j,*] = profile[j,*] * 0.5d0 * (1.d0+COS(dpi/FSR[i]*lam+phase[i]+tuningJS[j,i]))
ENDFOR


; WE PLOT THE FILTERS
WINDOW,0,RETAIN=2
plot,lam,profile[0,*],thick=3,xst=1,tit='!17',charsize=1.5,xrange=[-.6,.6],yrange=[0,0.7]
oplot,lam,profile[1,*],linestyle=2
oplot,lam,profile[2,*],linestyle=3
oplot,lam,profile[3,*],linestyle=4
oplot,lam,profile[4,*],linestyle=2,thick=2
oplot,lam,blocker0

READ,PAUSE

imx = DBLARR(nx,nx,ndop*nseq)
dl  = 0.069d0

cste = 2.99792458d8/6173.3433d0/(2.d0*!dpi/nseq/dl)



;header       = STRARR(N,91)     ; 90 is the number of keywords in the filtergram headers
;imx[*,*,*] = READFILE_SCAN(header,nx,base,list) ; SUBROUTINE WRITTEN BY S. COUVIDAT
RESTORE,'SEQUENCE_LINDOP'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

FOR ii=0,ndop-1 DO BEGIN
    PRINT,ii
    FOR jjj=0,nx-1 DO BEGIN
       ;PRINT,jjj
        FOR iii=0,nx-1 DO BEGIN

            f1 = 0.d0
            FOR i=0,nseq-1 DO f1 = f1 + exp(DCOMPLEX(0,1)*2.d0*!dpi*i/DOUBLE(nseq))*imx[iii,jjj,i+ii*nseq]
            lam0g[iii,jjj,ii] = ATAN(IMAGINARY(f1),DOUBLE(f1))*cste
            
        ENDFOR
    ENDFOR

ENDFOR

center=[130,128]
distance = SHIFT(DIST(nx,nx),center[0],center[1])*0.5d0*4096.d0/nx
a = WHERE(distance LT 935.d0,COMPLEMENT=b) ;!!!WARNING!!!
distance[a]=1.d0
distance[b]=0.d0
FOR i=0,ndop-1 DO lam0g[*,*,i] = lam0g[*,*,i] * distance


SAVE,lam0g,FILE='RESULTS/RESULTS_LINDOP_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

draw:
RESTORE,'RESULTS/RESULTS_LINDOP_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

power=ABS(FFT(lam0g))^2.0
puissance=FLTARR(nx/2-1,ndop)
distance=DIST(nx,nx)
FOR i=0,nx/2-2 DO BEGIN
a=where(distance GT i AND distance LE i+1)
FOR j=0,ndop-1 DO puissance[i,j]=TOTAL(power[*,*,j]*distance[a])
ENDFOR

WINDOW,0,RETAIN=2
tvim,puissance,aspect=1

READ,pause

END
