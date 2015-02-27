; THIS PROGRAM RECONSTRUCTS A DETUNE27 SEQUENCE FROM THE
; RESULTS OBTAINED WITH HMI_sun1.pro, TO CHECK THE
; ACCURACY OF THESE RESULTS

PRO reconstruction_calmode,draw

nx      = 128
list    ='listSun060224_225806'

nseq    = 27
nlam    = 2250               ; number of wavelengths
dlam    = 3.6d0/1.d3         ; resolution in wavelength
lam     = dlam*(DINDGEN(nlam)-(nlam-1)/2.)
lamref  = 6173.3433d0
dlamdv  = lamref/2.99792458d8 

IF draw EQ 1 THEN GOTO,draw

; RELATIVE-PHASE MAPS
;------------------------------------------------------------

;RESTORE,'RESULTS/RESULTS_listSun060224_225806_256.BIN'; phases in calmode
;RESTORE,'RESULTS/RESULTS_listSun060224_225806_128.BIN' ; idem but more recent and better dark removal, and contrast=0.98 for everything
RESTORE,'RESULTS/RESULTS_listSun060616_205620_128.BIN'
RESTORE,'RESULTS/RESULTS_listLaser060620_232532_128.BIN'
lam0=374.27737*dlamdv
;lam0     = 633.691*dlamdv  ; !!!!!!!!!!!!!!!!!!!
phig=phig0-180.d0
linewidth[*,*]=0.071
linedepth[*,*]=continuum[*,*]*0.51



Phig     = REBIN(Phig,nx,nx,3)*!dpi/180.d0
LINEDEPTH= REBIN(linedepth,nx,nx,3)
LINEWIDTH= REBIN(linewidth,nx,nx,3)
CONTINUUM= REBIN(continuum,nx,nx,3)
contrast = 0.98d0; contrasts assumed to obtained the phase map wiht HMI_sun1.pro
phase    = 0.d0  ; phases assumed for the non-tunable elements

; DEFINITION OF THE DETUNE27 SEQUENCE
;------------------------------------------------------------

tuning        = DBLARR(3,nseq)
tuning[*, 0]  = [         0.d0,          0.d0 ,         0.d0]
tuning[*, 1]  = [        80.d0,          0.d0 ,         0.d0]
tuning[*, 2]  = [       160.d0,          0.d0 ,         0.d0]
tuning[*, 3]  = [         0.d0,         80.d0,          0.d0]
tuning[*, 4]  = [        80.d0,         80.d0,          0.d0]
tuning[*, 5]  = [       160.d0,         80.d0,          0.d0]
tuning[*, 6]  = [         0.d0,        160.d0,          0.d0]
tuning[*, 7]  = [        80.d0,        160.d0,          0.d0]
tuning[*, 8]  = [       160.d0,        160.d0,          0.d0]
tuning[*, 9]  = [         0.d0,          0.d0,         80.d0]
tuning[*,10]  = [        80.d0,          0.d0,         80.d0]
tuning[*,11]  = [       160.d0,          0.d0,         80.d0]
tuning[*,12]  = [         0.d0,         80.d0,         80.d0] 
tuning[*,13]  = [        80.d0,         80.d0,         80.d0]
tuning[*,14]  = [       160.d0,         80.d0,         80.d0]
tuning[*,15]  = [         0.d0,        160.d0,         80.d0]
tuning[*,16]  = [        80.d0,        160.d0,         80.d0]
tuning[*,17]  = [       160.d0,        160.d0,         80.d0]
tuning[*,18]  = [         0.d0,          0.d0,        160.d0]
tuning[*,19]  = [        80.d0,          0.d0,        160.d0]
tuning[*,20]  = [       160.d0,          0.d0,        160.d0]
tuning[*,21]  = [         0.d0,         80.d0,        160.d0]
tuning[*,22]  = [        80.d0,         80.d0,        160.d0]
tuning[*,23]  = [       160.d0,         80.d0,        160.d0]
tuning[*,24]  = [         0.d0,        160.d0,        160.d0]
tuning[*,25]  = [        80.d0,        160.d0,        160.d0]
tuning[*,26]  = [       160.d0,        160.d0,        160.d0]

FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0

; THE FSRs
;------------------------------------------------------------

FSR           = DBLARR(7)    ; FSR in Angstrom
FSR[0]        = 0.172457     ; for the narrow-band Michelson
FSR[1]        = 0.344242     ; for the broad-band  Michelson
FSR[2]        = 0.702        ; for E1
FSR[3]        = 1.405        ; for E2
FSR[4]        = 2.779        ; for E3
FSR[5]        = 5.682        ; for E4
FSR[6]        = 11.354       ; for E5


; RESTORE THE ASSUMED FRONT WINDOW+BLOCKER FILTER PROFILE
;------------------------------------------------------------

RESTORE,'frontwindow.bin'
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')  
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam);,/LSQUADRATIC)

; WE RECONSTRUCT THE INTENSITIES
;------------------------------------------------------------

Inten_rec     = FLTARR(nx,nx,nseq) ; intensities

profilef  = blocker0
FOR i=3,6 DO profilef = profilef*(1.d0+contrast*COS(2.d0*!dpi/FSR[i]*lam+phase))/2.d0

WINDOW,0,RETAIN=2,XSIZE=1000,YSIZE=900

FOR j=0,nx-1 DO BEGIN
PRINT,'j=',j
!P.MULTI=[0,5,6]
FOR k=0,nseq-1 DO TVIM,Inten_rec[*,*,k],/scale

    FOR i=0,nx-1 DO BEGIN
        FOR k=0,nseq-1 DO BEGIN
            profile = profilef
            FOR l=0,2 DO profile = profile * (1.d0+contrast*COS(2.d0*!dpi/FSR[l]*lam+Phig[i,j,l]+tuning[l,k]))/2.d0
            
            line = continuum[i,j] - linedepth[i,j] * EXP(-(lam-lam0)^2.d0/linewidth[i,j]^2.d0)

            Inten_rec[i,j,k] = TOTAL(line * profile)*dlam

        ENDFOR
    ENDFOR


ENDFOR

;SAVE RESULT OF THE RECONSTRUCTION
;---------------------------------------------------------------

SAVE,Inten_rec,FILE='RESULTS/RECONSTRUCTION_CAL_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

; COMPARE RECONSTRUCTION WITH MEASURED INTENSITIES
;---------------------------------------------------------------

draw:

RESTORE,FILE='RESULTS/RECONSTRUCTION_CAL_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;Inten_rec = Inten

RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(256)),1)+'.BIN'
Inten = REBIN(imx,nx,nx,nseq) ;!!!!!!!!!!!!
imx   = 0.0                   ;!!!!!!!!!!!!

center=[65,63];nx/2] ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
distance=SHIFT(DIST(nx,nx),center[0],center[1])*0.5d0*4096.d0/nx
a=WHERE(distance LT 930.d0,COMPLEMENT=b)
distance[a]=1.d0
distance[b]=0.d0
FOR k=0,nseq-1 DO BEGIN
    Inten[*,*,k]=Inten[*,*,k]*distance
    Inten_rec[*,*,k]=Inten_rec[*,*,k]*distance
ENDFOR

b = WHERE(Inten EQ 0.0)
Inten[b] = 0.00000000001d0


SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xsize=23,ysize=28,bits=24,/color,xoffset=-1,yoffset=0
LOADCT,3
!P.MULTI=[0,2,4]
FOR k=0,nseq-1 DO BEGIN
    temp = REFORM((Inten_rec[*,*,k]-Inten[*,*,k])/Inten[*,*,k])
    TVIM,temp,/scale,range=[-.04,.04],barwidth=0.5,pcharsize=1.5,tit='!17'+STRING(k)
ENDFOR
FOR k=0,nseq-1 DO BEGIN
    temp = REFORM((Inten_rec[*,*,k]-Inten[*,*,k])/Inten[*,*,k])
    aa   = WHERE(temp LT 0.06 AND temp GT -0.06,na,COMPLEMENT=b)
    hist = histogram(temp[aa],binsize=0.001,min=-0.04,max=0.04)
    plot,FINDGEN(N_ELEMENTS(hist))*0.001-0.04,hist/FLOAT(na),psym=10,tit='!7r!17='+STRING(SIGMA(temp[aa]))+'/!7l!17='+STRING(MEAN(temp[aa])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5
    PRINT,MIN(temp[aa]),MAX(temp[aa]),SIGMA(temp[aa]),MEAN(temp[aa])
ENDFOR

DEVICE,/close
SET_PLOT,'x'
!P.MULTI=0


READ,PAUSE

END
