; PROGRAM TO PERFORM PART OF THE ANGULAR DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW THE PHASES AND CONTRASTS OF THE TUNABLE ELEMENTS
; AND TO SEE HOW THAT VARIES AT DIFFERENT WAVELENGTHS AND ANGLES
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_jspcudetune180,draw

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0 ; solar Fe I 6173 line central wavelength in air
; when phases are at 0, we are centered on this wavelength

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 128              ; number of rows
ny          = 128              ; number of columns

Bg0         = FLTARR(nx,ny,3,6)  ; contrasts of the tunable elements
Phig0       = FLTARR(nx,ny,3,6)  ; phases of the tunable elements
nseq        = 27               ; number of positions in the de-tune sequence
tuning      = FLTARR(3,nseq)   ; tuning positions for the detune sequence

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;--------------------------------------------

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172457         ; for the narrow-band Michelson
FSR[1]      = 0.344242         ; for the broad-band  Michelson
FSR[2]      = 0.702            ; for E1
FSR[3]      = 1.405            ; for E2
FSR[4]      = 2.779            ; for E3
FSR[5]      = 5.682            ; for E4
FSR[6]      = 11.354           ; for E5


; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
;                NB MICHELSON   WB MICHELSON         E1 LYOT
;-----------------------------------------------------------


tuning[*,0]  = [         0.d0,          0.d0,          0.d0]
tuning[*,1]  = [        80.d0,          0.d0,          0.d0]
tuning[*,2]  = [       160.d0,          0.d0,          0.d0]
tuning[*,3]  = [         0.d0,         80.d0,          0.d0]
tuning[*,4]  = [        80.d0,         80.d0,          0.d0]
tuning[*,5]  = [       160.d0,         80.d0,          0.d0]
tuning[*,6]  = [         0.d0,        160.d0,          0.d0]
tuning[*,7]  = [        80.d0,        160.d0,          0.d0]
tuning[*,8]  = [       160.d0,        160.d0,          0.d0]
tuning[*,9]  = [         0.d0,          0.d0,         80.d0]
tuning[*,10] = [        80.d0,          0.d0,         80.d0]
tuning[*,11] = [       160.d0,          0.d0,         80.d0]
tuning[*,12] = [         0.d0,         80.d0,         80.d0] 
tuning[*,13] = [        80.d0,         80.d0,         80.d0]
tuning[*,14] = [       160.d0,         80.d0,         80.d0]
tuning[*,15] = [         0.d0,        160.d0,         80.d0]
tuning[*,16] = [        80.d0,        160.d0,         80.d0]
tuning[*,17] = [       160.d0,        160.d0,         80.d0]
tuning[*,18] = [         0.d0,          0.d0,        160.d0]
tuning[*,19] = [        80.d0,          0.d0,        160.d0]
tuning[*,20] = [       160.d0,          0.d0,        160.d0]
tuning[*,21] = [         0.d0,         80.d0,        160.d0]
tuning[*,22] = [        80.d0,         80.d0,        160.d0]
tuning[*,23] = [       160.d0,         80.d0,        160.d0]
tuning[*,24] = [         0.d0,        160.d0,        160.d0]
tuning[*,25] = [        80.d0,        160.d0,        160.d0]
tuning[*,26] = [       160.d0,        160.d0,        160.d0]


FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0

; INFORMATION TO BE ADDED BY THE USER
;----------------------------------------

;list='list_detune180_060303_181502'
list='list_detune180_060303_042821'
;list='list_detune180_060303_003721'
;list='list_detune180_060303_023201'
;list='list_detune180_060303_230038'

IF(draw EQ 1) THEN GOTO,draw


GOTO,dejalu
readimages,list,images,headers,time,nbin=nx,iover=iover
SAVE,images,time,FILE='SEQUENCE0_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;RESTORE,'SEQUENCE0_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

; TO SUPPRESS THE SPIKES IN THE DARK FRAME DRIFT
corner             = REBIN(images(0:nx/64-1,0:nx/64-1,*),1,1,180)
FOR i=0,179 DO images[*,*,i] = images[*,*,i]-corner[i]

Inten              = FLTARR(nx,ny,162);measured output intensities (measured on a HMI CCD)
Darks              = FLTARR(nx,ny,12)
Calmode            = FLTARR(nx,nx,6)
timeimages         = FLTARR(162)

Darks[*,*,0]       = images[*,*,1]
Darks[*,*,1]       = images[*,*,29]
Darks[*,*,2]       = images[*,*,31]
Darks[*,*,3]       = images[*,*,59]
Darks[*,*,4]       = images[*,*,61]
Darks[*,*,5]       = images[*,*,89]
Darks[*,*,6]       = images[*,*,91]
Darks[*,*,7]       = images[*,*,119]
Darks[*,*,8]       = images[*,*,121]
Darks[*,*,9]       = images[*,*,149]
Darks[*,*,10]      = images[*,*,151]
Darks[*,*,11]      = images[*,*,179]

IF (list EQ 'list_detune180_060303_023201') THEN BEGIN  ; bad darks
Darks[*,*,5]       = Darks[*,*,4]
ENDIF
IF (list EQ 'list_detune180_060303_230038') THEN BEGIN  ; bad darks
Darks[*,*,2]       = Darks[*,*,3]
ENDIF

timedarks          = [time[1],time[29],time[31],time[59],time[61],time[89],time[91],time[119],time[121],time[149],time[151],time[179]]

Calmode[*,*,0]     = images[*,*,0]
Calmode[*,*,1]     = images[*,*,30]
Calmode[*,*,2]     = images[*,*,60]
Calmode[*,*,3]     = images[*,*,90]
Calmode[*,*,4]     = images[*,*,120]
Calmode[*,*,5]     = images[*,*,150]

Inten[*,*,0:26]    = images[*,*,2:28]
Inten[*,*,27:53]   = images[*,*,32:58]
Inten[*,*,54:80]   = images[*,*,62:88]
Inten[*,*,81:107]  = images[*,*,92:118]
Inten[*,*,108:134] = images[*,*,122:148]
Inten[*,*,135:161] = images[*,*,152:178]
timeimages[0:26]   = time[2:28]
timeimages[27:53]  = time[32:58]
timeimages[54:80]  = time[62:88]
timeimages[81:107] = time[92:118]
timeimages[108:134]= time[122:148]
timeimages[135:161]= time[152:178]

SET_PLOT,'x'
WINDOW,0,RETAIN=2
PLOT,timedarks,TOTAL(TOTAL(darks[*,*,*],1),1),xst=1,yst=1

; WE REMOVE THE DARK FRAME
PRINT,'DARK FRAME REMOVAL'
    FOR j=0,ny-1 DO BEGIN
        FOR jj=0,nx-1 DO BEGIN
            Inten[jj,j,*] = Inten[jj,j,*] - INTERPOL(REFORM(Darks[jj,j,*]),timedarks,timeimages[*])
        ENDFOR
    ENDFOR
a    = WHERE(Inten LT 0.d0)   
IF (a[0] NE -1) THEN Inten[a]= 0.d0


SAVE,Inten,FILE='SEQUENCE_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
images = 0.0
time   = 0.0
dejalu:
RESTORE,'SEQUENCE_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

; TO GET RID OF DUBIOUS IMAGES
;----------------------------------------

IF (list EQ 'list_detune180_060303_181502') THEN BEGIN
    a=[12, 49, 56, 68, 101, 118, 119, 125] ; bad exposures and overscans
    Inten[*,*,a]=-1.0
ENDIF

IF (list EQ 'list_detune180_060303_042821') THEN BEGIN
    a=[42, 151, 13, 130]
    Inten[*,*,a]=-1.0
ENDIF

IF (list EQ 'list_detune180_060303_003721') THEN BEGIN
    a=[95,116,133,135,139]
    Inten[*,*,a]=-1.0
ENDIF

IF (list EQ 'list_detune180_060303_023201') THEN BEGIN
    a=[14,118,146]
    Inten[*,*,a]=-1.0
ENDIF

IF (list EQ 'list_detune180_060303_230038') THEN BEGIN
    a=[62,117-11,139-14,140-14]
    Inten[*,*,a]=-1.0
ENDIF

IF (list EQ 'list_detune180_060303_181502') THEN lam0          = 6173.467
IF (list EQ 'list_detune180_060303_042821') THEN lam0          = 6173.394
IF (list EQ 'list_detune180_060303_003721') THEN lam0          = 6173.628;15
IF (list EQ 'list_detune180_060303_023201') THEN lam0          = 6173.210
IF (list EQ 'list_detune180_060303_230038') THEN lam0          = 6173.077


; BEGINNING OF THE ITERATIONS
;----------------------------------------

mat      = DBLARR(3,2,2)
vec      = DBLARR(1,2)
FOR i=0,2 DO BEGIN
    mat[i,0,0] =  COS(dpi/FSR[i]*lam0)
    mat[i,1,0] = -SIN(dpi/FSR[i]*lam0)
    mat[i,0,1] =  COS(dpi/FSR[i]*lam0 + dpi/3.d0)
    mat[i,1,1] = -SIN(dpi/FSR[i]*lam0 + dpi/3.d0)
    mat[i,*,*] =  LA_INVERT(REFORM(mat[i,*,*]),/DOUBLE)
ENDFOR

distance = SHIFT(DIST(nx,ny),nx/2,ny/2)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center

In = FLTARR(nseq)

FOR hole=0,5 DO BEGIN

FOR jjj=0,ny-1 DO BEGIN
    FOR iii=0,nx-1 DO BEGIN
  ; GUESS OF THE PHASES AND CONTRASTS OF TUNABLE ELEMENTS 
  ; ASSUMING THE LASER LINE PROFILE IS A DELTA FUNCTION
  ; we use the fact that cos(a)+cos(a+2pi/3)+cos(a+4pi/3) = 0
  ; I0 = [ cos(2\pi l0/FSR0)         -sin(2\pi l0/FSR0)         ] * B cos(\phi) * T/2 + T/2
  ; I1   [ cos(2\pi l0/FSR0+2\pi/3)  -sin(2\pi l0/FSR0+2\pi/3)  ]   B sin(\phi) * T/2 + T/2
  ; with 3 T/2 = I0+I1+I2
  ; USING THE MEASURED INTENSITIES
  ; THE FOLLOWING PROCEDURE IS WRITTEN FOR [1.5,1.5,1.5]
  ; IF ONE OF THE COEFFICIENT IS -1.5, THEN I1=I2 AND I2=I1 
  ;----------------------------------------------------------

        IF(distance[iii,jjj] LE 2024.d0) THEN BEGIN
           FOR i=0,2 DO BEGIN
           CASE i OF
               0: BEGIN ; CHANGED FOR [1.5,1.5,-1.5]
                   In[*] = Inten[iii,jjj,hole*nseq:(hole+1)*nseq-1]
                 ; TO TAKE CARE OF THE DUBIOUS IMAGES
                   IF(In[0] EQ -1.0 OR In[1] EQ -1.0 OR In[2] EQ -1.0) THEN BEGIN
                       In[0] = 0.0
                       In[1] = 0.0
                       In[2] = 0.0
                   ENDIF
                   IF(In[3] EQ -1.0 OR In[4] EQ -1.0 OR In[5] EQ -1.0) THEN BEGIN
                       In[3] = 0.0
                       In[4] = 0.0
                       In[5] = 0.0
                   ENDIF
                   IF(In[6] EQ -1.0 OR In[7] EQ -1.0 OR In[8] EQ -1.0) THEN BEGIN
                       In[6] = 0.0
                       In[7] = 0.0
                       In[8] = 0.0
                   ENDIF
                   IF(In[9] EQ -1.0 OR In[10] EQ -1.0 OR In[11] EQ -1.0) THEN BEGIN
                       In[9] = 0.0
                       In[10] = 0.0
                       In[11] = 0.0
                   ENDIF
                   IF(In[12] EQ -1.0 OR In[13] EQ -1.0 OR In[14] EQ -1.0) THEN BEGIN
                       In[12] = 0.0
                       In[13] = 0.0
                       In[14] = 0.0
                   ENDIF
                   IF(In[15] EQ -1.0 OR In[16] EQ -1.0 OR In[17] EQ -1.0) THEN BEGIN
                       In[15] = 0.0
                       In[16] = 0.0
                       In[17] = 0.0
                   ENDIF
                   IF(In[18] EQ -1.0 OR In[19] EQ -1.0 OR In[20] EQ -1.0) THEN BEGIN
                       In[18] = 0.0
                       In[19] = 0.0
                       In[20] = 0.0
                   ENDIF
                   IF(In[21] EQ -1.0 OR In[22] EQ -1.0 OR In[23] EQ -1.0) THEN BEGIN
                       In[21] = 0.0
                       In[22] = 0.0
                       In[23] = 0.0
                   ENDIF
                   IF(In[24] EQ -1.0 OR In[25] EQ -1.0 OR In[26] EQ -1.0) THEN BEGIN
                       In[24] = 0.0
                       In[25] = 0.0
                       In[26] = 0.0
                   ENDIF
                   I00= (In[0]+In[3]+In[6]+In[9]+In[12]+In[15]+In[18]+In[21]+In[24])
                   I1 = (In[1]+In[4]+In[7]+In[10]+In[13]+In[16]+In[19]+In[22]+In[25])
                   I2 = (In[2]+In[5]+In[8]+In[11]+In[14]+In[17]+In[20]+In[23]+In[26])
               END
               1: BEGIN
                   In[*] = Inten[iii,jjj,hole*nseq:(hole+1)*nseq-1]
                   IF(In[0] EQ -1.0 OR In[3] EQ -1.0 OR In[6] EQ -1.0) THEN BEGIN
                       In[0] = 0.0
                       In[3] = 0.0
                       In[6] = 0.0
                   ENDIF
                   IF(In[1] EQ -1.0 OR In[4] EQ -1.0 OR In[7] EQ -1.0) THEN BEGIN
                       In[1] = 0.0
                       In[4] = 0.0
                       In[7] = 0.0
                   ENDIF
                   IF(In[2] EQ -1.0 OR In[5] EQ -1.0 OR In[8] EQ -1.0) THEN BEGIN
                       In[2] = 0.0
                       In[5] = 0.0
                       In[8] = 0.0
                   ENDIF
                   IF(In[9] EQ -1.0 OR In[12] EQ -1.0 OR In[15] EQ -1.0) THEN BEGIN
                       In[9] = 0.0
                       In[12] = 0.0
                       In[15] = 0.0
                   ENDIF
                   IF(In[10] EQ -1.0 OR In[13] EQ -1.0 OR In[16] EQ -1.0) THEN BEGIN
                       In[10] = 0.0
                       In[13] = 0.0
                       In[16] = 0.0
                   ENDIF
                   IF(In[11] EQ -1.0 OR In[14] EQ -1.0 OR In[17] EQ -1.0) THEN BEGIN
                       In[11] = 0.0
                       In[14] = 0.0
                       In[17] = 0.0
                   ENDIF
                   IF(In[18] EQ -1.0 OR In[21] EQ -1.0 OR In[24] EQ -1.0) THEN BEGIN
                       In[18] = 0.0
                       In[21] = 0.0
                       In[24] = 0.0
                   ENDIF
                   IF(In[19] EQ -1.0 OR In[22] EQ -1.0 OR In[25] EQ -1.0) THEN BEGIN
                       In[19] = 0.0
                       In[22] = 0.0
                       In[25] = 0.0
                   ENDIF
                   IF(In[20] EQ -1.0 OR In[23] EQ -1.0 OR In[26] EQ -1.0) THEN BEGIN
                       In[20] = 0.0
                       In[23] = 0.0
                       In[26] = 0.0
                   ENDIF
                   I00= (In[0]+In[1]+In[2]+In[9]+In[10]+In[11]+In[18]+In[19]+In[20])
                   I1 = (In[3]+In[4]+In[5]+In[12]+In[13]+In[14]+In[21]+In[22]+In[23])
                   I2 = (In[6]+In[7]+In[8]+In[15]+In[16]+In[17]+In[24]+In[25]+In[26])
               END
               2: BEGIN
                   In[*] = Inten[iii,jjj,hole*nseq:(hole+1)*nseq-1]
                   IF(In[0] EQ -1.0 OR In[9] EQ -1.0 OR In[18] EQ -1.0) THEN BEGIN
                       In[0] = 0.0
                       In[9] = 0.0
                       In[18] = 0.0
                   ENDIF
                   IF(In[1] EQ -1.0 OR In[10] EQ -1.0 OR In[19] EQ -1.0) THEN BEGIN
                       In[1] = 0.0
                       In[10] = 0.0
                       In[19] = 0.0
                   ENDIF
                   IF(In[2] EQ -1.0 OR In[11] EQ -1.0 OR In[20] EQ -1.0) THEN BEGIN
                       In[2] = 0.0
                       In[11] = 0.0
                       In[20] = 0.0
                   ENDIF
                   IF(In[3] EQ -1.0 OR In[12] EQ -1.0 OR In[21] EQ -1.0) THEN BEGIN
                       In[3] = 0.0
                       In[12] = 0.0
                       In[21] = 0.0
                   ENDIF
                   IF(In[4] EQ -1.0 OR In[13] EQ -1.0 OR In[22] EQ -1.0) THEN BEGIN
                       In[4] = 0.0
                       In[13] = 0.0
                       In[22] = 0.0
                   ENDIF
                   IF(In[5] EQ -1.0 OR In[14] EQ -1.0 OR In[23] EQ -1.0) THEN BEGIN
                       In[5] = 0.0
                       In[14] = 0.0
                       In[23] = 0.0
                   ENDIF
                   IF(In[6] EQ -1.0 OR In[15] EQ -1.0 OR In[24] EQ -1.0) THEN BEGIN
                       In[6] = 0.0
                       In[15] = 0.0
                       In[24] = 0.0
                   ENDIF
                   IF(In[7] EQ -1.0 OR In[16] EQ -1.0 OR In[25] EQ -1.0) THEN BEGIN
                       In[7] = 0.0
                       In[16] = 0.0
                       In[25] = 0.0
                   ENDIF
                   IF(In[8] EQ -1.0 OR In[17] EQ -1.0 OR In[26] EQ -1.0) THEN BEGIN
                       In[8] = 0.0
                       In[17] = 0.0
                       In[26] = 0.0
                   ENDIF
                   I00= (In[0]+In[1]+In[2]+In[3]+In[4]+In[5]+In[6]+In[7]+In[8])
                   I2 = (In[9]+In[10]+In[11]+In[12]+In[13]+In[14]+In[15]+In[16]+In[17])
                   I1 = (In[18]+In[19]+In[20]+In[21]+In[22]+In[23]+In[24]+In[25]+In[26])
               END
           ENDCASE

                tot      = I00+I1+I2
                vec[0,0] = I00*3.d0/tot-1.d0
                vec[0,1] = I1 *3.d0/tot-1.d0
                vec      = REFORM(mat[i,*,*])##vec
                Bg0[iii,jjj,i,hole]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)

               ;WE CLEAN THE SOLUTION
                IF(Bg0[iii,jjj,i,hole] GT 1.6d0) THEN Bg0[iii,jjj,i,hole] = 1.6d0
                Phig0[iii,jjj,i,hole] = ATAN(vec[0,1],vec[0,0])
            ENDFOR
            
        ENDIF ELSE BEGIN

            Bg0[iii,jjj,*,hole]=0.d0
            Phig0[iii,jjj,*,hole]=0.d0

        ENDELSE
        
      ; TO RECONSTRUCT THE INTENSITIES AT THE WAVELENGTH lam2
       ; FOR j=0,nseq-1 DO BEGIN
       ;     profile = blocker0
       ;     FOR i=0,2 DO profile = profile * (1.d0+Bg0[iii,jjj,i]*COS(2.d0*!dpi/FSR[i]*lam2+Phig0[iii,jjj,i]+tuning[i,j]))/2.d0
       ;    FOR i=3,6 DO profile = profile * (1.d0+0.99d0*COS(2.d0*!dpi/FSR[i]*lam2))/2.d0
       ;     Intenrecons[iii,jjj,j] = profile
       ; ENDFOR
                
    ENDFOR
ENDFOR

ENDFOR

SAVE,Bg0,Phig0,nx,FILE='RESULTS/RESULTS_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'TOTAL TIME ELAPSED',SYSTIME(1)-time0

draw:
RESTORE,'RESULTS/RESULTS_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'


; WE APPLY A MASK TO KEEP ONLY THE PHASES AND CONTRASTS ON THE ELEMENT
; APERTURE
xcenter=nx/2
ycenter=nx/2
distance = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx
moy=DBLARR(3,6)
mini=DBLARR(3,6)
maxi=DBLARR(3,6)
mini2=DBLARR(3,6)
maxi2=DBLARR(3,6)

FOR hole=0,5 DO BEGIN
a=WHERE(distance GE 1024,COMPLEMENT=b)

FOR i=0,2 DO BEGIN
    temp=REFORM(Phig0[*,*,i,hole])
    temp[a]=-10000.d0
    phig0[*,*,i,hole]=temp[*,*]
    temp=REFORM(Bg0[*,*,i,hole])
    temp[a]=-10000.d0
    Bg0[*,*,i,hole]=temp[*,*]
   ;Phig0[0:10,*,i]=-10000.d0
   ;Phig0[*,490:511,i]=-10000.d0
   ;Bg0[0:10,*,i]=-10000.d0
   ;Bg0[*,490:511,i]=-10000.d0
ENDFOR

a = WHERE(FINITE(Phig0) EQ 0)
IF(a[0] NE -1) THEN Phig0[a]=-10000.d0
a = WHERE(FINITE(Bg0) EQ 0)
IF(a[0] NE -1) THEN Bg0[a]=-10000.d0

; WE COMPUTE THE AVERAGE VALUE

FOR i=0,2 DO BEGIN
    temp=REFORM(Phig0[*,*,i,hole])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    moy[i,hole]=MEAN(temp[a])
    mini[i,hole]=MIN(temp[a])*180./!dpi
    maxi[i,hole]=MAX(temp[a])*180./!dpi
    temp[b]=moy[i,hole]
    Phig0[*,*,i,hole]=temp
    temp=REFORM(Bg0[*,*,i,hole])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    temp[b]=MEAN(temp[a])
    mini2[i,hole]=MIN(temp[a])
    maxi2[i,hole]=MAX(temp[a])
    Bg0[*,*,i,hole]=temp
ENDFOR

ENDFOR

; WE PLOT THE RESULT
SET_PLOT,'ps'
!p.multi=[0,1,2]
device,file='yo.ps',bits=24,xoffset=0,yoffset=1,xsize=15,ysize=22,/color
loadct,3
FOR hole=0,5 DO BEGIN
tvim,phig0[*,*,0,hole]*180.d0/!dpi,/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[0,hole],maxi[0,hole]]
tvim,phig0[*,*,1,hole]*180.d0/!dpi,/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[1,hole],maxi[1,hole]]
tvim,phig0[*,*,2,hole]*180.d0/!dpi,/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[2,hole],maxi[2,hole]]
tvim,Bg0[*,*,0,hole],/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[0,hole],maxi2[0,hole]]
tvim,Bg0[*,*,1,hole],/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[1,hole],maxi2[1,hole]]
tvim,Bg0[*,*,2,hole],/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[2,hole],maxi2[2,hole]]
ENDFOR

DEVICE,/close
PRINT,'AVERAGES'

FOR hole=0,5 DO FOR i=0,2 DO PRINT,moy[i,hole]*180.d0/!dpi



SET_PLOT,'x'
!p.multi=0
READ,pause

END
