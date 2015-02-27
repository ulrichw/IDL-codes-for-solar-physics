; PROGRAM TO PERFORM THE ANGULAR DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE CO-TUNE SEQUENCE
; OBJECTIVE: TO DERIVE THE ANGULAR DEPENDENCE OF THE TRANSMISSION
; PROFILE
; IN CALMODE
; 1) WE FIT FOR THE PHASES OF THE TUNABLE ELEMENTS USING CONTRASTS
;    AND PHASES OF THE NON-TUNABLE PART ALREADY OBTAINED WITH
;    HMI_LASER1.PRO AND HMI_LASER2.PRO
; 2) WE FIT FOR THE SHIFT OF THE NON-TUNABLE PART
;
; NB: YOU NEED TO RUN THE CODE ONCE PER INCIDENCE ANGLE
;----------------------------------------------------------------------

PRO HMI_laser3

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 128              ; number of rows
ny          = 128              ; number of columns

Bg0         = DBLARR(nx,ny,3)  ; contrasts of the tunable elements
Phig0       = DBLARR(nx,ny,3)  ; phases of the tunable elements
I0          = DBLARR(nx,ny)    ; relative laser intensity
nseq        = 27               ; number of positions in the de-tune sequence
Inten       = DBLARR(nx,ny,nseq);measured output intensities (measured on a HMI CCD)
Inten_tunable=DBLARR(nseq)
tuning      = DBLARR(nseq,3)   ; tuning positions for the detune sequence
nlam        = 8000             ; number of wavelengths           
dlam        = 1d0/1.0d3      ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelength values
nl0         = 5      ; number of laser wavelengths
Inten_shifted= DBLARR(nl0)
Inten_nontunable=DBLARR(nl0)
nshift      = 300
lshift      = (DINDGEN(nshift)-(nshift-1)/2.)*dlam ; shifts between -149 and +149 mA with a precision of 1 mA
Residual    = DBLARR(nshift)    ; residual for the least-squares fit

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


; MEASURED INTENSITY MAPS
; AND DYE LASER CENTRAL WAVELENGTH
; result of the calibration test
;----------------------------------------


; list containing the names of filtergrams
list = 'XXXXXXXXXX'
base = 'XXXXXXXXXX'

FOR i=0,nl0-1 DO BEGIN
    header       = STRARR(nseq,91)
    lam1         = DBLARR(nseq) ; laser wavelengths in Angstrom
    Inten[*,*,i*nseq:(i+1)*nseq-1] = READFILE2(header,nx,ny,lam1,base,list)
    IF(TOTAL(lam1[1:nseq-1]-lam1[0:nseq-2]) NE 0.d0) THEN STOP
    lam0[i]      = lam1[0]
ENDFOR


; FIRST PART: OBTAINING THE TUNABLE PART
; RELATIVE PHASES
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

FOR jjj=0,ny-1 DO BEGIN
    FOR iii=0,nx-1 DO BEGIN

  ; GUESS OF THE PHASES AND CONTRASTS OF TUNABLE ELEMENTS 
  ; ASSUMING THE LASER LINE PROFILE IS A DELTA FUNCTION
  ; we use the fact that cos(a)+cos(a+2pi/3)+cos(a+4pi/3) = 0
  ; I0 = [ cos(2\pi l0/FSR0)         -sin(2\pi l0/FSR0)         ] * B cos(\phi) * T/2 + T/2
  ; I1   [ cos(2\pi l0/FSR0+2\pi/3)  -sin(2\pi l0/FSR0+2\pi/3)  ]   B sin(\phi) * T/2 + T/2
  ; with 3 T/2 = I0+I1+I2
  ; USING THE MEASURED INTENSITIES
  ;----------------------------------------------------------


        IF(distance[iii,jjj] LE 960.d0) THEN BEGIN ; XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
            FOR i=0,2 DO BEGIN
                CASE i OF
                    0: BEGIN
                        I00=  Inten[iii,jjj,0]+Inten[iii,jjj,1]+Inten[iii,jjj,2]+Inten[iii,jjj,3]+Inten[iii,jjj,4]+Inten[iii,jjj,9]+Inten[iii,jjj,10]+Inten[iii,jjj,11]+Inten[iii,jjj,12]
                        I1 =  Inten[iii,jjj,5]+Inten[iii,jjj,7]+Inten[iii,jjj,13]+Inten[iii,jjj,14]+Inten[iii,jjj,15]+Inten[iii,jjj,16]+Inten[iii,jjj,17]+Inten[iii,jjj,18]+Inten[iii,jjj,19]
                        I2 =  Inten[iii,jjj,6]+Inten[iii,jjj,8]+Inten[iii,jjj,20]+Inten[iii,jjj,21]+Inten[iii,jjj,22]+Inten[iii,jjj,23]+Inten[iii,jjj,24]+Inten[iii,jjj,25]+Inten[iii,jjj,26]
                    END
                    1: BEGIN
                        I00=  Inten[iii,jjj,0]+Inten[iii,jjj,1]+Inten[iii,jjj,2]+Inten[iii,jjj,5]+Inten[iii,jjj,6]+Inten[iii,jjj,13]+Inten[iii,jjj,14]+Inten[iii,jjj,20]+Inten[iii,jjj,21]
                        I1 =  Inten[iii,jjj,3]+Inten[iii,jjj,7]+Inten[iii,jjj,9]+Inten[iii,jjj,10]+Inten[iii,jjj,15]+Inten[iii,jjj,17]+Inten[iii,jjj,22]+Inten[iii,jjj,24]+Inten[iii,jjj,25]
                        I2 =  Inten[iii,jjj,4]+Inten[iii,jjj,8]+Inten[iii,jjj,11]+Inten[iii,jjj,12]+Inten[iii,jjj,16]+Inten[iii,jjj,18]+Inten[iii,jjj,19]+Inten[iii,jjj,23]+Inten[iii,jjj,26]               
                    END
                    2: BEGIN
                        I00=  Inten[iii,jjj,0]+Inten[iii,jjj,3]+Inten[iii,jjj,4]+Inten[iii,jjj,5]+Inten[iii,jjj,6]+Inten[iii,jjj,15]+Inten[iii,jjj,16]+Inten[iii,jjj,22]+Inten[iii,jjj,23]
                        I1 =  Inten[iii,jjj,1]+Inten[iii,jjj,7]+Inten[iii,jjj,9]+Inten[iii,jjj,11]+Inten[iii,jjj,13]+Inten[iii,jjj,18]+Inten[iii,jjj,20]+Inten[iii,jjj,24]+Inten[iii,jjj,26]
                        I2 =  Inten[iii,jjj,2]+Inten[iii,jjj,8]+Inten[iii,jjj,10]+Inten[iii,jjj,12]+Inten[iii,jjj,14]+Inten[iii,jjj,17]+Inten[iii,jjj,19]+Inten[iii,jjj,21]+Inten[iii,jjj,25]            
                    END
                ENDCASE
                tot      = I00+I1+I2
                vec[0,0] = I00*3.d0/tot-1.d0
                vec[0,1] = I1 *3.d0/tot-1.d0
                vec      = REFORM(mat[i,*,*])##vec
                Bg0[iii,jjj,i]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)
                IF(Bg0[iii,jjj,i] GT 1.d0) THEN Bg0[iii,jjj,i] = 1.d0
                Phig0[iii,jjj,i] = ATAN(vec[0,1],vec[0,0])
            ENDFOR
            
        ENDIF
        
        PRINT,'PHASES'   ,phig0[iii,jjj,*]
        PRINT,'CONTRASTS',Bg0  [iii,jjj,*]
        PRINT,iii,jjj
        
    ENDFOR
ENDFOR

;----------------------------------------------------------------------
;
; SECOND PART: OBTAINING THE NON-TUNABLE PART
; WAVELENGTH DRIFT
;
;----------------------------------------------------------------------

RESTORE,'SUNTEST_LASER2.BIN' ; Bg And Phig

; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;-----------------------------------------------

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

;blocker0     = DBLARR(nlam)+1.d0 
q            = READFITS('blocker11.fits')           ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam);,/LSQUADRATIC)

; non-tunable part transmission profile for a normal incidence ray
profile      = blocker0
FOR i=3,6 DO profile = profile * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam+Phig[iii,jjj,i]) )/2.d0 

a0 = WHERE(ABS(lam-lam0[0]) EQ MIN(ABS(lam-lam0[0])) )
a1 = WHERE(ABS(lam-lam0[1]) EQ MIN(ABS(lam-lam0[1])) )
a2 = WHERE(ABS(lam-lam0[2]) EQ MIN(ABS(lam-lam0[2])) )
a3 = WHERE(ABS(lam-lam[lam0[3]]) EQ MIN(ABS(lam-lam[lam0[3]])) )
a4 = WHERE(ABS(lam-lam[lam0[4]]) EQ MIN(ABS(lam-lam[lam0[4]])) )
a  = [a0[0],a1[0],a2[0],a3[0],a4[0]]

FOR jjj=0,ny-1 DO BEGIN
    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] LE 960.d0) THEN BEGIN ; XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          ; Estimate of I0
            FOR k=0,nl-1 DO BEGIN
                aaa       = WHERE(ABS(lam-lam[a[k]]) EQ MIN(ABS(lam-lam[a[k]]))) 
                profileg  = blocker0[aaa[0]] ; interpolation instead ?
                FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam[a[k]]+Phig[iii,jjj,i]) )/2.d0
                I0[iii,jjj] = TOTAL(Inten[iii,jjj,*]/27.d0*8.d0/profileg)+I0[iii,jjj]
            ENDFOR
            I0[iii,jjj]=I0[iii,jjj]/DOUBLE(nl) ; averaging over the entire detune sequence

           
            FOR i=0,nl0-1 DO BEGIN
                FOR j=0,nseq-1 DO BEGIN    
                    profileg = DBLARR(nlam)+1.d0
                    FOR ii=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,ii]*COS(dpi/FSR[ii]*lam[a[i]]+Phig[iii,jjj,ii]+tuning[ii,j]))/2.d0
                    Inten_tunable[j] = profileg
                ENDFOR
                temp                = Inten[iii,jjj,*] / Inten_tunable[*]
                Inten_nontunable[i] = MEAN(temp) ; MEDIAN(temp)
            ENDFOR
            
            FOR j=0,nshift-1 DO BEGIN
                FOR k=0,nl-1 DO BEGIN
                    aaa = WHERE(ABS(lam-lam[a[k]]-lshift[j]) EQ MIN(ABS(lam-lam[a[k]]-lshift[j]))) 
                    profile_shifted         = blocker0[aaa[0]]
                    FOR i=3,6 DO profile_shifted = profile_shifted * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*(lam[aaa[0]]-lshift[j])+Phig[iii,jjj,i]) )/2.d0

                    Inten_shifted[k] = profile_shifted * I0[iii,jjj]
                ENDFOR
                Residual[j] = TOTAL( (Inten_shifted-Inten_nontunable)^2.d0 )
            ENDFOR
            
            temp = WHERE(Residual EQ MIN(Residual))
            shift[iii,jjj] = lshift[temp[0]]

        ENDIF

    ENDFOR
ENDFOR

SAVE,Bg0,Phig0,shift,FILE='SUNTEST_LASER3.BIN'
PRINT,'TOTAL TIME ELAPSED',SYSTIME(1)-time0

READ,pause

END
