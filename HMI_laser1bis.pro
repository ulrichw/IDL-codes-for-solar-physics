; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE DE-TUNE SEQUENCE
; OBJECTIVE: TO LEARN HOW TO CO-TUNE THE ELEMENTS
; AND THE PHASES AND CONTRASTS OF THE TUNABLE ELEMENTS
; WITH 2 LASER WAVELENGTHS TO SUPPRESS THE SIGN AMBIGUITY
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser1bis

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0 ; solar Fe I 6173 line central wavelength in air

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 512              ; number of rows
ny          = 512              ; number of columns

Bg0         = DBLARR(nx,ny,3)  ; contrasts of the tunable elements
Phig0       = DBLARR(nx,ny,3)  ; phases of the tunable elements
nseq        = 27               ; number of positions in the de-tune sequence
Inten1      = DBLARR(nx,ny,nseq);measured output intensities (measured on a HMI CCD)
Inten2      = DBLARR(nx,ny,nseq)
tuning      = DBLARR(3,nseq)   ; tuning positions for the detune sequence
Intenrecons1= FLTARR(nx,ny,nseq)
Intenrecons2= FLTARR(nx,ny,nseq)

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

base         = '/tmp20/schou/hmi060130/'
list         = 'listLaser1'
;list         = 'listLaser2'


;header       = STRARR(nseq,91)
;lam          = DBLARR(nseq) ; laser wavelengths in Angstrom
;Inten[*,*,*] = READFILE2(header,nx,ny,lam,base,list)
;IF(TOTAL(lam[1:nseq-1]-lam[0:nseq-2]) NE 0.d0) THEN STOP
;lam0         = lam[0]
;RESTORE,'TMPREADFILE2.BIN' ; 1st detune
;RESTORE,'TMPREADFILE2_BIS.BIN' ; 2nd detune
;RESTORE,'TMPREADFILE2_TRI.BIN' ; 1st detune with dark frame computed line by line
RESTORE,'TMPREADFILE2_FIV_512.BIN'     ; 1st detune with dark frame computed line by line and for each position
Inten1       = imx
RESTORE,'TMPREADFILE2_BIS_FIV_512.BIN' ; 2nd detune with dark frame computed line by line and for each positio
Inten2       = imx

Inten1[*,*,5] = 0.d0 ; because Inten[*,*,5] had bad exposure
Inten1[*,*,21]= 0.d0
lam1          = 6173.439d0-lamref ; ACCORDING TO JESPER !!!!!!!!!!!!
lam2          = 6173.538d0-lamref

; BEGINNING OF THE ITERATIONS
;----------------------------------------

vec       = DBLARR(1,2)
mat1      = DBLARR(3,2,2)
FOR i=0,2 DO BEGIN
    mat1[i,0,0] =  COS(dpi/FSR[i]*lam1)
    mat1[i,1,0] = -SIN(dpi/FSR[i]*lam1)
    mat1[i,0,1] =  COS(dpi/FSR[i]*lam1 + dpi/3.d0)
    mat1[i,1,1] = -SIN(dpi/FSR[i]*lam1 + dpi/3.d0)
   ;mat1[i,*,*] =  LA_INVERT(REFORM(mat1[i,*,*]),/DOUBLE)
ENDFOR
mat2      = DBLARR(3,2,2)
FOR i=0,2 DO BEGIN
    mat2[i,0,0] =  COS(dpi/FSR[i]*lam2)
    mat2[i,1,0] = -SIN(dpi/FSR[i]*lam2)
    mat2[i,0,1] =  COS(dpi/FSR[i]*lam2 + dpi/3.d0)
    mat2[i,1,1] = -SIN(dpi/FSR[i]*lam2 + dpi/3.d0)
   ;mat2[i,*,*] =  LA_INVERT(REFORM(mat2[i,*,*]),/DOUBLE)
ENDFOR

distance = SHIFT(DIST(nx,ny),nx/2,ny/2)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center

q            = READFITS('blocker11.fits')           ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker01     = INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam1);,/LSQUADRATIC)
blocker02     = INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam2);,/LSQUADRATIC)


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


        IF(distance[iii,jjj] LE 960.d0) THEN BEGIN ; XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
            FOR i=0,2 DO BEGIN
;                CASE i OF
;                    0: BEGIN ; CHANGED FOR [1.5,1.5,-1.5]
;                        I0a= (Inten1[iii,jjj,0]+Inten1[iii,jjj,3]+Inten1[iii,jjj,6]+Inten1[iii,jjj,9]+Inten1[iii,jjj,12]+Inten1[iii,jjj,15]+Inten1[iii,jjj,18]+Inten1[iii,jjj,21]+Inten1[iii,jjj,24])
;                        I1a = (Inten1[iii,jjj,1]+Inten1[iii,jjj,4]+Inten1[iii,jjj,7]+Inten1[iii,jjj,10]+Inten1[iii,jjj,13]+Inten1[iii,jjj,16]+Inten1[iii,jjj,19]+Inten1[iii,jjj,22]+Inten1[iii,jjj,25])
;                        I2a = (Inten1[iii,jjj,2]+Inten1[iii,jjj,5]+Inten1[iii,jjj,8]+Inten1[iii,jjj,11]+Inten1[iii,jjj,14]+Inten1[iii,jjj,17]+Inten1[iii,jjj,20]+Inten1[iii,jjj,23]+Inten1[iii,jjj,26])
;                        I0b= (Inten2[iii,jjj,0]+Inten2[iii,jjj,3]+Inten2[iii,jjj,6]+Inten2[iii,jjj,9]+Inten2[iii,jjj,12]+Inten2[iii,jjj,15]+Inten2[iii,jjj,18]+Inten2[iii,jjj,21]+Inten2[iii,jjj,24])
;                        I1b = (Inten2[iii,jjj,1]+Inten2[iii,jjj,4]+Inten2[iii,jjj,7]+Inten2[iii,jjj,10]+Inten2[iii,jjj,13]+Inten2[iii,jjj,16]+Inten2[iii,jjj,19]+Inten2[iii,jjj,22]+Inten2[iii,jjj,25])
;                        I2b = (Inten2[iii,jjj,2]+Inten2[iii,jjj,5]+Inten2[iii,jjj,8]+Inten2[iii,jjj,11]+Inten2[iii,jjj,14]+Inten2[iii,jjj,17]+Inten2[iii,jjj,20]+Inten2[iii,jjj,23]+Inten2[iii,jjj,26])
;                    END
;                    1: BEGIN
;                        I0a= (Inten1[iii,jjj,0]+Inten1[iii,jjj,1]+Inten1[iii,jjj,2]+Inten1[iii,jjj,9]+Inten1[iii,jjj,10]+Inten1[iii,jjj,11]+Inten1[iii,jjj,18]+Inten1[iii,jjj,19]+Inten1[iii,jjj,20])
;                        I1a= (Inten1[iii,jjj,3]+Inten1[iii,jjj,4]+Inten1[iii,jjj,5]+Inten1[iii,jjj,12]+Inten1[iii,jjj,13]+Inten1[iii,jjj,14]+Inten1[iii,jjj,21]+Inten1[iii,jjj,22]+Inten1[iii,jjj,23])
;                        I2a= (Inten1[iii,jjj,6]+Inten1[iii,jjj,7]+Inten1[iii,jjj,8]+Inten1[iii,jjj,15]+Inten1[iii,jjj,16]+Inten1[iii,jjj,17]+Inten1[iii,jjj,24]+Inten1[iii,jjj,25]+Inten1[iii,jjj,26])
;                        I0b= (Inten2[iii,jjj,0]+Inten2[iii,jjj,1]+Inten2[iii,jjj,2]+Inten2[iii,jjj,9]+Inten2[iii,jjj,10]+Inten2[iii,jjj,11]+Inten2[iii,jjj,18]+Inten2[iii,jjj,19]+Inten2[iii,jjj,20])
;                        I1b= (Inten2[iii,jjj,3]+Inten2[iii,jjj,4]+Inten2[iii,jjj,5]+Inten2[iii,jjj,12]+Inten2[iii,jjj,13]+Inten2[iii,jjj,14]+Inten2[iii,jjj,21]+Inten2[iii,jjj,22]+Inten2[iii,jjj,23])
;                        I2b= (Inten2[iii,jjj,6]+Inten2[iii,jjj,7]+Inten2[iii,jjj,8]+Inten2[iii,jjj,15]+Inten2[iii,jjj,16]+Inten2[iii,jjj,17]+Inten2[iii,jjj,24]+Inten2[iii,jjj,25]+Inten2[iii,jjj,26])
;                    END
;                    2: BEGIN
;                        I0a= (Inten1[iii,jjj,0]+Inten1[iii,jjj,1]+Inten1[iii,jjj,2]+Inten1[iii,jjj,3]+Inten1[iii,jjj,4]+Inten1[iii,jjj,5]+Inten1[iii,jjj,6]+Inten1[iii,jjj,7]+Inten1[iii,jjj,8])
;                        I2a= (Inten1[iii,jjj,9]+Inten1[iii,jjj,10]+Inten1[iii,jjj,11]+Inten1[iii,jjj,12]+Inten1[iii,jjj,13]+Inten1[iii,jjj,14]+Inten1[iii,jjj,15]+Inten1[iii,jjj,16]+Inten1[iii,jjj,17])
;                        I1a= (Inten1[iii,jjj,18]+Inten1[iii,jjj,19]+Inten1[iii,jjj,20]+Inten1[iii,jjj,21]+Inten1[iii,jjj,22]+Inten1[iii,jjj,23]+Inten1[iii,jjj,24]+Inten1[iii,jjj,25]+Inten1[iii,jjj,26])
;                        I0b= (Inten2[iii,jjj,0]+Inten2[iii,jjj,1]+Inten2[iii,jjj,2]+Inten2[iii,jjj,3]+Inten2[iii,jjj,4]+Inten2[iii,jjj,5]+Inten2[iii,jjj,6]+Inten2[iii,jjj,7]+Inten2[iii,jjj,8])
;                        I2b= (Inten2[iii,jjj,9]+Inten2[iii,jjj,10]+Inten2[iii,jjj,11]+Inten2[iii,jjj,12]+Inten2[iii,jjj,13]+Inten2[iii,jjj,14]+Inten2[iii,jjj,15]+Inten2[iii,jjj,16]+Inten2[iii,jjj,17])
;                        I1b= (Inten2[iii,jjj,18]+Inten2[iii,jjj,19]+Inten2[iii,jjj,20]+Inten2[iii,jjj,21]+Inten2[iii,jjj,22]+Inten2[iii,jjj,23]+Inten2[iii,jjj,24]+Inten2[iii,jjj,25]+Inten2[iii,jjj,26])
;
;                    END
;                ENDCASE

                                ; FOR THE CASE WHERE INTEN[*,*,5] AND
                                ; INTEN[*,*,21] ARE CRAP
                CASE i OF
                    0: BEGIN ; CHANGED FOR [-1.5,1.5,1.5]
                        I0a= (Inten1[iii,jjj,0]+Inten1[iii,jjj,6]+Inten1[iii,jjj,9]+Inten1[iii,jjj,12]+Inten1[iii,jjj,15]+Inten1[iii,jjj,18]+Inten1[iii,jjj,24])
                        I1a = (Inten1[iii,jjj,1]+Inten1[iii,jjj,7]+Inten1[iii,jjj,10]+Inten1[iii,jjj,13]+Inten1[iii,jjj,16]+Inten1[iii,jjj,19]+Inten1[iii,jjj,25])
                        I2a = (Inten1[iii,jjj,2]+Inten1[iii,jjj,8]+Inten1[iii,jjj,11]+Inten1[iii,jjj,14]+Inten1[iii,jjj,17]+Inten1[iii,jjj,20]+Inten1[iii,jjj,26])
                        I0b= (Inten2[iii,jjj,0]+Inten2[iii,jjj,3]+Inten2[iii,jjj,6]+Inten2[iii,jjj,9]+Inten2[iii,jjj,12]+Inten2[iii,jjj,15]+Inten2[iii,jjj,18]+Inten2[iii,jjj,21]+Inten2[iii,jjj,24])
                        I1b = (Inten2[iii,jjj,1]+Inten2[iii,jjj,4]+Inten2[iii,jjj,7]+Inten2[iii,jjj,10]+Inten2[iii,jjj,13]+Inten2[iii,jjj,16]+Inten2[iii,jjj,19]+Inten2[iii,jjj,22]+Inten2[iii,jjj,25])
                        I2b = (Inten2[iii,jjj,2]+Inten2[iii,jjj,5]+Inten2[iii,jjj,8]+Inten2[iii,jjj,11]+Inten2[iii,jjj,14]+Inten2[iii,jjj,17]+Inten2[iii,jjj,20]+Inten2[iii,jjj,23]+Inten2[iii,jjj,26])
                    END
                    1: BEGIN
                        I0a= (Inten1[iii,jjj,0]+Inten1[iii,jjj,1]+Inten1[iii,jjj,9]+Inten1[iii,jjj,10]+Inten1[iii,jjj,11]+Inten1[iii,jjj,19]+Inten1[iii,jjj,20])
                        I1a= (Inten1[iii,jjj,3]+Inten1[iii,jjj,4]+Inten1[iii,jjj,12]+Inten1[iii,jjj,13]+Inten1[iii,jjj,14]+Inten1[iii,jjj,22]+Inten1[iii,jjj,23])
                        I2a= (Inten1[iii,jjj,6]+Inten1[iii,jjj,7]+Inten1[iii,jjj,15]+Inten1[iii,jjj,16]+Inten1[iii,jjj,17]+Inten1[iii,jjj,25]+Inten1[iii,jjj,26])
                        I0b= (Inten2[iii,jjj,0]+Inten2[iii,jjj,1]+Inten2[iii,jjj,2]+Inten2[iii,jjj,9]+Inten2[iii,jjj,10]+Inten2[iii,jjj,11]+Inten2[iii,jjj,18]+Inten2[iii,jjj,19]+Inten2[iii,jjj,20])
                        I1b= (Inten2[iii,jjj,3]+Inten2[iii,jjj,4]+Inten2[iii,jjj,5]+Inten2[iii,jjj,12]+Inten2[iii,jjj,13]+Inten2[iii,jjj,14]+Inten2[iii,jjj,21]+Inten2[iii,jjj,22]+Inten2[iii,jjj,23])
                        I2b= (Inten2[iii,jjj,6]+Inten2[iii,jjj,7]+Inten2[iii,jjj,8]+Inten2[iii,jjj,15]+Inten2[iii,jjj,16]+Inten2[iii,jjj,17]+Inten2[iii,jjj,24]+Inten2[iii,jjj,25]+Inten2[iii,jjj,26])
                    END
                    2: BEGIN
                        I0a= (Inten1[iii,jjj,0]+Inten1[iii,jjj,1]+Inten1[iii,jjj,2]+Inten1[iii,jjj,4]+Inten1[iii,jjj,6]+Inten1[iii,jjj,7]+Inten1[iii,jjj,8])
                        I2a= (Inten1[iii,jjj,9]+Inten1[iii,jjj,10]+Inten1[iii,jjj,11]+Inten1[iii,jjj,13]+Inten1[iii,jjj,15]+Inten1[iii,jjj,16]+Inten1[iii,jjj,17])
                        I1a= (Inten1[iii,jjj,18]+Inten1[iii,jjj,19]+Inten1[iii,jjj,20]+Inten1[iii,jjj,22]+Inten1[iii,jjj,24]+Inten1[iii,jjj,25]+Inten1[iii,jjj,26])
                        I0b= (Inten2[iii,jjj,0]+Inten2[iii,jjj,1]+Inten2[iii,jjj,2]+Inten2[iii,jjj,3]+Inten2[iii,jjj,4]+Inten2[iii,jjj,5]+Inten2[iii,jjj,6]+Inten2[iii,jjj,7]+Inten2[iii,jjj,8])
                        I2b= (Inten2[iii,jjj,9]+Inten2[iii,jjj,10]+Inten2[iii,jjj,11]+Inten2[iii,jjj,12]+Inten2[iii,jjj,13]+Inten2[iii,jjj,14]+Inten2[iii,jjj,15]+Inten2[iii,jjj,16]+Inten2[iii,jjj,17])
                        I1b= (Inten2[iii,jjj,18]+Inten2[iii,jjj,19]+Inten2[iii,jjj,20]+Inten2[iii,jjj,21]+Inten2[iii,jjj,22]+Inten2[iii,jjj,23]+Inten2[iii,jjj,24]+Inten2[iii,jjj,25]+Inten2[iii,jjj,26])

                    END
                ENDCASE


                tota      = I0a+I1a+I2a
                totb      = I0b+I1b+I2b
                mat=REFORM(mat1[i,*,*]*tota/3.d0+mat2[i,*,*]*totb/3.d0)
                mat=LA_INVERT(mat,/DOUBLE)
                vec[0,0]  = (I0a+I0b)-(tota+totb)/3.d0
                vec[0,1]  = (I1a+I1b)-(tota+totb)/3.d0
                vec       =  mat[*,*]##vec
                Bg0[iii,jjj,i]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)
              ; TO CLEAN THE RESULTS I FORCE THE CONTRASTS TO 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF(Bg0[iii,jjj,i] GT 1.0d0) THEN Bg0[iii,jjj,i] = 1.0d0
                Phig0[iii,jjj,i] = ATAN(vec[0,1],vec[0,0])
            ENDFOR
            
        ENDIF ELSE BEGIN

            Bg0[iii,jjj,*]  =1.d0
            Phig0[iii,jjj,*]=0.d0

        ENDELSE
        

        FOR j=0,nseq-1 DO BEGIN
            profile = [blocker01,blocker02]
            FOR i=0,2 DO profile = profile * (1.d0+Bg0[iii,jjj,i]*COS(2.d0*!dpi/FSR[i]*[lam1,lam2]+Phig0[iii,jjj,i]+tuning[i,j]))/2.d0 ; non-tunable part
            FOR i=3,6 DO profile = profile * (1.d0+0.99d0*COS(2.d0*!dpi/FSR[i]*[lam1,lam2]))/2.d0 ; tunable part
            Intenrecons1[iii,jjj,j] = profile[0]
            Intenrecons2[iii,jjj,j] = profile[1]
        ENDFOR
        
    ENDFOR
ENDFOR

temp = REFORM(phig0[*,*,0])
a = WHERE(abs(temp*180.d0/!dpi) gt 100,n)
temp[a] = 0.d0
phig0[*,*,0]=temp[*,*]



SAVE,Bg0,Phig0,FILE='SUNTEST_LASER1_BIS.BIN'
PRINT,'TOTAL TIME ELAPSED',SYSTIME(1)-time0

SET_PLOT,'ps'
!p.multi=[0,1,2]
device,file='yo.ps',bits=24,xoffset=0,yoffset=1,xsize=15,ysize=22,/color
loadct,4
tvim,phig0[*,*,0]*180.d0/!dpi,/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)'
tvim,phig0[*,*,1]*180.d0/!dpi,/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)'
tvim,phig0[*,*,2]*180.d0/!dpi,/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)'
tvim,Bg0[*,*,0],/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast'
tvim,Bg0[*,*,1],/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast'
tvim,Bg0[*,*,2],/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast'
moy=DBLARR(3)

;FOR i=0,2 DO moy[i]=MEAN(phig0[880/(4096./nx):2320/(4096./nx),2080/(4096./nx):3120/(4096./nx),i])*180.d0/!dpi
FOR i=0,2 DO moy[i]=MEAN(phig0[150:300,250:350,i])*180.d0/!dpi

tvim,phig0[*,*,0]*180.d0/!dpi-moy[0],/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase after average removed',range=[-30,35]
tvim,phig0[*,*,1]*180.d0/!dpi-moy[1],/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase after average removed',range=[-50,40]
tvim,phig0[*,*,2]*180.d0/!dpi-moy[2],/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase after average removed',range=[-50,60]
DEVICE,/close
PRINT,'AVERAGES'
FOR i=0,2 DO PRINT,moy[i]

!p.multi=[0,1,2]
device,file='yo2.ps',bits=24,xoffset=0,yoffset=0,xsize=15,ysize=22
moncul=DBLARR(27)

for i=0,26 do moncul[i]=MEAN(intenrecons1[1200/(4096./nx):1760/(4096./nx),2400/(4096./nx):3040/(4096./nx),i])
moncul[5]=0.d0
moncul[21]=0.d0
plot,moncul/MAX(moncul),tit='Intensities',xst=1,ytit='Norm. Int. for lam1',xtit='Position in the detune sequence',yst=1
for i=0,26 do moncul[i]=MEAN(inten1[1200/(4096./nx):1760/(4096./nx),2400/(4096./nx):3040/(4096./nx),i])
oplot,moncul/MAX(moncul),linestyle=2

for i=0,26 do moncul[i]=MEAN(intenrecons2[1200/(4096./nx):1760/(4096./nx),2400/(4096./nx):3040/(4096./nx),i])
plot,moncul/MAX(moncul),tit='Reconstructed Intensities',xst=1,ytit='Norm. Int. for lam2',xtit='Position in the detune sequence',yst=1
for i=0,26 do moncul[i]=MEAN(inten2[1200/(4096./nx):1760/(4096./nx),2400/(4096./nx):3040/(4096./nx),i])
oplot,moncul/MAX(moncul),linestyle=2
device,/close

; FIND THE MOTOR POSITIONS FOR THE CO-TUNE SEQUENCE

; DEFINITION OF A CO-TUNE IN RADIANS
!p.multi=0
nseq2=20 ; NUMBER OF POSITION IN THE CO-TUNE
tuning2=DBLARR(3,nseq2)
;FOR i=0,nseq2-1 DO tuning2[*,i] = [ (DOUBLE(i)/DOUBLE(nseq2)*2.d0*!dpi-!dpi)*FSR[2]/FSR[0],(DOUBLE(i)/DOUBLE(nseq2)*2.d0*!dpi-!dpi)*FSR[2]/FSR[1],(DOUBLE(i)/DOUBLE(nseq2)*2.d0*!dpi-!dpi)*FSR[2]/FSR[2]]

tuning2[*,0] =[ -30,   0,   0]
tuning2[*,1] =[ -27,  -6, -12]
tuning2[*,2] =[ -24, -12, -24]
tuning2[*,3] =[ -21, -18,  24]
tuning2[*,4] =[ -18, -24,  12]
tuning2[*,5] =[ -15, -30,   0]
tuning2[*,6] =[ -12,  24, -12]
tuning2[*,7] =[  -9,  18, -24]
tuning2[*,8] =[  -6,  12,  24]
tuning2[*,9] =[  -3,   6,  12]
tuning2[*,10]=[  0 ,  0 ,   0]
tuning2[*,11]=[  3 , -6 , -12]
tuning2[*,12]=[  6 ,-12 , -24]
tuning2[*,13]=[  9 ,-18 ,  24]
tuning2[*,14]=[ 12 ,-24 ,  12]
tuning2[*,15]=[ 15 ,-30 ,   0]
tuning2[*,16]=[ 18 , 24 , -12]
tuning2[*,17]=[ 21 , 18 , -24]
tuning2[*,18]=[ 24 , 12 ,  24]
tuning2[*,19]=[ 27 ,  6 ,  12]
FOR i=0,nseq2-1 DO tuning2[*,i] = tuning2[*,i] + [111,100,67]

tuning2=REVERSE(tuning2,1)
FOR i=0,nseq2-1 DO tuning2[*,i] = REFORM(tuning2[*,i])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]


nlam        = 5000    ; number of wavelengths            
dlam        = 2.d0/1.d3        ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
SET_PLOT,'PS'
DEVICE,file='yo3.ps',bits=24,/color
loadct,4
FOR i=0,nseq2-1 DO BEGIN & profile=DBLARR(nlam)+1.d0 & FOR j=0,2 DO profile = profile * (1.d0+COS(2.d0*!dpi/FSR[j]*lam+tuning2[j,i]+moy[j]*!dpi/180.d0))/2.d0 & FOR j=3,6 DO profile = profile * (1.d0+COS(2.d0*!dpi/FSR[j]*lam))/2.d0 & IF i EQ 0 THEN plot,lam,profile,yrange=[0,1],yst=1,xrange=[-0.6,0.6],xst=1,col=i*254/(nseq2-1),thick=3 ELSE oplot,lam,profile,col=i*254/(nseq2-1) & ENDFOR
DEVICE,/close


SET_PLOT,'x'
!p.multi=0
READ,pause

END
