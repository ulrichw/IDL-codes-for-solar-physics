; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE LAMP AS A SOURCE, AND DE-TUNE SEQUENCE
; TO DERIVE THE FIXED TRANSMISSION PROFILE IN THE RANGE [-2 FSR[0],+2 FSR[0]]
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_lamp,draw

lamref      = 6173.3433;reference central wavelength for the FeI 6173 line
                      ; IN AIR
                      ; used to center the averaged blocker
                      ; and front window transmission profiles
                      ; value from Stenflo & Lindegren (1977)
                      ; NB: IF THIS VALUE IS NOT THE ACTUAL ONE, THEN
                      ; WE HAVE AN AVERAGE VELOCITY ON THE SOLAR DISK
                      ; DIFFERENT FROM 0. TO BE CORRECTED


time0       = SYSTIME(1)
dpi         = 2.d0*!dpi

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 128              ; number of rows
ny          = 128              ; number of columns
mu0         = 0.d0;0.013       ; regularization parameter for the least-squares fit. TO ADJUST DEPENDING ON THE NOISE LEVEL
nparam      = 15               ; number of parameters to fit for
nseq        = 27               ; number of positions in the detune sequence
; !!! WARNING: if nseq is not equal to 27, then the guess value of I0g is
; not correct and the code must be modified
nlam        = 2250             ; number of wavelengths           
dlam        = 3.6d0/1.0d3      ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelength values
Inten       = DBLARR(nx,ny,nseq) ; measured output intensities (measured on a HMI CCD)
Inten0      = Inten            ; output intensities with a perfect blocker (no fringes)
Inten0_res  = Inten            ; output intensities in the wavelength range [-2 FSR[0],+2 FSR[0]]
tuning      = DBLARR(3,nseq)   ; tuning positions for the detune sequence
Residual    = DBLARR(nseq)     ; residual of the least-squares fit
Jac         = DBLARR(nparam,nseq);Jacobian matrix of the least-squares fit
proxy1      = DBLARR(nx,ny,(nparam-1)/2);coefficients of the Fourier series expansion
proxy2      = DBLARR(nx,ny,(nparam-1)/2)
proxy0      = DBLARR(nx,ny)
center      = DBLARR(nx,ny)    ; central wavelength of the fixed transmission profile
I0g         = DBLARR(nx,ny)
dIdproxy1   = DBLARR((nparam-1)/2)
dIdproxy2   = DBLARR((nparam-1)/2)    
maxsteps    = 40;100
history     = DBLARR(maxsteps,nparam)
err         = DBLARR(maxsteps)

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; actual FSRs
; e-mail 10/25/2005
; e-mail 11/18/2005
; NB: the code is formally written for FSRs that are multiples of the FSR
; of the narrow Michelson
;-----------------------------------


FSR     = DBLARR(7)    ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]  = 0.1724570    ; for the narrow-band Michelson in Angstrom
FSR[1]  = 0.344242     ; for the broad-band  Michelson
FSR[2]  = 0.702        ; for E1
FSR[3]  = 1.405        ; for E2
FSR[4]  = 2.779        ; for E3
FSR[5]  = 5.682        ; for E4
FSR[6]  = 11.354       ; for E5

; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
;                NB MICHELSON   WB MICHELSON         E1 LYOT
;-----------------------------------------------------------

tuning[*,0]  = [         0.d0,          0.d0 ,         0.d0]
tuning[*,1]  = [        80.d0,          0.d0 ,         0.d0]
tuning[*,2]  = [       160.d0,          0.d0 ,         0.d0]
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


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;-----------------------------------------------

; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

;blocker0    = DBLARR(nlam)+1.d0
q            = READFITS('blocker11.fits')           ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam);,/LSQUADRATIC); I center the profile


; CONTRASTS AND PHASES OBTAINED BY HMI_laserX.pro and HMI_sunX.pro
;---------------------------------------------------------------

;RESTORE,file='SUNTEST_LASER2.BIN' ; Bg and Phig !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;Bg   = DBLARR(128,128,7)+1.d0
;Phig = DBLARR(128,128,7);+30*!dpi/180.d0  ; CAREFUL WITH THE UNITS OF Phig

;RESTORE,'RESULTS/RESULTS_listSun060224_225806_256.BIN'; phases in
;calmode from sunlight: contrast = 1 for non-tunable part, = 0.98 for
;tunable part
RESTORE,'RESULTS/RESULTS_listSun060224_225806_128.BIN' ; idem but more recent and better dark removal, and contrast=0.98 for everything
Phig0= REBIN(Phig,nx,nx,3)/ 180.d0 * !dpi ; convert from degrees to radians
Phig = DBLARR(128,128,7)
Phig[*,*,0:2] = Phig0
RESTORE,'RESULTS/RESULTS_DETUNE180_list_detune180_060303_003721_128.BIN'
Bg   = DBLARR(nx,nx,7)
Bg[*,*,0:2] = Bg0[*,*,0:2,0]   ; NB: write here what is used in HMI_sun1.pro
Bg[*,*,3:6] = 0.98d0  ;1.d0    ; NB: write here what is used in HMI_sun1.pro


; MEASURED INTENSITIES
;-------------------------------------------------------

list = 'list_060304_031619'

IF draw EQ 1 THEN GOTO,draw

;GOTO,dejalu
;readimages,list,images,headers,time,29,nbin=nx,iover=iover
;SAVE,images,time,FILE='SEQUENCE0_DETUNE27_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
RESTORE,'SEQUENCE0_DETUNE27_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

images[*,*,21]     = images[*,*,21]/4502.34d0*3824.67d0 ; attempt to correct for bad exposure

corner             = REBIN(images(0:nx/64-1,0:nx/64-1,*),1,1,29)
FOR i              = 0,28 DO images[*,*,i] = images[*,*,i]-corner[i]

Darks              = FLTARR(nx,ny,2)
Darks[*,*,0]       = images[*,*,0]
Darks[*,*,1]       = images[*,*,28]
Inten              = FLTARR(nx,nx,nseq)
Inten[*,*,*]       = images[*,*,1:27]
; WE GET RID OF THE BAD EXPOSURES
;Inten[*,*,0:19]    = images[*,*,1:20]
;Inten[*,*,20:25]   = images[*,*,22:27]
timeimages         = FLTARR(nseq)
timeimages         = time[1:27]
;timeimages[0:19]   = time[1:20]
;timeimages[20:25]  = time[22:27]
timedarks          = [time[0],time[28]]


; WE REMOVE THE DARK FRAME
PRINT,'DARK FRAME REMOVAL'
    FOR j=0,ny-1 DO BEGIN
        FOR jj=0,nx-1 DO BEGIN
            Inten[jj,j,*] = Inten[jj,j,*] - INTERPOL(REFORM(Darks[jj,j,*]),timedarks,timeimages[*])
        ENDFOR
    ENDFOR
a    = WHERE(Inten LT 0.d0)   
IF (a[0] NE -1) THEN Inten[a]= 0.d0

SAVE,Inten,FILE='SEQUENCE_DETUNE27_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

dejalu:
RESTORE,'SEQUENCE_DETUNE27_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

; GUESS PARAMETERS
;---------------------------

proxy1[*,*,*] = 0.005d0 ; COSINE coefficients of the Fourier series expansion
proxy2[*,*,*] = 0.005d0 ;   SINE coefficients of the Fourier series expansion
        
WINDOW,0,retain=2,xsize=700,ysize=900
!p.multi=[0,1,2]

distance = SHIFT(DIST(nx,ny),nx/2,ny/2)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


FOR jjj=0,ny-1 DO BEGIN

    PRINT,jjj
    TVIM,proxy0,/scale
    TVIM,I0g,/scale

    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] LE 1024.d0) THEN BEGIN  ; XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

            profileg0 = blocker0
            FOR i=3,6 DO profileg0  = profileg0 * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam+Phig[iii,jjj,i]) )/2.d0
        
            I0g[iii,jjj] = 8.d0/DOUBLE(nseq)*TOTAL(Inten[iii,jjj,*])/TOTAL(profileg0)/dlam ; GUESS OF I0g because cos(a)+cos(a+120)+cos(a+240)=0

          ; INTENSITIES IN THE RANGE -2 FSR[0],+2 FSR[0] FOR A PERFECT BLOCKER
          ; ---------------------------------------------------------------

          ; we want to center the [-2 FSR[0],+2 FSR[0]] interview
          ; WE COMPUTE AN ESTIMATE OF THE FIXED TRANSMISSION PROFILE
            
            zero            = WHERE(profileg0 EQ MAX(profileg0))
            center[iii,jjj] = lam[zero[0]]
            lam2            = lam-center[iii,jjj]
        
         ; the Fourier series expansion supposes that the blocker+front
         ; window + non tunable Lyot profile is periodic with a period of 4*FSR[0]
            a = WHERE(ABS(lam2) LE 2.d0*FSR[0])
            na= DOUBLE(N_ELEMENTS(a))
            
            FOR j=0,nseq-1 DO BEGIN
                profile = profileg0
                FOR i=0,2 DO profile = profile * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam+Phig[iii,jjj,i]+tuning[i,j]))/2.d0
                Inten0[j]     = TOTAL(profile[*])*dlam;*I0g[iii,jjj]
                Inten0_res[j] = TOTAL(profile[a])*dlam;*I0g[iii,jjj]
            ENDFOR

         ; WE CORRECT THE INTENSITIES
         ; to have an approximate value of the intensities in the range [-2 FSR[0],+2 FSR[0]]
         ;-------------------------------------------------------------------------------

            Inten[iii,jjj,*]= Inten[iii,jjj,*] * Inten0_res[*]/Inten0[*]
;           Inten[iii,jjj,*]= Inten[iii,jjj,*] -(Inten0[*]-Inten0_res[*])
            
            converg = 0
            jj      = 0
           ;mu      = mu0
            minimum = 1.d10

            int_profileg0  = TOTAL(profileg0[a])/DOUBLE(na) ; Fourier coefficient a0          
   
            WHILE (converg EQ 0) DO BEGIN
                
                FOR j=0,nseq-1   DO BEGIN
                    
                  ; Fourier series of the blocker+front window+non-tunable Lyot transmission profile
                    profileg = DBLARR(na) + int_profileg0
                    FOR i=1,(nparam-1)/2  DO profileg = profileg + proxy1[iii,jjj,i-1]*COS(dpi*lam[a]*i/4.d0/FSR[0])+proxy2[iii,jjj,i-1]*SIN(dpi*lam[a]*i/4.d0/FSR[0])
                    
                    profileg1  = DBLARR(na) + 1.d0
                    FOR i=0,2 DO profileg1  = profileg1 * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam[a]+Phig[iii,jjj,i]+tuning[i,j]))/2.d0
                    profileg   = profileg * profileg1
                    
                    Residual[j] = Inten[iii,jjj,j] - TOTAL(profileg)*dlam*I0g[iii,jjj]
                    
                    FOR i=1,(nparam-1)/2 DO BEGIN
                        dIdproxy1[i-1] = TOTAL(profileg1 *  COS(dpi*lam[a]*i/4.d0/FSR[0])  )*dlam*I0g[iii,jjj]
                        dIdproxy2[i-1] = TOTAL(profileg1 *  SIN(dpi*lam[a]*i/4.d0/FSR[0])  )*dlam*I0g[iii,jjj]
                    ENDFOR
                    dIdproxy0 = TOTAL(profileg1)*dlam*I0g[iii,jjj]
                   ;dIdI0     = TOTAL(profileg) *dlam

                    FOR i=0       ,(nparam-1)/2-1 DO Jac[i,j] = dIdproxy1[i]
                    FOR i=nparam/2, nparam-2      DO Jac[i,j] = dIdproxy2[i-nparam/2]
                    Jac[nparam-1,j] = dIdproxy0
                   ;Jac[nparam,j]   = dIdI0

                ENDFOR
                
                LA_SVD,Jac,W,U,V,/DOUBLE
               ;filter = W*W/(W*W+mu^2.d0) ; mu is the regularization parameter
               ;Dx     = V##DIAG_MATRIX(1.d0/W*filter)##TRANSPOSE(U)##Residual
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual

                temp   = TRANSPOSE(Dx)/[REFORM(proxy1[iii,jjj,*]),REFORM(proxy2[iii,jjj,*]),int_profileg0];,I0g[iii,jjj]]
               ;err = TOTAL(err^2.d0)/DOUBLE(nparam)
                err[jj]= MAX(ABS(temp))
                
                ;IF(jj GE maxsteps-20) THEN IF err LT minimum THEN minimum = err ELSE converg = 2
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7)  THEN converg = 2
                ;IF(err LT  5.d-6)                     THEN mu      = 0.8*mu
                ;IF(err GE  5.d-6 AND err LT 1.d- 1)   THEN mu      = mu0 
                ;IF(FINITE(err) EQ 0 AND jj GT 0 )     THEN mu      = 2.5d0*mu
                IF(err[jj] LE 1.d-7)                       THEN converg = 1
                

                IF(converg EQ 2) THEN BEGIN
                    j       = WHERE(err EQ MIN(err))
                    FOR i=0,(nparam-1)/2-1 DO BEGIN
                        proxy1[iii,jjj,i] = history[j[0],i+1]
                        proxy2[iii,jjj,i] = history[j[0],i+(nparam-1)/2+1]
                    ENDFOR
                    int_profileg0 = history[j[0],0]
                   ;I0g[iii,jjj]  = history[j[0],nparam]
                ENDIF



                IF(converg NE 2) THEN BEGIN
                    FOR i=0,(nparam-1)/2-1 DO BEGIN
                        proxy1[iii,jjj,i] = proxy1[iii,jjj,i]+Dx[i]
                        IF(ABS(proxy1[iii,jjj,i]) GT 1.d0) THEN proxy1[iii,jjj,i] = 0.005d0
                        proxy2[iii,jjj,i] = proxy2[iii,jjj,i]+Dx[i+(nparam-1)/2]
                        IF(ABS(proxy2[iii,jjj,i]) GT 1.d0) THEN proxy2[iii,jjj,i] = 0.005d0
                    ENDFOR
                    int_profileg0 = int_profileg0 + Dx[nparam-1]
                    IF(int_profileg0 GT 1.d0 OR int_profileg0 LT 0.d0) THEN int_profileg0 = TOTAL(profileg0[a])/DOUBLE(na)
                   ;I0g[iii,jjj]  = I0g[iii,jjj]  + Dx[nparam]
               ENDIF

                history[jj,*] = [int_profileg0,REFORM(proxy1[iii,jjj,0:(nparam-1)/2-1]),REFORM(proxy2[iii,jjj,0:(nparam-1)/2-1])];,I0g[iii,jjj]]

                ;IF(jj EQ 50 OR jj EQ 100 OR jj EQ 150 AND converg NE 1) THEN BEGIN
                ;    proxy1[iii,jjj,*] = RANDOMU(seed,(nparam-1)/2)*0.01d0-0.005
                ;    proxy2[iii,jjj,*] = RANDOMU(seed,(nparam-1)/2)*0.01d0-0.005
                ;ENDIF
                
               ;PRINT,jj,err
               jj = jj+1
                
                
            ENDWHILE
            
            ;PRINT,'parameters='
            ;FOR i=0,(nparam-1)/2-1 DO PRINT,'Amp0=',proxy1[iii,jjj,i]
            ;FOR i=0,(nparam-1)/2-1 DO PRINT,'Amp1=',proxy2[iii,jjj,i]
            ;PRINT,'A0=',int_profileg0
            proxy0[iii,jjj] = int_profileg0
            ;PRINT,'I0=',I0g[iii,jjj]

            ;profileg = DBLARR(nlam)+proxy0[iii,jjj] ; reconstructed non-tunable profile
            ;FOR i=1,(nparam-1)/2  DO profileg = profileg + proxy1[iii,jjj,i-1]*COS(dpi*lam*i/4.d0/FSR[0])+proxy2[iii,jjj,i-1]*SIN(dpi*lam*i/4.d0/FSR[0])
            ;plot,lam,profileg,xst=1
            ;oplot,[center[iii,jjj]-2.d0*FSR[0],center[iii,jjj]-2.d0*FSR[0]],[0,1],line=2
            ;oplot,[center[iii,jjj]+2.d0*FSR[0],center[iii,jjj]+2.d0*FSR[0]],[0,1],line=2

        ENDIF
            
    ENDFOR
ENDFOR

SAVE,proxy1,proxy2,proxy0,center,I0g,FILE='RESULTS/RESULTS_LAMP_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'TIME=',SYSTIME(1)-time0

draw:

RESTORE,'RESULTS/RESULTS_LAMP_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
distance = SHIFT(DIST(nx,ny),nx/2,ny/2)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center

;a = WHERE(proxy0 GT 0.431 AND proxy0 LT 0.4321 AND distance LT 600) ;WARNING distance!!!!!!!!!!!!!!
;proxy0 = MEAN(proxy0[a])
;tempa  = DBLARR(7)
;tempb  = DBLARR(7)
;FOR i=0,6 DO BEGIN
;    temp = REFORM(proxy1[*,*,i])
;    tempa[i] = MEAN(temp[a])
;    temp = REFORM(proxy2[*,*,i])
;    tempb[i] = MEAN(temp[a])
;ENDFOR

;profileg0 = blocker0
;FOR i=3,6 DO profileg0  = profileg0 * (1.d0+MEAN(Bg[*,*,i])*COS(dpi/FSR[i]*lam) )/2.d0 ;!!!!WARNING CONTRASTS AND PHASES
;profileg = DBLARR(nlam)+proxy0 ; reconstructed non-tunable profile
;FOR i=1,(nparam-1)/2  DO profileg = profileg + tempa[i-1]*COS(dpi*lam*i/4.d0/FSR[0])+tempb[i-1]*SIN(dpi*lam*i/4.d0/FSR[0])

profileg0 = blocker0
FOR i=3,6 DO profileg0  = profileg0 * (1.d0+MEAN(Bg[*,*,i])*COS(dpi/FSR[i]*lam) )/2.d0
profileg  = DBLARR(nlam)
ncount    = 0.d0
FOR iii=0,nx-1 DO BEGIN
    FOR jjj=0,nx-1 DO BEGIN
        IF(distance[iii,jjj] LT 750) THEN BEGIN
            FOR i=1,(nparam-1)/2  DO profileg = profileg + proxy1[iii,jjj,i-1]*COS(dpi*lam*i/4.d0/FSR[0])+proxy2[iii,jjj,i-1]*SIN(dpi*lam*i/4.d0/FSR[0])
            profileg = profileg + proxy0[iii,jjj]
            ncount = ncount+1.d0
        ENDIF
    ENDFOR
ENDFOR
profileg = profileg/ncount


SET_PLOT,'ps'
!P.MULTI=0
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=15,ysize=10
plot,lam,profileg/profileg0-1.0,xrange=[-2.d0*FSR[0],2.d0*FSR[0]],xst=1,charsize=1.5
DEVICE,/close

SET_PLOT,'x'

READ,pause

END
