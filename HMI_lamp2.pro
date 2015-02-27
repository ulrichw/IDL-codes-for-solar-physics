; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE LAMP AS A SOURCE, AND DE-TUNE SEQUENCE
; TO DERIVE THE FIXED TRANSMISSION PROFILE IN THE RANGE [-2 FSR[0],+2
; FSR[0]]

; UNLIKE HMI_lamp.pro, WE FIT 3 FREQUENCIES RELATED TO THE FRONT WINDOW
; THICKNESSES INSTEAD OF 7

;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_lamp2,draw

lamref      = 6173.3433d0;reference central wavelength for the FeI 6173 line
                      ; IN AIR
                      ; used to center the averaged blocker
                      ; and front window transmission profiles
                      ; value from Stenflo & Lindegren (1977)
                      ; NB: IF THIS VALUE IS NOT THE ACTUAL ONE, THEN
                      ; WE HAVE AN AVERAGE VELOCITY ON THE SOLAR DISK
                      ; DIFFERENT FROM 0. TO BE CORRECTED


time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
anglim      = 950.

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 256              ; number of rows
ny          = 256              ; number of columns
mu0         = 0.d0;0.013       ; regularization parameter for the least-squares fit. TO ADJUST DEPENDING ON THE NOISE LEVEL
nparam      = 6                ; number of parameters to fit for
nseq        = 27               ; number of positions in the detune sequence
nlam        = 2250             ; number of wavelengths           
dlam        = 3.6d0/1.0d3      ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelength values
Inten       = DBLARR(nx,ny,nseq) ; measured output intensities (measured on a HMI CCD)
tuning      = DBLARR(3,nseq)   ; tuning positions for the detune sequence
Residual    = DBLARR(nseq)     ; residual of the least-squares fit
Jac         = DBLARR(nparam+1,nseq);Jacobian matrix of the least-squares fit
proxy1      = DBLARR(nx,ny,nparam/2);coefficients of the Fourier series expansion
proxy2      = DBLARR(nx,ny,nparam/2)
I0g         = DBLARR(nx,ny)   
dIdp1       = DBLARR(nparam/2,nseq)
dIdp2       = DBLARR(nparam/2,nseq)
maxsteps    = 100
history     = DBLARR(maxsteps,nparam+1)
err         = DBLARR(maxsteps)
freqs       = DBLARR(nparam/2)

; NB: with 4 freqs instead of 3, the code does not converge. It's
; because there is nothing at 0.015d0
;freqs      = ([.003d0,.006d0,.009d0,.015d0]*1.d10*2.d0*1.516d0+lamref)/lamref^2.d0 ; in A^{-1} DUE TO FRONT WINDOW
freqs       = ([.003d0,.006d0,.009d0]*1.d10*2.d0*1.516d0+lamref)/lamref^2.d0 ; in A^{-1}

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; actual FSRs
; e-mail 10/25/2005
; e-mail 11/18/2005
; NB: the code is formally written for FSRs that are multiples of the FSR
; of the narrow Michelson
;-----------------------------------


FSR     = DBLARR(7)    ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
; FSRs CHANGED IN JUNE 2009
FSR[0]      = 0.17094240d0;0.1710098d0            ; for the narrow-band Michelson
FSR[1]      = 0.34192317d0;0.3421506d0            ; for the broad-band  Michelson
FSR[2]      = 0.69348347d0;0.6943613d0            ; for E1
FSR[3]      = 1.407d0;*0.993d0                    ; for E2
FSR[4]      = 2.779d0;*0.993d0                    ; for E3
FSR[5]      = 5.682d0;*0.993d0                    ; for E4
FSR[6]      = 11.354d0;*0.993d0                   ; for E5


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
tuning[*,9 ] = [         0.d0,          0.d0,         80.d0]
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

transmission = DBLARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength   = REFORM(transmission[0,*])
transmission = REFORM(transmission[1,*])
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-lamref+2.6d0,lam)

; CONTRASTS AND PHASES OBTAINED BY HMI_laserX.pro and HMI_sunX.pro
;---------------------------------------------------------------

RESTORE,'CPT/CPT_laser_551228_phasesmodified.bin' ;Obsmode

Phig1 = DBLARR(nx,nx,7)
Phig1[*,*,0:2] = Phig0
Bg1   = DBLARR(nx,nx,7)
Bg1[*,*,0:2] = Bg0

RESTORE,'RESULTS/RESULTS_June09_710660_CAL_256_SIDE.BIN' ;NON-TUNABLE ELEMENTS PHASE+CONTRASTS
Phig1[*,*,3:6]=Phig[*,*,3:6]*!dpi/180.d0
Bg1[*,*,3:6]=Bg[*,*,3:6]

Phig=Phig1
Bg=Bg1

center   = [127.5,127.5]  ; !!!WARNING, depends on the images
distance=FLTARR(nx,nx)
for i=0,nx-1 do for j=0,nx-1 do distance[i,j]=SQRT((i-center[0])^2.d0 + (j-center[1])^2.d0)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center

; MEASURED INTENSITIES
;-------------------------------------------------------


IF draw EQ 1 THEN GOTO,draw

RESTORE,'./CPT/SEQUENCE_listSun070905_561378_256.BIN' ;CALMODE
Inten=Imf[*,*,2:28] ; we use front camera only
Imf=0.0

; GUESS PARAMETERS
;---------------------------

proxy1[*,*,*] = 0.005d0 ; COSINE coefficients of the Fourier series expansion
proxy2[*,*,*] = 0.005d0 ;   SINE coefficients of the Fourier series expansion
        
WINDOW,0,retain=2,xsize=1300,ysize=700

; BEGINNING OF THE FIT
;---------------------------

FOR jjj=0,ny-1 DO BEGIN

    PRINT,jjj
    !p.multi=[0,4,2]
    FOR i=0,nparam/2-1 DO BEGIN
        TVIM,proxy1[*,*,i],/scale
        TVIM,proxy2[*,*,i],/scale
    ENDFOR

    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] LE anglim) THEN BEGIN

            profileg0 = blocker0
            FOR i=3,6 DO profileg0  = profileg0 * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam+Phig[iii,jjj,i]) )/2.d0
            I0g[iii,jjj] = 8.d0/DOUBLE(nseq)*TOTAL(Inten[iii,jjj,*])/TOTAL(profileg0)/dlam ; GUESS OF I0g because cos(a)+cos(a+120)+cos(a+240)=0
            profileg1  = FLTARR(nlam,nseq)
            FOR j=0,nseq-1 DO profileg1[*,j]  = profileg0
            FOR j=0,nseq-1 DO BEGIN 
                FOR i=0,2  DO profileg1[*,j]  = profileg1[*,j] * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam+Phig[iii,jjj,i]+tuning[i,j]))/2.d0       
                FOR ii=0,nparam/2-1 DO BEGIN
                    dIdp1[ii,j] = TOTAL(profileg1[*,j] *  COS(dpi*lam*freqs[ii]) )
                    dIdp2[ii,j] = TOTAL(profileg1[*,j] *  SIN(dpi*lam*freqs[ii]) )
                ENDFOR
            ENDFOR



            converg = 0
            jj      = 0
            minimum = 1.d10
   
            WHILE (converg EQ 0) DO BEGIN
                
                FOR j=0,nseq-1   DO BEGIN
                    
                    profileg = DBLARR(nlam)+1.d0
                    FOR i=0,nparam/2-1  DO profileg = profileg + proxy1[iii,jjj,i]*COS(dpi*lam*freqs[i])+proxy2[iii,jjj,i]*SIN(dpi*lam*freqs[i])
                    
                    profileg   = profileg * profileg1[*,j]
                    temp       = TOTAL(profileg)*dlam
                    Residual[j]= Inten[iii,jjj,j] - temp*I0g[iii,jjj]
                   
                    FOR i=0       ,nparam/2-1     DO Jac[i,j] = dIdp1[i,j]*dlam*I0g[iii,jjj]
                    FOR i=nparam/2,nparam-1       DO Jac[i,j] = dIdp2[i-nparam/2,j]*dlam*I0g[iii,jjj]
                    Jac[nparam,j] = temp
                ENDFOR
                
                LA_SVD,Jac,W,U,V,/DOUBLE

                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual

                temp   = TRANSPOSE(Dx)/[REFORM(proxy1[iii,jjj,*]),REFORM(proxy2[iii,jjj,*]),I0g[iii,jjj]]
                err[jj]= MAX(ABS(temp))
                
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7)  THEN converg = 2
                IF(err[jj] LE 1.d-7)                       THEN converg = 1
                

                IF(converg EQ 2) THEN BEGIN
                    PRINT,'NO CONVERGENCE'
                    j       = WHERE(err EQ MIN(err))
                    FOR i=0,nparam/2-1 DO BEGIN
                        proxy1[iii,jjj,i] = history[j[0],i]
                        proxy2[iii,jjj,i] = history[j[0],i+nparam/2]
                    ENDFOR
                    I0g[iii,jjj]  = history[j[0],nparam]
                ENDIF



                IF(converg NE 2) THEN BEGIN
                    FOR i=0,nparam/2-1 DO BEGIN
                        proxy1[iii,jjj,i] = proxy1[iii,jjj,i]+Dx[i]
                        IF(ABS(proxy1[iii,jjj,i]) GT 1.d0) THEN proxy1[iii,jjj,i] = 0.005d0
                        proxy2[iii,jjj,i] = proxy2[iii,jjj,i]+Dx[i+nparam/2]
                        IF(ABS(proxy2[iii,jjj,i]) GT 1.d0) THEN proxy2[iii,jjj,i] = 0.005d0
                    ENDFOR
                    I0g[iii,jjj]  = I0g[iii,jjj]  + Dx[nparam]
               ENDIF

               history[jj,*] = [REFORM(proxy1[iii,jjj,0:nparam/2-1]),REFORM(proxy2[iii,jjj,0:nparam/2-1]),I0g[iii,jjj]]

               jj = jj+1
                
                
            ENDWHILE
            

        ENDIF
            
    ENDFOR
ENDFOR

;SAVE,proxy1,proxy2,center,I0g,FILE='RESULTS/RESULTS_LAMP2_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
SAVE,proxy1,proxy2,center,I0g,FILE='RESULTS/RESULTS_LAMP2_561378_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'TIME=',SYSTIME(1)-time0

draw:

RESTORE,'RESULTS/RESULTS_LAMP2_561378_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

aa = WHERE(distance LT anglim, COMPLEMENT=b)
distance[aa] = 1.d0
distance[b]  = 0.d0

SET_PLOT,'ps'
!p.multi=[0,2,3]
loadct,4
DEVICE,file='yo.ps',xoffset=-0.8,yoffset=0,xsize=22,ysize=27,/color,bits=24
FOR i=0,nparam/2-1 DO BEGIN
    TVIM,proxy1[*,*,i]*distance,range=[-0.025,0.02],/scale,tit=STRING(freqs[i])
    TVIM,proxy2[*,*,i]*distance,range=[-0.025,0.02],/scale,tit=STRING(freqs[i])
ENDFOR
DEVICE,/close

profileg0 = blocker0
Inten0    = DBLARR(nx,nx,nseq)
FOR i=3,6 DO profileg0  = profileg0 * (1.d0+MEAN(Bg[*,*,i])*COS(dpi/FSR[i]*lam) )/2.d0
FOR iii=0,nx-1 DO BEGIN
    PRINT,iii
    FOR jjj=0,nx-1 DO BEGIN
        IF(distance[iii,jjj] EQ 1.d0) THEN BEGIN
            profileg  = DBLARR(nlam)
            FOR i=0,nparam/2-1  DO profileg = profileg + proxy1[iii,jjj,i]*COS(dpi*lam*freqs[i])+proxy2[iii,jjj,i]*SIN(dpi*lam*freqs[i])
            profileg = profileg+1.d0
            profileg = profileg * profileg0
            FOR i=0,nseq-1 DO Inten0[iii,jjj,i] = TOTAL(profileg * (1.d0+Bg[iii,jjj,0]*COS(dpi/FSR[0]*lam+Phig[iii,jjj,0]+tuning[0,i]))/2.d0 * (1.d0+Bg[iii,jjj,1]*COS(dpi/FSR[1]*lam+Phig[iii,jjj,1]+tuning[1,i]))/2.d0 * (1.d0+Bg[iii,jjj,2]*COS(dpi/FSR[2]*lam+Phig[iii,jjj,2]+tuning[2,i]))/2.d0)*I0g[iii,jjj]*dlam    
        ENDIF
    ENDFOR
ENDFOR

RESTORE,'./CPT/SEQUENCE_listSun070905_561378_256.BIN' ;CALMODE
a = WHERE(Inten EQ 0.d0)
IF (a[0] NE -1) THEN Inten[a] = 0.0000000000001d0

!P.MULTI=[0,2,3]
loadct,4
DEVICE,file='yo2.ps',xoffset=-0.8,yoffset=0,xsize=22,ysize=27,/color,bits=24
FOR i=0,nseq-1 DO TVIM,(Inten0[*,*,i]-Inten[*,*,i])/Inten[*,*,i],range=[-0.03,0.03],/scale,tit=STRING(i+1),stit='Relative error',pcharsize=1.5
;plot,lam,profileg-1.d0,xrange=[-2.d0*FSR[0],2.d0*FSR[0]],xst=1,charsize=1.5
FOR i=0,nseq-1 DO BEGIN
    temp = REFORM((Inten0[*,*,i]-Inten[*,*,i])/Inten[*,*,i])
    aa   = WHERE(temp LT 0.04 AND temp GT -0.04,na,COMPLEMENT=b)
    hist = histogram(temp[aa],binsize=0.001,min=-0.03,max=0.03)
    plot,FINDGEN(N_ELEMENTS(hist))*0.001-0.03,hist/FLOAT(na),psym=10,tit='!7r!17='+STRING(SIGMA(temp[aa]))+'/!7l!17='+STRING(MEAN(temp[aa])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5
    PRINT,MIN(temp[aa]),MAX(temp[aa]),SIGMA(temp[aa]),MEAN(temp[aa])
ENDFOR
DEVICE,/close

SET_PLOT,'x'

READ,pause

END
