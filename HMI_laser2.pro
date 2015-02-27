; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE CO-TUNE SEQUENCES
; OBJECTIVE: TO DERIVE THE NON-TUNABLE TRANSMISSION PROFILE
;----------------------------------------------------------------------

FUNCTION reconstruct,moy,moyb

; moy = average phases of non-tunable elements in degrees
; moyb= average contrasts

lamref      = 6173.3433d0
dpi         = 2.d0*!dpi
FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172d0;0.17076972d0         ; for the narrow-band Michelson
FSR[1]      = 0.344d0;0.34158512d0         ; for the broad-band  Michelson
FSR[2]      = 0.693d0            ; for E1
FSR[3]      = 1.407d0            ; for E2
FSR[4]      = 2.779d0            ; for E3
FSR[5]      = 5.682d0            ; for E4
FSR[6]      = 11.354d0           ; for E5

; RECONSTRUCT THE SPATIALLY-AVERAGED TRANSMISSION PROFILE OF THE NON TUNABLE ELEMENTS:
lam         = FINDGEN(20000)/19999.*30.-15.
RESTORE,'frontwindow.bin'
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q           = READFITS('blocker11.fits')
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lamref,lam) ; I center the profile ;+3.61328

profileg    = blocker0
FOR i=3,6 DO profileg = profileg * (1.d0+moyb[i-3]*COS(dpi/FSR[i]*lam+moy[i-3]*!pi/180.))/2.d0
;SET_PLOT,'x'
;WINDOW,0,RETAIN=2
!P.MULTI=0
SET_PLOT,'PS'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=16,/color
LOADCT,3
a=where(lam ge -2 and lam le 2)
PLOT,lam,profileg/max(profileg[a]),xst=1,xrange=[-5.25,4.],yrange=[0,1],yst=1,charsize=1.5,tit='!17',xtit='Wavelength (A)',ytit='Transmittance'
oplot,lam,profileg/max(profileg[a]),col=180
oplot,[-7.5,7.5],[0.05,0.05],linestyle=2
;oplot,[-0.5,-0.5],[0,1],linestyle=2
;oplot,[0.5,0.5],[0,1],linestyle=2
oplot,[0,0],[0,1],linestyle=2
oplot,[-7.5,7.5],[0.5,0.5],linestyle=2
DEVICE,/CLOSE
SET_PLOT,'x'
res  = FLTARR(2,20000)
res[0,*] = lam
res[1,*] = profileg/max(profileg[a])

RETURN,res

END


;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser2,draw

;SET_PLOT,'X'
;WINDOW,0,RETAIN=2,xsize=900,ysize=900
;!P.MULTI    = [0,2,2]
time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0        ; Fe I 6173 line central wavelength in air

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 64               ; number of rows
ny          = 64               ; number of columns
factor      = 100.
nseq        = 27               ; number of positions in the co-tune sequence
nl0         = 16.              ; number of wavelength positions for the dye laser
nparam      = 9                ; number of filter parameters we fit for with the cotune: 4 phases, 4 contrasts, laser intensity
anglim      = 950
xcenter     = nx/2
ycenter     = nx/2
distance    = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;------------------------------------------------------

FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172d0;0.17076972d0         ; for the narrow-band Michelson
FSR[1]      = 0.344d0;0.34158512d0         ; for the broad-band  Michelson
FSR[2]      = 0.693d0            ; for E1
FSR[3]      = 1.407d0            ; for E2
FSR[4]      = 2.779d0            ; for E3
FSR[5]      = 5.682d0            ; for E4
FSR[6]      = 11.354d0           ; for E5


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;------------------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin'
blocker0     = transmission/100.d0;INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
lam          = wavelength*10.d0-lamref
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lamref,lam) ; I center the profile ;+3.61328

; MEASURED INTENSITIES
;-------------------------------------------------------

; list containing the names of filtergrams
; USE readfile3_128 TO PRODUCE THE INTENSITIES


list    = STRARR(nl0)
expo    = FLTARR(nl0)
lam0    = FLTARR(nl0)
relI    = FLTARR(nl0)

list[0 ]= 'listLaser070905_549188'
list[1 ]= 'listLaser070905_549324'
list[2 ]= 'listLaser070905_549596'
list[3 ]= 'listLaser070905_549732'
list[4 ]= 'listLaser070905_549868'
list[5 ]= 'listLaser070905_550004'
list[6 ]= 'listLaser070905_550412'
list[7 ]= 'listLaser070905_550548'
list[8 ]= 'listLaser070905_550820'
list[9 ]= 'listLaser070905_550956'
list[10]= 'listLaser070905_551092'
list[11]= 'listLaser070905_551364'
list[12]= 'listLaser070905_551500'
list[13]= 'listLaser070905_551636'
list[14]= 'listLaser070905_551772'
list[15]= 'listLaser070905_551908'
lam0[0 ] = 6173.0491;6172.3566
lam0[1 ] = 6172.9141;6172.2250
lam0[2 ] = 6172.4215;6171.7294
lam0[3 ] = 6172.1727;6171.4854
lam0[4 ] = 6171.6754;6170.9934
lam0[5 ] = 6171.1905;6170.5000
lam0[6 ] = 6169.6988;6169.0295
lam0[7 ] = 6173.4235;6172.7279
lam0[8 ] = 6173.6538;6172.9519
lam0[9 ] = 6173.7778;6173.0754
lam0[10] = 6173.9083;6173.2033
lam0[11] = 6174.4087;6173.6985
lam0[12] = 6174.6634;6173.9552
lam0[13] = 6174.8915;6174.1862
lam0[14] = 6175.1511;6174.4400
lam0[15] = 6175.6438;6174.9383
lam0     = lam0-lamref;+0.693d0  !!!WARNING!!!
relI[0 ] = 4.56688 ;4.60826
relI[1 ] = 2.95802; 2.96647
relI[2 ] = 3.43708; 3.43317
relI[3 ] = 4.29308; 4.31731
relI[4 ] = 3.13255; 3.12851
relI[5 ] = 2.49317; 2.49997
relI[6 ] = 4.04419; 4.05959
relI[7 ] = 3.97266; 3.98424
relI[8 ] = 3.1608; 3.15832
relI[9 ] = 3.26378; 3.27738
relI[10] = 2.8725; 2.86727
relI[11] = 2.6067; 2.57554
relI[12] = 2.96883; 2.95981
relI[13] = 3.15652; 3.13853
relI[14] = 2.93554; 2.93016
relI[15] = 2.81568; 2.80612
expo[0]  = 1005.
expo[1]  = 1005.
expo[2]  = 3755.
expo[3]  = 3755.
expo[4]  = 3755.
expo[5]  = 3755.
expo[6]  = 3755.
expo[7]  = 1005.
expo[8]  = 1005.
expo[9]  = 1005.
expo[10] = 3755.
expo[11] = 3755.
expo[12] = 3755.
expo[13] = 3755.
expo[14] = 3755.
expo[15] = 3755.


;list[0] = 'listLaser070905_549256'
;list[1] = 'listLaser070905_549392'
;list[2] = 'listLaser070905_549664'
;list[3] = 'listLaser070905_549800'
;list[4] = 'listLaser070905_549936'
;list[5] = 'listLaser070905_550072'
;list[6] = 'listLaser070905_550480'
;list[7] = 'listLaser070905_550616'
;list[8 ]= 'listLaser070905_550888'
;list[9 ]= 'listLaser070905_551024'
;list[10]= 'listLaser070905_551160'
;list[11]= 'listLaser070905_551432'
;list[12]= 'listLaser070905_551568'
;list[13]= 'listLaser070905_551704'
;list[14]= 'listLaser070905_551840'
;list[15]= 'listLaser070905_551976'
;lam0[0] = 6173.0489            ; 0.0245983
;lam0[1] = 6172.9138            ; 0.0301045
;lam0[2] = 6172.4212            ; 0.0246679
;lam0[3] = 6172.1727            ; 0.0227419
;lam0[4] = 6171.6749            ; 0.0671509
;lam0[5] = 6171.1903            ; 0.0570870
;lam0[6] = 6169.6986            ;  0.102417
;lam0[7] = 6173.4235            ; 0.0336450
;lam0[8] = 6173.6535            ; 0.0538418
;lam0[9] = 6173.7776            ; 0.0525744
;lam0[10]= 6173.9083            ; 0.0870561
;lam0[11]= 6174.4092            ; 0.0642272
;lam0[12]= 6174.6632            ; 0.0770916
;lam0[13]= 6174.8915            ; 0.0963618
;lam0[14]= 6175.1514            ;  0.115435
;lam0[15]= 6175.6437            ;  0.115777
;lam0     = lam0-lamref
;relI[*] = 0.d0
;expo[0]  = 1005.
;expo[1]  = 1005.
;expo[2]  = 3755.
;expo[3]  = 3755.
;expo[4]  = 3755.
;expo[5]  = 3755.
;expo[6]  = 3755.
;expo[7]  = 1005.
;expo[8]  = 1005.
;expo[9]  = 1005.
;expo[10] = 3755.
;expo[11] = 3755.
;expo[12] = 3755.
;expo[13] = 3755.
;expo[14] = 3755.
;expo[15] = 3755.


nx2     = 256
ny2     = 256
Inten   = DBLARR(nx2,ny2,nseq*nl0) ; measured output intensities (measured on a HMI CCD) with the cotune
FOR i=0,nl0-1 DO BEGIN
    RESTORE,'CPT/SEQUENCE_'+STRTRIM(list[i],1)+'_'+STRTRIM(STRING(LONG(nx2)),1)+'.BIN'
    Inten[*,*,nseq*i:nseq*(i+1)-1] = imx[*,*,2:28]/factor
ENDFOR

RESTORE,'CORRECTION_INTENSITY_549188.BIN'; OBTAINED BY intensity_correction.pro RECOMPUTE EACH TIME !!!!

; WE NORMALIZE THE INTENSITIES
FOR i=0,nl0-1 DO BEGIN
    FOR j=0,nseq-1 DO BEGIN
       ;Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/relI[i]*MEAN(relI)/expo[i]*MEAN(expo)
        Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/INTENSITYdetunes[i*nseq+j]*MEAN(INTENSITYdetunes)/expo[i]*MEAN(expo)
    ENDFOR
ENDFOR

; WE REQUIRE INTENSITIES TO BE LARGER THAN 0 !!! BIAS !!!
;temp = REFORM(REBIN(inten,1,1,nseq*nl0))
;IF(MIN(temp) LT 0.0) THEN inten= inten-MIN(temp)

prof=fltarr(nl0)
for i=0,nl0-1 do prof[i]=rebin(Inten[*,*,nseq*i:nseq*(i+1)-1],1,1,1)
prof=prof[sort(lam0)]
moncul=lam0[sort(lam0)]
set_plot,'x'
!p.multi=0
window,0,retain=2
plot,moncul+lamref,prof,xst=1
read,p

Inten   = REBIN(Inten,nx,ny,nseq*nl0)

; VARIABLE DEFINITION
;--------------------------------------------------------

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"
Residual    = DBLARR(nseq*nl0) ; residual for the least-squares fit
dIntendB3   = DBLARR(nseq*nl0) ; derivatives to compute the Jacobian matrix
dIntendB4   = DBLARR(nseq*nl0)
dIntendB5   = DBLARR(nseq*nl0)
dIntendB6   = DBLARR(nseq*nl0)
dIntendPhi3 = DBLARR(nseq*nl0)
dIntendPhi4 = DBLARR(nseq*nl0)
dIntendPhi5 = DBLARR(nseq*nl0)
dIntendPhi6 = DBLARR(nseq*nl0)
dIntendInt  = DBLARR(nseq*nl0)
lateconv    = INTARR(nx,ny)    ; to discard spurious results
maxsteps    = 100.;5555
Jac         = DBLARR(nparam,nseq*nl0) ; Jacobian matrix of the least-squares fit
contrasts   = 0.95d0


; DEFINITION OF THE CO-TUNE SEQUENCE
; co-tuning in the range [-2 FSR[E1],+2 FSR[E1]]
; TABLE PROVIDED BY JESPER
;----------------------------------------------------------

tuningJS       = DBLARR(nseq,3)
tuningJS[0 ,*] = [         0.d0,          0.d0,          0.d0]
tuningJS[1 ,*] = [        80.d0,          0.d0,          0.d0]
tuningJS[2 ,*] = [       160.d0,          0.d0,          0.d0]
tuningJS[3 ,*] = [         0.d0,         80.d0,          0.d0]
tuningJS[4 ,*] = [        80.d0,         80.d0,          0.d0]
tuningJS[5 ,*] = [       160.d0,         80.d0,          0.d0]
tuningJS[6 ,*] = [         0.d0,        160.d0,          0.d0]
tuningJS[7 ,*] = [        80.d0,        160.d0,          0.d0]
tuningJS[8 ,*] = [       160.d0,        160.d0,          0.d0]
tuningJS[9 ,*] = [         0.d0,          0.d0,         80.d0]
tuningJS[10,*] = [        80.d0,          0.d0,         80.d0]
tuningJS[11,*] = [       160.d0,          0.d0,         80.d0]
tuningJS[12,*] = [         0.d0,         80.d0,         80.d0] 
tuningJS[13,*] = [        80.d0,         80.d0,         80.d0]
tuningJS[14,*] = [       160.d0,         80.d0,         80.d0]
tuningJS[15,*] = [         0.d0,        160.d0,         80.d0]
tuningJS[16,*] = [        80.d0,        160.d0,         80.d0]
tuningJS[17,*] = [       160.d0,        160.d0,         80.d0]
tuningJS[18,*] = [         0.d0,          0.d0,        160.d0]
tuningJS[19,*] = [        80.d0,          0.d0,        160.d0]
tuningJS[20,*] = [       160.d0,          0.d0,        160.d0]
tuningJS[21,*] = [         0.d0,         80.d0,        160.d0]
tuningJS[22,*] = [        80.d0,         80.d0,        160.d0]
tuningJS[23,*] = [       160.d0,         80.d0,        160.d0]
tuningJS[24,*] = [         0.d0,        160.d0,        160.d0]
tuningJS[25,*] = [        80.d0,        160.d0,        160.d0]
tuningJS[26,*] = [       160.d0,        160.d0,        160.d0]

FOR i=0,nseq-1 DO tuningJS[i,*]  = tuningJS[i,*] * [6.d0,6.d0,-6.d0]*!dpi/180.d0


IF draw EQ 1 THEN GOTO,draw


; GUESS VALUES
;----------------------------------------------------

RESTORE,'CPT/CPT_laser_front_549732.BIN' ; OBSMODE
;RESTORE,'CPT/CPT_laser_front_549800.BIN' ; CALMODE


a=WHERE(FINITE(Bg0) EQ 0)
IF (a[0] NE -1) THEN Bg0[a]  =1.d0
a=WHERE(FINITE(Phig0) EQ 0)
IF (a[0] NE -1) THEN Phig0[a]=0.d0

;FOR iii=0,nx-1 DO FOR jjj=0,nx-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]+[-106.26-169.22735,42.35+4.2600189,-140.32+164.03892]*!pi/180.d0 ; 79336
;FOR iii=0,nx-1 DO FOR jjj=0,nx-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]+[-106.26-170.95412,42.35+2.8555888,-140.32+158.35755]*!pi/180.d0 ; 79916

Bg0 = REBIN(Bg0,nx,nx,3)
Phig0=REBIN(Phig0,nx,nx,3)

Bg[*,*,0:2]   = Bg0[*,*,*]         ; contrasts and phases of the tuning elements
Phig[*,*,0:2] = Phig0[*,*,*]


 RESTORE,'I0g_549188.BIN' ; WARNING WARNING WARNING WARNING WARNING !!!!!!!!!!!!! for I0g in Obsmode
;RESTORE,'temp2.bin' ; WARNING WARNING WARNING WARNING WARNING !!!!!!!!!!!!! for I0g in Calmode



FOR jjj=0,ny-1 DO BEGIN

   ;TVIM,Phig[*,*,3],/scale 
   ;TVIM,Phig[*,*,4],/scale
   ;TVIM,Phig[*,*,5],/scale
   ;TVIM,Phig[*,*,6],/scale
    TVIM,lateconv

    FOR iii=0,nx-1 DO BEGIN

   ; BEGINNING OF THE LEAST-SQUARES FIT
   ;-----------------------------------


        IF(distance[iii,jjj] LE anglim) THEN BEGIN

           ;Bg[iii,jjj,3:6]  = [0.9,0.99,1.3,0.6] ; contrasts
           ;Phig[iii,jjj,3:6]= [-6.7,0.08,-1.45,177.]/180.d0*!dpi 
            Phig[iii,jjj,3:6]= [0.0,0.0,0.0,0.]/180.d0*!dpi
            Bg[iii,jjj,3:6]  = [1.0,1.,1.,1.]

            ; guess of the initial intensity
           ;sum = 0.d0
           ;FOR j=0,nseq-1 DO BEGIN
           ;    FOR ii=0,nl0-1 DO BEGIN
           ;        profileg  = INTERPOL(blocker0,lam,lam0[ii])
           ;        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
           ;        FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0
           ;        sum = sum+profileg
           ;    ENDFOR
           ;ENDFOR
           ;thresh       = TOTAL(Inten[iii,jjj,*])/sum
           ;I0g[iii,jjj] = thresh


            ;sum = FLTARR(nl0)
            ;FOR ii=0,nl0-1 DO BEGIN
            ;ii=10
            ;    profileg  = INTERPOL(blocker0,lam,lam0[ii])
            ;    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0
            ;    sum = 8.d0*TOTAL(Inten[iii,jjj,*])/27.d0/profileg
           ;ENDFOR
            ;thresh       = 900.;sum/nl0
            thresh=I0g[iii,jjj];1050.;4000.
            I0g[iii,jjj] = thresh
            ;PRINT,'yo',SIGMA(sum)


            converg = 0
            jj      = 0

            WHILE converg EQ 0 DO BEGIN

                FOR j=0,nseq-1 DO BEGIN
                    FOR ii=0,nl0-1 DO BEGIN

                        profileg  = INTERPOL(blocker0,lam,lam0[ii],/QUADRATIC)
                        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                        profilegf = profileg
                        FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0

                        Residual[ii*nseq+j]  = (Inten[iii,jjj,ii*nseq+j] - profileg*I0g[iii,jjj])
            
                        dIntendB3[ii*nseq+j] = profilegf*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3])*(1.d0+Bg[iii,jjj,4]$
                       *COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendB4[ii*nseq+j] = profilegf*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4])*(1.d0+Bg[iii,jjj,3]$
                       *COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))*(1.d0+Bg[iii,jjj,5]*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendB5[ii*nseq+j] = profilegf*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5])*(1.d0+Bg[iii,jjj,4]$
                       *COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))$
                       *(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendB6[ii*nseq+j] = profilegf*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6])*(1.d0+Bg[iii,jjj,4]$
                       *COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))/16.d0*I0g[iii,jjj]

                        dIntendPhi3[ii*nseq+j]= profilegf*(-Bg[iii,jjj,3]*SIN(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))$
                       *(1.d0+Bg[iii,jjj,4]*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]$
                       *COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))*(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendPhi4[ii*nseq+j]= profilegf*(-Bg[iii,jjj,4]*SIN(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))$
                       *(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))*(1.d0+Bg[iii,jjj,5]$
                       *COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))*(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendPhi5[ii*nseq+j]= profilegf*(-Bg[iii,jjj,5]*SIN(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,4]*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,3]$
                       *COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))*(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                       dIntendPhi6[ii*nseq+j]= profilegf*(-Bg[iii,jjj,6]*SIN(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))$
                       *(1.d0+Bg[iii,jjj,4]*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]$
                       *COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))*(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))/16.d0*I0g[iii,jjj]

                        dIntendInt[ii*nseq+j] = profileg

                    ENDFOR
                ENDFOR



              ; Computation of the Jacobian matrix
                FOR i= 0,nseq*nl0-1 DO BEGIN
                    Jac[0,i]  = dIntendB3[i]
                    Jac[1,i]  = dIntendB4[i]
                    Jac[2,i]  = dIntendB5[i]
                    Jac[3,i]  = dIntendB6[i]
                    Jac[4,i]  = dIntendPhi3[i]
                    Jac[5,i]  = dIntendPhi4[i]
                    Jac[6,i]  = dIntendPhi5[i]
                    Jac[7,i]  = dIntendPhi6[i]
                    Jac[8,i]  = dIntendInt[i]
                ENDFOR

                IF (iii EQ 32 AND JJJ EQ 10) THEN READ,pause 

                LA_SVD,Jac,W,U,V,/DOUBLE ; Singular Value Decomposition using the LAPACK algorithm

                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual

                err = TRANSPOSE(Dx)/[REFORM(Bg[iii,jjj,3]),REFORM(Bg[iii,jjj,4]),REFORM(Bg[iii,jjj,5]),REFORM(Bg[iii,jjj,6])$
,REFORM(Phig[iii,jjj,3]),REFORM(Phig[iii,jjj,4]),REFORM(Phig[iii,jjj,5]),REFORM(Phig[iii,jjj,6]),I0g[iii,jjj]]

                err = MAX(ABS(err))

                IF(jj EQ maxsteps-1 AND err GT 1.d-3 )      THEN          converg         = 2 ;d-7
                IF(err LE 1.d-3)                            THEN          converg         = 1

                Bg[iii,jjj,3]   = Bg[iii,jjj,3]  +Dx[0]
                IF(Bg[iii,jjj,3] LT 0.d0 OR Bg[iii,jjj,3] GT 10.4d0) THEN Bg[iii,jjj,3]   = contrasts ;10.4
                Bg[iii,jjj,4]   = Bg[iii,jjj,4]  +Dx[1]
                IF(Bg[iii,jjj,4] LT 0.d0 OR Bg[iii,jjj,4] GT 10.4d0) THEN Bg[iii,jjj,4]   = contrasts
                Bg[iii,jjj,5]   = Bg[iii,jjj,5]  +Dx[2]
                IF(Bg[iii,jjj,5] LT 0.d0 OR Bg[iii,jjj,5] GT 10.4d0) THEN Bg[iii,jjj,5]   = contrasts
                Bg[iii,jjj,6]   = Bg[iii,jjj,6]  +Dx[3]
                IF(Bg[iii,jjj,6] LT 0.d0 OR Bg[iii,jjj,6] GT 10.4d0) THEN Bg[iii,jjj,6]   = contrasts
                Phig[iii,jjj,3] = Phig[iii,jjj,3]+Dx[4]
                IF(ABS(Phig[iii,jjj,3]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,3] = 0.d0
                Phig[iii,jjj,4] = Phig[iii,jjj,4]+Dx[5]
                IF(ABS(Phig[iii,jjj,4]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,4] = 0.d0
                Phig[iii,jjj,5] = Phig[iii,jjj,5]+Dx[6]
                IF(ABS(Phig[iii,jjj,5]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,5] = 0.d0
                Phig[iii,jjj,6] = Phig[iii,jjj,6]+Dx[7]
                IF(ABS(Phig[iii,jjj,6]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,6] = 0.d0
                I0g[iii,jjj]    = I0g[iii,jjj]   +Dx[8]
                IF(I0g[iii,jjj] LT 0.2*thresh OR I0g[iii,jjj] GT 4.0*thresh) THEN I0g[iii,jjj] = thresh

                lateconv[iii,jjj] = converg
                jj = jj+1

            ENDWHILE

            PRINT,' '
            PRINT,'PIXEL',iii,jjj
            PRINT,'PARAMETERS='
            PRINT,'CONTRASTS'
            FOR i=3,6 DO PRINT,'Bgi =',Bg[iii,jjj,i]
            PRINT,'PHASES'
            Phig[iii,jjj,*] = Phig[iii,jjj,*]*180.d0/!dpi  ; in degrees
            FOR i=3,6 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'LASER INTENSITY'
            PRINT,'I0g=',I0g[iii,jjj]
            PRINT,'lateconv=',converg            

        ENDIF

    ENDFOR
ENDFOR

SAVE,Bg,Phig,I0g,lateconv,FILE='RESULTS/RESULTS2_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'ELAPSED TIME=',SYSTIME(1)-time0

draw:

;-----------------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS2_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a=WHERE(distance GT anglim OR lateconv NE 1,COMPLEMENT=b) ; depends on the size of the target ;963 obsmode; 925 calmode

FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & temp[a]=-10000.d0 & phig[*,*,i]=temp[*,*] & temp=REFORM(Bg[*,*,i]) & temp[a]=-10000.d0 & Bg[*,*,i]=temp[*,*] & phig[0:1,*,i]=-10000.d0 & Bg[0:1,*,i]=-10000.d0 & ENDFOR

a = WHERE(FINITE(Phig) EQ 0)
IF(a[0] NE -1) THEN Phig[a]= -10000.d0
a = WHERE(FINITE(Bg) EQ 0)
IF(a[0] NE -1) THEN Bg[a]  = -10000.d0


; ADD 180 DEGREES TO PHASES < -180
FOR i=3,6 DO BEGIN
temp = REFORM(Phig[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        Phig[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        Phig[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(4)
moyb=moy
FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moy[i-3]=MEAN(temp[a]) & temp[b]=moy[i-3] & Phig[*,*,i]=temp & temp=REFORM(Bg[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyb[i-3]=MEAN(temp[a]) & temp[b]=moyb[i-3] & Bg[*,*,i]=temp  & ENDFOR

; WE PLOT THE RESULT
SET_PLOT,'ps'
!P.MULTI=[0,2,2]
device,file='yo.ps',bits=24,xoffset=-0.5,yoffset=1,xsize=22.5,ysize=19,/color
LOADCT,4
tvim,phig[*,*,3],/scale,tit='!17E2 '+STRTRIM(STRING(moy[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,3]),MAX(phig[*,*,3])]
tvim,phig[*,*,4],/scale,tit='!17E3 '+STRTRIM(STRING(moy[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,4]),MAX(phig[*,*,4])]
tvim,phig[*,*,5],/scale,tit='!17E4 '+STRTRIM(STRING(moy[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,5]),MAX(phig[*,*,5])]
tvim,phig[*,*,6],/scale,tit='!17E5 '+STRTRIM(STRING(moy[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,6]),MAX(phig[*,*,6])]
tvim,Bg[*,*,3]  ,/scale,tit='!17E2',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,3]),1.02]
tvim,Bg[*,*,4]  ,/scale,tit='!17E3',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,4]),1.04]
tvim,Bg[*,*,5]  ,/scale,tit='!17E4',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[1.1,1.4]
tvim,Bg[*,*,6]  ,/scale,tit='!17E5',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,6]),1.02]
tvim,I0g,/scale,tit='!17Intensity',barwidth=0.5,range=[MIN(I0g),MAX(I0g)]
DEVICE,/close
PRINT,'PHASE AVERAGES'
FOR i=0,3 DO PRINT,moy[i]
PRINT,'CONTRAST AVERAGES'
FOR i=0,3 DO PRINT,moyb[i]

; RECONSTRUCTION OF THE SEQUENCE

Intenrecons = FLTARR(nx,ny,nseq*nl0)
a=WHERE(distance GT anglim,COMPLEMENT=b)
distance[a] = 0.0
distance[b] = 1.0
FOR i=0,nseq*nl0-1 DO Inten[*,*,i] = Inten[*,*,i]*distance

FOR jjj=0,ny-1 DO BEGIN
    PRINT,jjj
    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] NE 0.0) THEN BEGIN
            
            FOR j=0,nseq-1 DO BEGIN
                FOR ii=0,nl0-1 DO BEGIN
                    
                    profileg  = INTERPOL(blocker0,lam,lam0[ii])
                    FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!pi/180.+tuningJS[j,i]))/2.d0
                    profilegf = profileg
                    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!pi/180.))/2.d0
                    
                    Intenrecons[iii,jjj,ii*nseq+j]  = profileg*I0g[iii,jjj]
                    
                ENDFOR
            ENDFOR
            
        ENDIF
        
    ENDFOR
ENDFOR

!P.MULTI=0
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='yo3.ps',xsize=35,ysize=25,xoffset=0,yoffset=0,/color
LOADCT,3
PLOT,REBIN(Inten,1,1,nseq*nl0),xst=1,tit='!17',xtit='position number',ytit='intensity',charsize=1.5,yst=1,thick=3;,/ylog
;OPLOT,REBIN(Inten,1,1,nseq*nl0),psym=2
OPLOT,REBIN(Intenrecons,1,1,nseq*nl0),color=150,linestyle=2,thick=2
device,/close
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten,1,1,nseq*nl0))-REFORM(REBIN(Intenrecons,1,1,nseq*nl0)))^2.0 )


a           = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0
DEVICE,file='yo2.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,bits=24,/color
LOADCT,3
!P.multi=[0,2,5]
FOR i=0,nseq*nl0-1 DO BEGIN & temp=REFORM((inten[*,*,i]-intenrecons[*,*,i])/inten[*,*,i])*distance & a=WHERE(FINITE(temp) EQ 0) & IF(a[0] NE -1) THEN temp[a]=0.0 & TVIM,temp,/scale,barwidth=0.5,range=[-0.1,0.1] & ENDFOR
DEVICE,/CLOSE




READ,pause

END

