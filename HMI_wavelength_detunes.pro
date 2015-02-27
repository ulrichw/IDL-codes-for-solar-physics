; PROGRAM TO FIT FOR THE LASER WAVELENGTH OF THE DETUNES TAKEN DURING
; THE WAVELENGTH DEPENDENCE TEST OF SEPTEMBER 2007, VACUUM CALIBRATION


;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
; NB: be aware of the fact that there might be a systematic error on the
; determination of the solar Doppler velocity due to the lack of
; absolute reference in the wavelength scale (we assume that we are
; centered on 6173.3433 A)
;
;----------------------------------------------------------------------

PRO HMI_wavelength_detunes

lamref  = 6173.3433d0;reference central wavelength for the FeI 6173 line
  

time0    = SYSTIME(1)
dpi      = 2.d0*!dpi

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------
nseq     = 27  ; number of positions in the sequence
nl0      = 33;30  ; number of wavelengths

; MEASURED INTENSITIES
;--------------------------------------------
list= LONARR(nl0)
lam0= DBLARR(nl0)

;GOTO,cal
; OBSMODE
list[0] = 710420
list[1] = 710600
list[2] = 710720
list[3] = 710900
list[4] = 711020
list[5] = 711140
list[6] = 711260
list[7] = 711380
list[8] = 711500
list[9] = 711620
list[10]= 711740
list[11]= 711860
list[12]= 711980
list[13]= 712100
list[14]= 712220
list[15]= 712408
list[16]= 712528
list[17]= 712708
list[18]= 712828
list[19]= 712948
list[20]= 713068
list[21]= 713188
list[22]= 713308
list[23]= 713428
list[24]= 713608
list[25]= 713728
list[26]= 713848
list[27]= 713968
list[28]= 714088
list[29]= 714208
list[30]= 710540
list[31]= 710840
list[32]= 712648

; GUESS WAVELENGTHS IN AIR
lam0[0] = 6173.3090d0
lam0[1] = 6173.4590d0
lam0[2] = 6173.6991d0
lam0[3] = 6173.8387d0
lam0[4] = 6174.0989d0
lam0[5] = 6174.2985d0
lam0[6] = 6174.6388d0
lam0[7] = 6174.8985d0
lam0[8] = 6175.8084d0
lam0[9] = 6176.0085d0
lam0[10]= 6176.3082d0
lam0[11]= 6176.5084d0
lam0[12]= 6178.6074d0
lam0[13]= 6179.0077d0
lam0[14]= 6175.4286d0
lam0[15]= 6173.2988d0
lam0[16]= 6173.1089d0
lam0[17]= 6172.9093d0
lam0[18]= 6172.7394d0
lam0[19]= 6172.5095d0
lam0[20]= 6172.2991d0
lam0[21]= 6172.0096d0
lam0[22]= 6171.7094d0
lam0[23]= 6170.8097d0
lam0[24]= 6170.6096d0
lam0[25]= 6170.3099d0
lam0[26]= 6168.0107d0
lam0[27]= 6167.6104d0
lam0[28]= 6171.3096d0
lam0[29]= 6173.6088d0
lam0[30]= 6173.4590d0
lam0[31]= 6173.8392d0
lam0[32]= 6172.9095d0

cal:
GOTO,rogert
; CALMODE
list[0 ] = 710480
list[1 ] = 710660
list[2 ] = 710780
list[3 ] = 710960
list[4 ] = 711080
list[5 ] = 711200
list[6 ] = 711320
list[7 ] = 711440
list[8 ] = 711560
list[9 ] = 711680
list[10] = 711800
list[11] = 711920
list[12] = 712040
list[13] = 712160
list[14] = 712280
list[15] = 712468
list[16] = 712588
list[17] = 712768
list[18] = 712888
list[19] = 713008
list[20] = 713128
list[21] = 713248
list[22] = 713368
list[23] = 713488
list[24] = 713668
list[25] = 713788
list[26] = 713908
list[27] = 714028
list[28] = 714148
list[29] = 714268

; GUESS WAVELENGTHS IN AIR
lam0[0] = 6173.3290d0
lam0[1] = 6173.4590d0
lam0[2] = 6173.6991d0
lam0[3] = 6173.8387d0
lam0[4] = 6174.0989d0
lam0[5] = 6174.2985d0
lam0[6] = 6174.6388d0
lam0[7] = 6174.8985d0
lam0[8] = 6175.8084d0
lam0[9] = 6176.0085d0
lam0[10]= 6176.3082d0
lam0[11]= 6176.5084d0
lam0[12]= 6178.6074d0
lam0[13]= 6179.0077d0
lam0[14]= 6175.4286d0
lam0[15]= 6173.2988d0
lam0[16]= 6173.1089d0
lam0[17]= 6172.9093d0
lam0[18]= 6172.7394d0
lam0[19]= 6172.5095d0
lam0[20]= 6172.2991d0
lam0[21]= 6172.0096d0
lam0[22]= 6171.7094d0
lam0[23]= 6170.8097d0
lam0[24]= 6170.6096d0
lam0[25]= 6170.3099d0
lam0[26]= 6168.0107d0
lam0[27]= 6167.6104d0
lam0[28]= 6171.3096d0
lam0[29]= 6173.6088d0
rogert:



lam0    = lam0-lamref
lam0g   = lam0 
I0g     = FLTARR(nl0)

nx      = 256

anglim  = 975               ; for Obsmode
xcenter = nx/2
ycenter = nx/2
distance= SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx
a       = WHERE(distance le anglim,COMPLEMENT=b)
distance[a]=1.d0
distance[b]=0.d0
 
SET_PLOT,'PS'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color
LOADCT,3
    
FOR iii=0,nl0-1 DO BEGIN

    FSN     = LONARR(nseq)
    FSN[0:nseq-1] = list[iii]+FINDGEN(nseq)*2+4 ; for front camera

    RESTORE,'CPT/SEQUENCE_listLaser071015_'+STRTRIM(list[iii],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
    Inten   = imf[*,*,2:28] ; front camera
                 
   ;RESTORE,'CORRECTION_INTENSITY_710420.BIN'
    RESTORE,'CORRECTION_INT_710420.BIN'
    INTENSITYdetunesf=INTENSITYdetunesf-0.03156d-9 ; front camera, we remove the background intensity

    FOR j=0,nseq-1 DO BEGIN
        a=WHERE(FSNdetunesf EQ FSN[j])
        Inten[*,*,j]=Inten[*,*,j]*distance/INTENSITYdetunesf[a[0]]*MEAN(INTENSITYdetunesf)
    ENDFOR

    Inten2 = FLTARR(nseq)
    FOR i=0,nseq-1 DO BEGIN
        temp    = REFORM(Inten[*,*,i])
        a       = WHERE(temp NE 0.0)
        Inten2[i] = TOTAL(temp[a])
    ENDFOR
    Inten = Inten2
   ;Inten = Inten-MIN(Inten)  ;WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!


; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------

    Bg          = DBLARR(3)     ; contrasts of the tunable elements
    Phig        = DBLARR(3)     ; relative phases
    thresh      = TOTAL(Inten)*0.8d0
    nparam      = 2             ; number of parameters we fit for
    Residual    = DBLARR(nseq)  ; residuals of the least-squares fit
    dIntendlam0 = DBLARR(nseq)  ; derivatives to compute the Jacobian matrix
    dIntendIc   = DBLARR(nseq)
    tuning      = DBLARR(3,nseq) ; tuning sequence
    maxsteps    = 150            ; maximum steps allowed for the iterative Least-Squares algorithm
    history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err         = DBLARR(maxsteps)  ; error estimate at each step    


; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

; VALUES OF JUNE 2009
FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.17094240d0;0.1710098d0            ; for the narrow-band Michelson
FSR[1]      = 0.34192317d0;0.3421506d0            ; for the broad-band  Michelson
FSR[2]      = 0.69348347d0;0.6943613d0            ; for E1
FSR[3]      = 1.407d0;*0.993d0                    ; for E2
FSR[4]      = 2.779d0;*0.993d0                    ; for E3
FSR[5]      = 5.682d0;*0.993d0                    ; for E4
FSR[6]      = 11.354d0;*0.993d0                   ; for E5


; DEFINITION OF THE DE-TUNE SEQUENCE
;----------------------------------------------------------

tuningJS       = DBLARR(nseq,3)
tuningJS[0 ,*] = [         0.d0,          0.d0 ,         0.d0]
tuningJS[1 ,*] = [        80.d0,          0.d0 ,         0.d0]
tuningJS[2 ,*] = [       160.d0,          0.d0 ,         0.d0]
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

    
; RELATIVE PHASES 
;------------------------------------------------

   ;Phig[0:2] = [-134.244,4.3051,-137.093] ; with 541606 ; CALMODE    
   ;Bg[0:2]   = [0.97839180,0.98647068,0.97871106]

   ;Phig[0:2] = [-133.66125,3.2847898,-138.94738] ; with 549732 ; OBSMODE    
   ;Bg[0:2]   = [0.96895287,0.98816689,0.96285676]
    
   ;Phig[0:2] = [-114.559d0,10.371d0,-138.415d0]   ; 707090 CALMODE
   ;Bg[0:2]   = [0.96895287,0.98816689,0.96285676] ; FROM CPT/phases_contrast.txt

   ;Phig[0:2] = [-116.0406,4.6753,-138.4487] ; 712408 OBSMODE
   ;Bg[0:2]   = [0.96895287,0.98816689,0.96285676]

   ;Phig[0:2] = [10.390520,-99.526381,-137.42685] ; 946934 OBSMODE
   ;Bg[0:2]   = [0.9695,0.9905,0.9764]

   ;Phig[0:2] = [54.083571,-61.131336,-137.40789] ; 566521 CALMODE
   ;Bg[0:2]   = [0.98203985,0.98528701,0.98451732]; for wavelength 6173.573 A

   ;Phig[0:2] = [-118.33074d0,3.6014040d0,-139.21310d0] ; FROM CPT/phases_contrast.txt
   ;Bg[0:2]   = [0.96803371d0,0.99558504d0,0.96608620d0]

    Phig[0:2] = [-113.287d0,5.26063d0,-138.755d0] ; FROM CPT/phase_contrast_june2009.txt (MEDIAN VALUES, OBSMODE, intensity correction)
    Bg[0:2]   = [0.9686d0,0.9939d0,0.9657d0]

   ;Phig[0:2] = [-112.853d0,6.20089d0,-135.252d0] ; FROM CPT/phase_contrast_june2009.txt (MEDIAN VALUES, CALMODE, intensity correction)
   ;Bg[0:2]   = [0.96928d0,0.991368d0,0.976951d0]



   ;Phig[0:2] = [-114.852d0,10.186d0,-139.284d0]  ; FROM HMI_SUN1.PRO ON CALMODE DETUNE FSN=707090; front camera SO CONTAMINATED BY CONVECTIVE BLUESHIFT AND LIMB EFFECT
   ;Bg[0:2]   = [0.96928d0,0.991368d0,0.976951d0]


    Phig[*]   = Phig[*] / 180.d0 * !dpi

; LEAST-SQUARES ALGORITHM
;-------------------------------------------------

;    WINDOW,0,RETAIN=2

    converg  = 0
    jj       = 0
    I0g[iii] = thresh
    
    WHILE (converg EQ 0) DO BEGIN
        

        FOR j=0,nseq-1 DO BEGIN
            
            profileg    = 1.d0
            FOR i=0,2 DO profileg  = profileg  * (1.d0+Bg[i]*COS(dpi/FSR[i]*lam0g[iii]+Phig[i]+tuningJS[j,i]))/2.d0
            
            dh = 0.0001d0
            profilegp2h = 1.d0
            profilegph  = 1.d0
            profilegm2h = 1.d0
            profilegmh  = 1.d0
            FOR i=0,2 DO profilegp2h = profilegp2h * (1.d0+Bg[i]*COS(dpi/FSR[i]*(lam0g[iii]+2.d0*dh)+Phig[i]+tuningJS[j,i]))/2.d0
            FOR i=0,2 DO profilegph  = profilegph  * (1.d0+Bg[i]*COS(dpi/FSR[i]*(lam0g[iii]+dh)+Phig[i]+tuningJS[j,i]))/2.d0
            FOR i=0,2 DO profilegm2h = profilegm2h * (1.d0+Bg[i]*COS(dpi/FSR[i]*(lam0g[iii]-2.d0*dh)+Phig[i]+tuningJS[j,i]))/2.d0
            FOR i=0,2 DO profilegmh  = profilegmh  * (1.d0+Bg[i]*COS(dpi/FSR[i]*(lam0g[iii]-dh)+Phig[i]+tuningJS[j,i]))/2.d0



            Residual[j]   =  Inten[j] - profileg*I0g[iii]
            dIntendIc[j]  =  profileg
            dIntendlam0[j]=  I0g[iii]*(-profilegp2h+8.d0*profilegph-8.d0*profilegmh+profilegm2h)/12.d0/dh ; 5-point numerical derivative



        ENDFOR

      ; Jacobian matrix
        Jac  = DBLARR(nparam,nseq)
        FOR i= 0,nseq-1 DO BEGIN
            
            Jac[0,i]  = dIntendIc[i]
            Jac[1,i]  = dIntendlam0[i]
            
        ENDFOR
        
        LA_SVD,Jac,W,U,V,/DOUBLE
        
        
        Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
        
        temp   = TRANSPOSE(Dx)/[I0g[iii],lam0g[iii]]
        err[jj]= MAX(ABS(temp))
        
        IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7)     THEN converg = 2 
        IF(err[jj] LE 1.d-7)                          THEN converg = 1
        
        IF(converg EQ 2) THEN BEGIN
            j     = WHERE(err EQ MIN(err))
            I0g[iii]   = history[j[0],0]
            lam0g[iii] = history[j[0],1]
        ENDIF
        
        IF(converg NE 2) THEN BEGIN
            I0g[iii] = I0g[iii] + Dx[0]
            IF(I0g[iii] LT 0.002*thresh OR I0g[iii] GT 100.0*thresh) THEN I0g[iii] = thresh
            lam0g[iii] = lam0g[iii] + Dx[1]
            IF(ABS(lam0g[iii]) GT  10.d0) THEN lam0g[iii]  = lam0[iii] ;-0.030
        ENDIF
        
        history[jj,*] = [I0g[iii],lam0g[iii]]
        
        jj = jj+1
        
    ENDWHILE

  ; RECONSTRUCTION
    Intenrecons=FLTARR(nseq)
    FOR j=0,nseq-1 DO BEGIN
        profileg=1.0
        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[i]*COS(dpi/FSR[i]*lam0g[iii]+Phig[i]+tuningJS[j,i]))/2.d0
        Intenrecons[j]  = profileg*I0g[iii]
    ENDFOR
   ;PRINT,list[iii],lam0g[iii]+lamref,lam0[iii]+lamref,SQRT(TOTAL((Intenrecons-Inten)^2.))/SQRT(TOTAL(Inten^2.0)),MEAN(Inten),I0g[iii],converg
    PRINT,list[iii],lam0g[iii]+lamref,lam0[iii]+lamref,SQRT(TOTAL((Intenrecons-Inten)^2.))/1000.d0,MEAN(Inten)/1000.d0,I0g[iii],converg
    
    PLOT,Inten,xst=1
    OPLOT,Intenrecons,linestyle=2,thick=2,col=180

   ;READ,pause

ENDFOR

DEVICE,/CLOSE
SET_PLOT,'X'

END

