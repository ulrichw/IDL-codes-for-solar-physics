; PROGRAM TO FIT FOR THE FSRs OF THE TUNABLE ELEMENTS USING THE
; COTUNES TAKEN WITH THE LASER AND ASSUMING THAT THE WAMETER HAS A
; CONSTANT ERROR THROUGHOUT THE WAVELENGTH RANGE, OF ABOUT 0.00405093 A
; FSR NB is fixed, we fit for FSR WB and FSR E1

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

PRO HMI_FSR

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

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------
nseq        = 10                ; number of positions in the sequence
nl0         = 27;30  ; number of wavelengths

; MEASURED INTENSITIES
;--------------------------------------------
 list= STRARR(nl0)
 lam0= DBLARR(nl0)

 list[0]= 'listLaser070223_78715'
 list[1]= 'listLaser070223_78754'
 list[2]= 'listLaser070223_78819' ; wavelength jumped a couple of times
 list[3]= 'listLaser070223_78845'
 list[4]= 'listLaser070223_78912'
 list[5]= 'listLaser070223_78938'
 list[6]= 'listLaser070223_78964' ; wavelength had a very short jump
 list[7]= 'listLaser070223_78990'
 list[8]= 'listLaser070223_79016'
 list[9]= 'listLaser070223_79094'
 list[10]= 'listLaser070223_79120'
 list[11]= 'listLaser070223_79146'
 list[12]= 'listLaser070223_79198'
 list[13]= 'listLaser070223_79237'
;list[]= 'listLaser070223_79366' ; etalon A and logs disagree for wavelength
 list[14]= 'listLaser070223_79405'
 list[15]= 'listLaser070223_79431'
;list[]= 'listLaser070223_79444' ; a few short wavelength jumps
 list[16]= 'listLaser070223_79483'
 list[17]= 'listLaser070223_79522'
 list[18]= 'listLaser070223_79548'
 list[19]= 'listLaser070223_79574'
 list[20]= 'listLaser070223_79600'
 list[21]= 'listLaser070223_79652'
 list[22]= 'listLaser070223_79704'
 list[23]= 'listLaser070223_79730'
 list[24]= 'listLaser070223_79782' ; wavelength may have jumped a bit
 list[25]= 'listLaser070223_79808'
 list[26]= 'listLaser070223_79834'
;list[]= 'listLaser070223_79431'
;list[]= 'listLaser070223_79860'
 lam0[0]= 6173.150d0
 lam0[1]= 6172.900d0
 lam0[2]= 6172.650d0
 lam0[3]= 6172.400d0
 lam0[4]= 6172.150d0
 lam0[5]= 6171.900d0
 lam0[6]= 6171.400d0
 lam0[7]= 6170.899d0
 lam0[8]= 6170.399d0
 lam0[9]= 6169.900d0
 lam0[10]= 6169.399d0 
 lam0[11]= 6168.900d0
 lam0[12]= 6168.399d0 
 lam0[13]= 6173.339d0;6173.342d0 
 lam0[14]= 6173.650d0
 lam0[15]= 6173.899d0;900d0
 lam0[16]= 6174.150d0 
 lam0[17]= 6174.400d0
 lam0[18]= 6174.650d0
 lam0[19]= 6174.900d0
 lam0[20]= 6175.400d0
 lam0[21]= 6175.907d0
 lam0[22]= 6176.400d0
 lam0[23]= 6176.900d0
 lam0[24]= 6177.400d0
 lam0[25]= 6177.901d0
 lam0[26]= 6178.400d0
;lam0= 6173.342d0 
;lam0= 6173.900d0
;lam0= 6173.340d0

;list[0]= 'listLaser070222_78702'
;list[1]= 'listLaser070222_78728'
;list[2]= 'listLaser070222_78767'
;list[3]= 'listLaser070222_78832'
;list[4]= 'listLaser070222_78858'
;list[5]= 'listLaser070222_78899'
;list[6]= 'listLaser070222_78925'
;list[7]= 'listLaser070222_78951'
;list[8]= 'listLaser070222_78977'
;list[9]= 'listLaser070222_79003'
;list[10]= 'listLaser070222_79029'
;list[11]= 'listLaser070222_79081'
;list[12]= 'listLaser070222_79107'
;list[13]= 'listLaser070222_79133'
;list[14]= 'listLaser070222_79159'
;list[15]= 'listLaser070222_79185'
;list[16]= 'listLaser070222_79263'
;list[17]= 'listLaser070222_79418'
;list[18]= 'listLaser070222_79457'
;list[19]= 'listLaser070222_79496'
;list[20]= 'listLaser070222_79535'
;list[21]= 'listLaser070222_79561'
;list[22]= 'listLaser070222_79587'
;list[23]= 'listLaser070222_79613'
;list[24]= 'listLaser070222_79665'
;list[25]= 'listLaser070222_79717'
;list[26]= 'listLaser070222_79743'
;list[27]= 'listLaser070222_79795'
;list[28]= 'listLaser070222_79821'
;list[29]= 'listLaser070222_79847'
;lam0[0] = 6173.403d0
;lam0[1] = 6173.150d0
;lam0[2] = 6172.900d0
;lam0[3] = 6172.647d0
;lam0[4] = 6172.400d0
;lam0[5] = 6172.399d0 
;lam0[6] = 6172.150d0
;lam0[7] = 6171.900d0
;lam0[8] = 6171.400d0
;lam0[9] = 6170.899d0
;lam0[10] = 6170.398d0
;lam0[11] = 6169.909d0
;lam0[12] = 6169.900d0
;lam0[13] = 6169.399d0
;lam0[14] = 6168.900d0
;lam0[15] = 6168.399d0
;lam0[16] = 6173.342d0
;lam0[17] = 6173.650d0
;lam0[18] = 6173.900d0
;lam0[19] = 6174.150d0
;lam0[20] = 6174.400d0
;lam0[21] = 6174.650d0
;lam0[22] = 6174.900d0
;lam0[23] = 6175.400d0
;lam0[24] = 6175.906d0
;lam0[25] = 6176.400d0
;lam0[26] = 6176.900d0
;lam0[27] = 6177.400d0
;lam0[28] = 6177.901d0
;lam0[29] = 6178.400d0

lam0= lam0-lamref+0.0040593d0
I0g         = FLTARR(nl0)

FOR iii=0,nl0-1 DO BEGIN

    nx          = 256
    RESTORE,'SEQUENCE_'+STRTRIM(list[iii],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
    Inten   = imx[*,*,2:11] ; for data after February 2007 obtained with readimages2.pro
    Inten   = REFORM(REBIN(Inten,1,1,nseq))


; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------

    Bg          = DBLARR(3)     ; contrasts of the tunable elements
    Phig        = DBLARR(3)     ; relative phases
    FSRg        = FLTARR(2)
    thresh      = TOTAL(Inten)*0.8
    nparam      = 3             ; number of parameters we fit for
    Residual    = DBLARR(nseq)  ; residuals of the least-squares fit
    dIntendFSR1 = DBLARR(nseq) ; derivatives to compute the Jacobian matrix
    dIntendFSR2 = DBLARR(nseq)
    dIntendIc   = DBLARR(nseq)
    tuning      = DBLARR(3,nseq) ; tuning sequence
    maxsteps    = 105 ; maximum steps allowed for the iterative Least-Squares algorithm
    history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err         = DBLARR(maxsteps)  ; error estimate at each step    


; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

    FSR         = DBLARR(7)     ; FSR in Angstrom
    FSR[0]      = 0.17110619;0.17076972d0  ; for the narrow-band Michelson
    FSR[1]      = 0.34158512d0  ; for the broad-band  Michelson
    FSR[2]      = 0.693d0       ; for E1
    FSR[3]      = 1.405d0       ; for E2
    FSR[4]      = 2.779d0       ; for E3
    FSR[5]      = 5.682d0       ; for E4
    FSR[6]      = 11.354d0      ; for E5

    FSRg[*]     = [FSR[1],FSR[2]]

; DEFINITION OF THE CO-TUNE SEQUENCE
; co-tuning in the range [-2 FSR[E1],+2 FSR[E1]]
; TABLE PROVIDED BY JESPER
;----------------------------------------------------------

    tuningJS      = DBLARR(nseq,3) ; co-tune sequence of Jesper
    
    tuningJS[0,*] = [ 7, 53, 78]
    tuningJS[1,*] = [13, 41, 54]
    tuningJS[2,*] = [19, 29, 90]
    tuningJS[3,*] = [25, 77, 66]
    tuningJS[4,*] = [31, 65,102]
    tuningJS[5,*] = [37, 53, 78]
    tuningJS[6,*] = [43, 41, 54]
    tuningJS[7,*] = [49, 29, 90]
    tuningJS[8,*] = [55, 77, 66]
    tuningJS[9,*] = [61, 65,102]
    
    tuningJS      = REVERSE(tuningJS,2)
    FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]
    
; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVDED BY LOCKHEED-MARTIN
;-----------------------------------------------

; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
    RESTORE,'frontwindow.bin'
    blocker0     = transmission/100.d0
    q            = READFITS('blocker11.fits',/SILENT) 
    blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+2.-lamref,wavelength*10.d0-lamref)
    

; RELATIVE PHASES OBTAINED BY HMI_sun1.pro
;------------------------------------------------

   ;Phig[0:2] = [170.94679,-2.8646375,-158.36210] ; with 79916
   ;Bg[0:2]   = [0.97,0.986,0.958]

   ;Phig[0:2] = [172.02769,-2.6257597,-159.19239] ; with 79306 ; CALMODE
   ;Bg[0:2] = [0.9711,0.9908,0.9697]

    Phig[0:2] = [173.46944,-2.1415770,-163.00307] ; with 79336 ; OBSMODE
    Bg[0:2] = [0.969337,0.986574,0.958388]


    Phig[*]   = Phig[*] / 180.d0 * !dpi

; LEAST-SQUARES ALGORITHM
;-------------------------------------------------

    WINDOW,0,RETAIN=2
    

    converg  = 0
    jj       = 0
    I0g[iii] = thresh
    
    WHILE (converg EQ 0) DO BEGIN
        

        FOR j=0,nseq-1 DO BEGIN
            
            profileg    = (1.d0+Bg[0]*COS(dpi/FSR[0]*lam0[iii]+Phig[0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[1]*COS(dpi/FSRg[0]*lam0[iii]+Phig[1]+tuningJS[j,1]))/2.d0*(1.d0+Bg[2]*COS(dpi/FSRg[1]*lam0[iii]+Phig[2]+tuningJS[j,2]))/2.d0

            profileg1p  = (1.d0+Bg[1]*COS(dpi/(FSRg[0]+0.002d0)*lam0[iii]+Phig[1]+tuningJS[j,1]))/2.d0*(1.d0+Bg[0]*COS(dpi/FSR[0]*lam0[iii]+Phig[0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[2]*COS(dpi/FSRg[1]*lam0[iii]+Phig[2]+tuningJS[j,2]))/2.d0
            profileg1m  = (1.d0+Bg[1]*COS(dpi/(FSRg[0]-0.002d0)*lam0[iii]+Phig[1]+tuningJS[j,1]))/2.d0*(1.d0+Bg[0]*COS(dpi/FSR[0]*lam0[iii]+Phig[0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[2]*COS(dpi/FSRg[1]*lam0[iii]+Phig[2]+tuningJS[j,2]))/2.d0
            profileg2p  = (1.d0+Bg[2]*COS(dpi/(FSRg[1]+0.002d0)*lam0[iii]+Phig[2]+tuningJS[j,2]))/2.d0*(1.d0+Bg[0]*COS(dpi/FSR[0]*lam0[iii]+Phig[0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[1]*COS(dpi/FSRg[0]*lam0[iii]+Phig[1]+tuningJS[j,1]))/2.d0
            profileg2m  = (1.d0+Bg[2]*COS(dpi/(FSRg[1]-0.002d0)*lam0[iii]+Phig[2]+tuningJS[j,2]))/2.d0*(1.d0+Bg[0]*COS(dpi/FSR[0]*lam0[iii]+Phig[0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[1]*COS(dpi/FSRg[0]*lam0[iii]+Phig[1]+tuningJS[j,1]))/2.d0


           
            Residual[j]   = Inten[j] - profileg*I0g[iii]
        ;PRINT,profileg*I0g,Inten[j]
            dIntendIc[j]  =  profileg
            dIntendFSR1[j]= (profileg1p-profileg1m)/0.004d0*I0g[iii]
            dIntendFSR2[j]= (profileg2p-profileg2m)/0.004d0*I0g[iii]

        ENDFOR
              ; Jacobian matrix
    ;read,p
        Jac  = DBLARR(nparam,nseq)
        FOR i= 0,nseq-1 DO BEGIN
            
            Jac[0,i]  = dIntendIc[i]
            Jac[1,i]  = dIntendFSR1[i]
            Jac[2,i]  = dIntendFSR2[i]
            
        ENDFOR
        
        LA_SVD,Jac,W,U,V,/DOUBLE
        
        
        Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
        
        temp   = TRANSPOSE(Dx)/[I0g[iii],FSRg[0],FSRg[1]]
        err[jj]= MAX(ABS(temp))
        
        IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7)     THEN converg = 2 
        IF(err[jj] LE 1.d-7)                          THEN converg = 1
        
        IF(converg EQ 2) THEN BEGIN
            j     = WHERE(err EQ MIN(err))
            I0g[iii]   = history[j[0],0]
            FSRg[0]    = history[j[0],1]
            FSRg[1]    = history[j[0],2]
            
        ENDIF
        
        IF(converg NE 2) THEN BEGIN
            I0g[iii] = I0g[iii] + Dx[0]
            IF(I0g[iii] LT 0.002*thresh OR I0g[iii] GT 100.0*thresh) THEN I0g[iii] = thresh
            FSRg[0] = FSRg[0] + Dx[1]
            IF(ABS(FSRg[0]-FSR[1])/FSR[1] GT  1.d0) THEN FSRg[0]  = FSR[1] ;-0.030
            FSRg[1] = FSRg[1] + Dx[2]
            IF(ABS(FSRg[1]-FSR[2])/FSR[2] GT  1.d0) THEN FSRg[1]  = FSR[2] ;-0.030
        ENDIF
        
        history[jj,*] = [I0g[iii],FSRg[0],FSRg[1]]
        
        jj = jj+1
        
    ENDWHILE
    ;PRINT,'PARAMETERS='
    ;PRINT,'LISTNAME',list[iii]
    ;PRINT,'Doppler Shift=',lam0g[iii]+lamref
    ;PRINT,'I0g=',I0g[iii]            
    ;PRINT,'lateconv=',converg
    
; RECONSTRUCTION
    Intenrecons=FLTARR(nseq)
    FOR j=0,nseq-1 DO BEGIN
                                ;profileg  = INTERPOL(blocker0,wavelength*10.d0-lamref,lam0g)
        profileg = (1.d0+Bg[0]*COS(dpi/FSR[0]*lam0[iii]+Phig[0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[1]*COS(dpi/FSRg[0]*lam0[iii]+Phig[1]+tuningJS[j,1]))/2.d0*(1.d0+Bg[2]*COS(dpi/FSRg[1]*lam0[iii]+Phig[2]+tuningJS[j,2]))/2.d0
        Intenrecons[j]  = profileg*I0g[iii]
    ENDFOR
    PRINT,list[iii],FSRg[0],FSRg[1],SQRT(TOTAL((Intenrecons-Inten)^2.))/SQRT(TOTAL(Intenrecons^2.0)),I0g[iii]
    
    PLOT,Inten,xst=1
    OPLOT,Intenrecons,linestyle=2,thick=2
    READ,pause

ENDFOR



END

