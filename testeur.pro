PRO moncul,X,A,F,pder

RESTORE,'phase.bin',VERBOSE=0

nseq        = 10
F           = DBLARR(nseq)
lamref      = 6173.3433d0
dpi         = 2.d0*!dpi
nlam        = 20000            ; number of wavelengths
dlam        = 1.d0/1.5d3       ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
RESTORE,'frontwindow.bin' ,VERBOSE=0
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q           = READFITS('blocker11.fits',/SILENT)      ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam);,/LSQUADRATIC) ; I center the profile

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172457         ; for the narrow-band Michelson
FSR[1]      = 0.344242         ; for the broad-band  Michelson
FSR[2]      = 0.702            ; for E1
FSR[3]      = 1.405            ; for E2
FSR[4]      = 2.779            ; for E3
FSR[5]      = 5.682            ; for E4
FSR[6]      = 11.354           ; for E5

tuningJS     = DBLARR(10,3) 
tuningJS[0,*]=[  111 , 100 ,   67]
tuningJS[1,*]=[  117 ,  88 ,   43]
tuningJS[2,*]=[  123 ,  76 ,   79]
tuningJS[3,*]=[  129 , 124 ,   55]
tuningJS[4,*]=[  135 , 112 ,   91]
tuningJS[5,*]=[   81 , 100 ,   67]
tuningJS[6,*]=[   87 ,  88 ,   43]
tuningJS[7,*]=[   93,   76 ,   79]
tuningJS[8,*]=[   99 , 124 ,   55]
tuningJS[9,*]=[  105 , 112 ,   91]

tuningJS=REVERSE(tuningJS,2)
FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]

profileg=INTERPOL(blocker0,lam,lam1+A[0])*ABS(A[1])
FOR ii=3,6 DO profileg = profileg * (1.d0+contrast[ii]*COS(dpi/FSR[ii]*(A[0]+lam1)+phase[ii]))/2.d0

FOR j=0,nseq-1 DO BEGIN
    profilef=profileg
    FOR ii=0,2 DO profilef = profilef * (1.d0+contrast[ii]*COS(dpi/FSR[ii]*(A[0]+lam1)+phase[ii]+tuningJS[j,ii]))/2.d0 
    F[j] = profilef
ENDFOR        


END


PRO testeur

nx          = 32              ; number of rows
ny          = 32              ; number of columns

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"
nseq        = 10
nlam        = 40000            ; number of wavelengths
dlam        = 1.d0/2.d3       ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
nl0         = 27
lamref      = 6173.3433d0
dpi         = 2.d0*!dpi
seed        = 1l
Inten       = DBLARR(nx,ny,nl0*nseq)

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172457         ; for the narrow-band Michelson
FSR[1]      = 0.344242         ; for the broad-band  Michelson
FSR[2]      = 0.702            ; for E1
FSR[3]      = 1.405            ; for E2
FSR[4]      = 2.779            ; for E3
FSR[5]      = 5.682            ; for E4
FSR[6]      = 11.354           ; for E5

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

;blocker0     = DBLARR(nlam)+1.d0
q           = READFITS('blocker11.fits')      ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam);,/LSQUADRATIC) ; I center the profile


lam0    = DBLARR(nl0)
lam0[0] = 6173.463d0
lam0[1] = 6173.647d0
lam0[2] = 6173.277d0
lam0[3] = 6173.904d0
lam0[4] = 6174.150d0
lam0[5] = 6174.400d0
lam0[6] = 6174.552d0
lam0[7] = 6174.913d0
lam0[8] = 6175.399d0 ; some bad images
lam0[9] = 6175.900d0
lam0[10]= 6176.906d0
lam0[11]= 6177.398d0
lam0[12]= 6177.901d0
lam0[13]= 6173.400d0
lam0[14]= 6173.147d0
lam0[15]= 6172.900d0
lam0[16]= 6172.668d0
lam0[17]= 6172.276d0 ;wavemeter jumped sometimes
lam0[18]= 6172.150d0
lam0[19]= 6171.902d0
lam0[20]= 6171.399d0
lam0[21]= 6170.900d0
lam0[22]= 6170.400d0
lam0[23]= 6169.896d0
lam0[24]= 6169.399d0
lam0[25]= 6168.899d0
lam0[26]= 6168.400d0

relI    = DBLARR(nl0)
relI[0] = 85.
relI[1] = 80.
relI[2] = 17.
relI[3] = 102.
relI[4] = 37.
relI[5] = 30.
relI[6] = 20.
relI[7] = 14.
relI[8] = 22.
relI[9] = 11.
relI[10]= 22.
relI[11]= 12.
relI[12]= 20.
relI[13]= 60.
relI[14]= 31.
relI[15]= 104.
relI[16]= 20.
relI[17]= 15.
relI[18]= 22.
relI[19]= 35.
relI[20]= 25.
relI[21]= 15.
relI[22]= 30.
relI[23]= 26.
relI[24]= 19.
relI[25]= 37.
relI[26]= 42.

tuningJS      = DBLARR(nseq,3)
tuningJS[0,*]=[  111 , 100 ,   67]
tuningJS[1,*]=[  117 ,  88 ,   43]
tuningJS[2,*]=[  123 ,  76 ,   79]
tuningJS[3,*]=[  129 , 124 ,   55]
tuningJS[4,*]=[  135 , 112 ,   91]
tuningJS[5,*]=[   81 , 100 ,   67]
tuningJS[6,*]=[   87 ,  88 ,   43]
tuningJS[7,*]=[   93,   76 ,   79]
tuningJS[8,*]=[   99 , 124 ,   55]
tuningJS[9,*]=[  105 , 112 ,   91]

tuningJS=REVERSE(tuningJS,2)
FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]

; DATA

shiftr  = DBLARR(nl0);RANDOMU(seed,nl0)*0.10-0.05   !!!!!!!!!!!!!!!!!!!!!
Phig0   =(RANDOMU(seed,nx,ny,3)*50.-25.)*!dpi/180.d0
Phig    =(RANDOMU(seed,nx,ny,7)*160.-80.)*!dpi/180.d0
Phig[*,*,0:2] = Phig0
Bg      = RANDOMU(seed,nx,ny,7)*0.1+0.9
I0g     = RANDOMU(seed,nx,ny)*5000.+5000.


; ARTIFICIAL INTENSITIES

lam0[*]    = lam0[*]-lamref

FOR iii=0,nx-1 DO BEGIN
    FOR jjj=0,ny-1 DO BEGIN
        FOR ii=0,nl0-1 DO BEGIN
            FOR j=0,nseq-1 DO BEGIN
                profileg                 = INTERPOL(blocker0,lam,(lam0[ii]+shiftr[ii]))
                profileg                 = profileg * I0g[iii,jjj]
                FOR i=0,2 DO profileg    = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*(lam0[ii]+shiftr[ii])+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                FOR i=3,6 DO profileg    = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*(lam0[ii]+shiftr[ii])+Phig[iii,jjj,i]))/2.d0
                Inten[iii,jjj,ii*nseq+j] = profileg
            ENDFOR
        ENDFOR
    ENDFOR        
ENDFOR

imx=Inten
FOR i=0,nl0-1 DO BEGIN
FOR j=0,nseq-1 DO BEGIN
imx[*,*,i*nseq+j]=imx[*,*,i*nseq+j]*relI[i]/MEAN(relI)
ENDFOR
ENDFOR
Bg0=Bg[*,*,0:2]
SAVE,Phig0,Bg0,Bg,Phig,FILE='testeur1.bin'
SAVE,imx,FILE='testeur2.bin'

;FOR i=0,nl0-1 DO IF( i NE 1 AND i NE 8 AND i NE 14 AND i NE 17 AND i NE 26) THEN PRINT,i,shiftr[i]
;FOR i=0,nl0-1 DO PRINT,i,shiftr[i]

estshift = DBLARR(11,11,nl0)

;FOR iii=12,18 DO BEGIN
;FOR jjj=12,18 DO BEGIN
;print,iii,jjj

;phase        = DBLARR(7)
;contrast     = DBLARR(7)+0.985d0
;phase[0:2]   = Phig[iii,jjj,0:2]
;contrast[0:2]= Bg[iii,jjj,0:2]

;FOR i=0,nl0-1 DO BEGIN
;    lam1 = lam0[i]
;    SAVE,phase,contrast,lam1,file='phase.bin'
;    intensity         = REFORM(Inten[iii,jjj,i*nseq:(i+1)*nseq-1])
;    A = [0.d0,5000.d0]
;    res =CURVEFIT(DINDGEN(10),intensity,DBLARR(10)+1.d0,A,FUNCTION_NAME='moncul',/NODERIV,ITMAX=1000,TOL=1.d-7)
;    estshift[iii-12,jjj-12,i] = A[0]
;    PRINT,A[0]
;ENDFOR

;ENDFOR
;ENDFOR

;FOR i=0,nl0-1 DO BEGIN
;    temp=REFORM(estshift[*,*,i])
;    a=WHERE(ABS(temp-MEAN(temp)) LT 2.*SIGMA(temp))
;    PRINT,i,MEAN(temp[a]),shiftr[i]
;ENDFOR


seuil = MAX(Inten)
Inten = Inten/MAX(Inten)
dblock= DERIV(lam,blocker0)


FOR iii=10,20 DO BEGIN
    FOR jjj=10,20 DO BEGIN
        print,iii,jjj
        
        phase        = DBLARR(7)
        contrast     = DBLARR(7)+0.985d0
        phase[0:2]     = Phig[iii,jjj,0:2]
        contrast[0:2]  = Bg[iii,jjj,0:2]
        
        ntun         = DBLARR(4)
        tun          = DBLARR(nseq,3)
        dntun        = DBLARR(4)
        dtun         = DBLARR(nseq,3)
        dIntendlam   = DBLARR(nseq)
        dIntendI     = DBLARR(nseq)
        Jac          = DBLARR(2,nseq)
        Residual     = DBLARR(nseq)
        
        FOR ii=0,nl0-1 DO BEGIN
            
            
            lam0g = lam0[ii]
            I0g   = 1.d0
            intensity         = REFORM(Inten[iii,jjj,ii*nseq:(ii+1)*nseq-1])
            
            
            converg = 0
            jj      = 0
            WHILE (converg EQ 0) DO BEGIN
                
                block  = INTERPOL(blocker0,lam,lam0g)
                ddblock= INTERPOL(dblock  ,lam,lam0g)
               ;blockp = INTERPOL(blocker0,lam,lam0g+dlam)
               ;blockm = INTERPOL(blocker0,lam,lam0g-dlam)

                FOR i=0,3 DO                   ntun[i]  = 0.5d0*(1.d0+contrast[i+3]*COS(dpi/FSR[i+3]*lam0g+phase[i+3]))
                FOR i=0,3 DO                   dntun[i] =-contrast[i+3]/2.d0*dpi/FSR[i+3]*SIN(dpi/FSR[i+3]*lam0g+phase[i+3])
                FOR j=0,nseq-1 DO FOR i=0,2 DO tun[j,i] = 0.5d0*(1.d0+contrast[i]*COS(dpi/FSR[i]*lam0g+phase[i]+tuningJS[j,i]))
                FOR j=0,nseq-1 DO FOR i=0,2 DO dtun[j,i]=-contrast[i]/2.d0*dpi/FSR[i]*SIN(dpi/FSR[i]*lam0g+phase[i]+tuningJS[j,i])
                
                FOR j=0,nseq-1 DO BEGIN
                    
                    dIntendlam[j] = I0g * block * (dtun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[1]*ntun[2]*ntun[3]+dtun[j,1]*tun[j,0]*tun[j,2]*ntun[0]*ntun[1]*ntun[2]*ntun[3]+dtun[j,2]*tun[j,0]*tun[j,1]*ntun[0]*ntun[1]*ntun[2]*ntun[3]+dntun[0]*tun[j,0]*tun[j,1]*tun[j,2]*ntun[1]*ntun[2]*ntun[3]+dntun[1]*tun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[2]*ntun[3]+dntun[2]*tun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[1]*ntun[3]+dntun[3]*tun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[1]*ntun[2])+I0g*tun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[1]*ntun[2]*ntun[3]*ddblock;(blockp-blockm)/(2.d0*dlam)
                    
                    dIntendI[j]   = tun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[1]*ntun[2]*ntun[3]*block
                    
                    Residual[j]   = intensity[j] - I0g * block * tun[j,0]*tun[j,1]*tun[j,2]*ntun[0]*ntun[1]*ntun[2]*ntun[3]
                    
                ENDFOR
                
                Jac[0,*] = dIntendlam[*]
                Jac[1,*] = dIntendI[*]
                
                LA_SVD,Jac,W,U,V,/DOUBLE
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
                
                err    = TRANSPOSE(Dx)/[lam0g,I0g]
                err = MAX(ABS(err))
                
                IF(jj EQ 25 OR err LT 1.d-7) THEN converg = 1
                
                lam0g = lam0g + Dx[0]
                IF(ABS(lam0g) GT 7.d0) THEN lam0g = lam0[ii]

                I0g   = I0g   + Dx[1]
                IF(I0g LT 0.d0) THEN I0g = 1.d0

                jj = jj+1
        
            ENDWHILE
    
            estshift[iii-10,jjj-10,ii] = lam0g-lam0[ii]
            PRINT,estshift[iii-10,jjj-10,ii],I0g*seuil
    
        ENDFOR

ENDFOR
ENDFOR

FOR i=0,nl0-1 DO BEGIN
    temp=REFORM(estshift[*,*,i])
    a=WHERE(ABS(temp-MEAN(temp)) LT 2.*SIGMA(temp))
    PRINT,i,MEDIAN(temp[a]),shiftr[i]
ENDFOR


READ,pause

END
