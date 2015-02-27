; code to test HMI_laser2.pro

PRO testHMI_laser2

lamref = 6173.3433d0
dpi    = 2.d0*!dpi

RESTORE,'frontwindow.bin'
blocker0     = transmission/100.d0;INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
lam          = wavelength*10.d0-lamref
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+2.d0-lamref,lam) ; I center the profile

nx=64
ny=64
anglim= 970.
nseq        = 10               ; number of positions in the co-tune sequence
nl0         = 27               ; number of wavelength positions for the dye laser
xcenter     = nx/2
ycenter     = nx/2
distance    = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx
FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172457d0         ; for the narrow-band Michelson
FSR[1]      = 0.344242d0         ; for the broad-band  Michelson
FSR[2]      = 0.690d0            ; for E1
FSR[3]      = 1.405d0            ; for E2
FSR[4]      = 2.779d0            ; for E3
FSR[5]      = 5.682d0            ; for E4
FSR[6]      = 11.354d0           ; for E5
Inten   = DBLARR(nx,ny,nseq*nl0)
lam0    = FLTARR(nl0)
lam0[0]=6173.4
lam0[1]=6173.65
lam0[2]=6173.90
lam0[3]=6174.15
lam0[4]=6174.4
lam0[5]=6174.65
lam0[6]=6174.9
lam0[7]=6175.4
lam0[8]=6175.9
lam0[9]=6176.4
lam0[10]=6176.9
lam0[11]=6177.4
lam0[12]=6177.9
lam0[13]=6178.4
lam0[14]=6173.15
lam0[15]=6172.9
lam0[16]=6172.65
lam0[17]=6172.4
lam0[18]=6172.15
lam0[19]=6171.9
lam0[20]=6171.4
lam0[21]=6170.9
lam0[22]=6170.4
lam0[23]=6169.9
lam0[24]=6169.4
lam0[25]=6168.9
lam0[26]=6168.4
lam0    = lam0-lamref

tuningJS      = DBLARR(nseq,3)  ; co-tune sequence of Jesper

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

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"

RESTORE,'RESULTS/RESULTS_listLaser070220_74548_256.BIN' ;OBSMODE (WARNING CANCEL nx)
a=WHERE(FINITE(Bg0) EQ 0)
IF (a[0] NE -1) THEN Bg0[a]  =1.d0
a=WHERE(FINITE(Phig0) EQ 0)
IF (a[0] NE -1) THEN Phig0[a]=0.d0

nx=64
ny=64
Bg0 = REBIN(Bg0,nx,nx,3)
Phig0=REBIN(Phig0,nx,nx,3)

Bg[*,*,0:2]   = Bg0[*,*,*]         ; contrasts and phases of the tuning elements
Phig[*,*,0:2] = Phig0[*,*,*]


seed        = 1l
FOR jjj=0,ny-1 DO BEGIN
    FOR iii=0,nx-1 DO BEGIN
        IF(distance[iii,jjj] LE anglim) THEN BEGIN
            Bg[iii,jjj,3:6] = RANDOMU(seed,4)*0.15+0.85
            Phig[iii,jjj,3:6] = RANDOMU(seed,4)*240.-120.
            I0g[iii,jjj] = RANDOMU(seed)*4000.+8000.
        ENDIF
    ENDFOR
ENDFOR
Phig = Phig*!dpi/180.d0

FOR jjj=0,ny-1 DO BEGIN
 FOR iii=0,nx-1 DO BEGIN
     IF(distance[iii,jjj] LE anglim) THEN BEGIN
         FOR j=0,nseq-1 DO BEGIN
             FOR ii=0,nl0-1 DO BEGIN
                 
                 profileg  = INTERPOL(blocker0,lam,lam0[ii])
                 FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]$
                                                                            +Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                 FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]$
                                                                            +Phig[iii,jjj,i]))/2.d0
                 Inten[iii,jjj,ii*nseq+j] = profileg*I0g[iii,jjj]
                 Inten[iii,jjj,ii*nseq+j] = Inten[iii,jjj,ii*nseq+j]+RANDOMN(seed)*SQRT(Inten[iii,jjj,ii*nseq+j])   
             ENDFOR
         ENDFOR
     ENDIF
 ENDFOR
ENDFOR

SAVE,Bg,Phig,I0g,Inten,file='simul.bin'

READ,pause

END
