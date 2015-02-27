; COMPUTE DOPPLERGRAMS FROM LINE SCANS


PRO dopplergram

nx      =  1024;512
;N       =  20.d0 ; number of filters
;dl      =  0.069d0 ;separation between filters
N       = 6.d0
dl      = 0.068d0
LCP     = [0,4, 8,12,16,20]
RCP     = [2,6,10,14,18,22]
Nelem   = 24 ; number of filtergrams in a sequence
Nseq1   = 32.
Nseq2   = 53.
GOTO,plot

imx     = DBLARR(nx,nx,N)

cste    = 2.99792458d8/6173.3433d0/(2.d0*!dpi)*N*dl
richard = 2d0*!dpi/(N*dl)

;list = 'listScan_060224_172051' ; no problem
;list = 'listScan_060225_171619' ; one bad exposure
;list = 'listSun_060630_213249'  ; no problem
;list = 'listSun060630_211004'   ; no problem

;READIMAGES,list,images,headers,time,N+2,wgood,nbin=nx,iover=iover
;CLEAN,images,imx,headers
;SAVE,imx,time,wgood,FILE='SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;read,pause

;RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
; IF THERE IS NO BAD EXPOSURE OR OVERSCAN:
;IF list EQ 'listSun_060630_213249' THEN BEGIN
;    Inten        = FLTARR(nx,nx,N)
;    Inten[*,*,0] = imx[*,*,7]
;    Inten[*,*,1] = imx[*,*,9]
;    Inten[*,*,2] = imx[*,*,11]
;    Inten[*,*,3] = imx[*,*,13]
;    Inten[*,*,4] = imx[*,*,15]
;;Inten = imx[*,*,6:16]
;ENDIF
;
;IF list EQ 'listScan_060224_172051'THEN BEGIN
;    Inten        = FLTARR(nx,nx,N)
;    Inten[*,*,0] = imx[*,*,2]
;    Inten[*,*,1] = imx[*,*,4]
;    Inten[*,*,2] = imx[*,*,6]
;    Inten[*,*,3] = imx[*,*,8]
;    Inten[*,*,4] = imx[*,*,10]
;ENDIF

;RESTORE,'CPT/Observables_706315_708589.BIN'
;RESTORE,'CPT/Observables_706315.bin'
RESTORE,'CPT/Observables_706315_708589_HiRes.BIN'
FOR i=0,2274 DO images[*,*,i]   = images[*,*,i]-images[*,*,777] ; remove dark frame
lam0g    = DBLARR(nx,nx,Nseq1+Nseq2)
linewidth= DBLARR(nx,nx,Nseq1+Nseq2)
linedepth= DBLARR(nx,nx,Nseq1+Nseq2)

beginning= 0
FOR j=0,Nseq1-1 DO BEGIN
    PRINT,j
    Inten=images[*,*,j*24.+beginning:(j+1)*24.-1.+beginning]
    FOR jjj=0,nx-1 DO BEGIN
                                ;PRINT,jjj
        FOR iii=0,nx-1 DO BEGIN
            
            f1LCP = 0.d0
            f2LCP = 0.d0
            f1RCP = 0.d0
            f2RCP = 0.d0
            FOR i=0,N-1 DO BEGIN
                f1LCP = f1LCP + exp(   DCOMPLEX(0,1)*2.d0*!dpi*i/N)*Inten[iii,jjj,LCP[i]] ;1st Fourier component
                f2LCP = f2LCP + exp(2.*DCOMPLEX(0,1)*2.d0*!dpi*i/N)*Inten[iii,jjj,LCP[i]] ;2nd Fourier component
                f1RCP = f1RCP + exp(   DCOMPLEX(0,1)*2.d0*!dpi*i/N)*Inten[iii,jjj,RCP[i]] ;1st Fourier component
                f2RCP = f2RCP + exp(2.*DCOMPLEX(0,1)*2.d0*!dpi*i/N)*Inten[iii,jjj,RCP[i]] ;2nd Fourier component
            ENDFOR
            
            lam0g[iii,jjj,j]     = (ATAN(IMAGINARY(f1LCP),DOUBLE(f1LCP))*cste+ATAN(IMAGINARY(f1RCP),DOUBLE(f1RCP))*cste);/2.d0
            linewidth[iii,jjj,j] = (SQRT(ABS(f1LCP)/ABS(f2LCP)*2./3.)/richard*SQRT(2.d0)+SQRT(ABS(f1RCP)/ABS(f2RCP)*2./3.)/richard*SQRT(2.d0))/2.d0 ; sqrt(2.) because of my definition of linewidth
            linedepth[iii,jjj,j] = (TOTAL(Inten[iii,jjj,LCP])-ABS(f1LCP)+TOTAL(Inten[iii,jjj,RCP])-ABS(f1RCP))/2.d0

        ENDFOR
    ENDFOR
    TVIM,lam0g[*,*,j],/scale
ENDFOR
beginning= 835
FOR j=0,Nseq2-1 DO BEGIN
    PRINT,j
    Inten=images[*,*,j*24.+beginning:(j+1)*24.-1.+beginning]
    FOR jjj=0,nx-1 DO BEGIN
                                ;PRINT,jjj
        FOR iii=0,nx-1 DO BEGIN
            
            f1LCP = 0.d0
            f2LCP = 0.d0
            f1RCP = 0.d0
            f2RCP = 0.d0
            FOR i=0,N-1 DO BEGIN
                f1LCP = f1LCP + exp(   DCOMPLEX(0,1)*2.d0*!dpi*i/(N-1.))*Inten[iii,jjj,LCP[i]] ;1st Fourier component
                f2LCP = f2LCP + exp(2.*DCOMPLEX(0,1)*2.d0*!dpi*i/(N-1.))*Inten[iii,jjj,LCP[i]] ;2nd Fourier component
                f1RCP = f1RCP + exp(   DCOMPLEX(0,1)*2.d0*!dpi*i/(N-1.))*Inten[iii,jjj,RCP[i]] ;1st Fourier component
                f2RCP = f2RCP + exp(2.*DCOMPLEX(0,1)*2.d0*!dpi*i/(N-1.))*Inten[iii,jjj,RCP[i]] ;2nd Fourier component
            ENDFOR
            
            lam0g[iii,jjj,j+Nseq1]     = (ATAN(IMAGINARY(f1LCP),DOUBLE(f1LCP))*cste+ATAN(IMAGINARY(f1RCP),DOUBLE(f1RCP))*cste);/2.d0
            linewidth[iii,jjj,j+Nseq1] = (SQRT(ABS(f1LCP)/ABS(f2LCP)*2./3.)/richard*SQRT(2.d0)+SQRT(ABS(f1RCP)/ABS(f2RCP)*2./3.)/richard*SQRT(2.d0))/2.d0 ; sqrt(2.) because of my definition of linewidth
            linedepth[iii,jjj,j+Nseq1] = (TOTAL(Inten[iii,jjj,LCP])-ABS(f1LCP)+TOTAL(Inten[iii,jjj,RCP])-ABS(f1RCP))/2.d0

        ENDFOR
    ENDFOR
    TVIM,lam0g[*,*,j+Nseq1],/scale
ENDFOR

a=where(finite(linewidth) eq 0)
linewidth[a]=0.07
SAVE,lam0g,file='dopplergram_HiRes.bin' ;dopplergram.bin ;magnetogram.bin; linewidth.bin
plot:
RESTORE,'dopplergram_HiRes.bin'

;center=[131,131] ;!!!WARNING!!! listScan_060224_172051
;center=[262,262]
;center=[133,128] ;!!!WARNING!!! listScan_060225_171619
;center   = [2077,1985];[259,247]
center   = [519,496]
distance = SHIFT(DIST(nx,nx),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center

a = WHERE(distance LT 920.d0,COMPLEMENT=b) ;!!!WARNING!!!
distance[a]=1.0d0
distance[b]=0.0d0
FOR i=0,Nseq1+Nseq2-1 DO lam0g[*,*,i] = lam0g[*,*,i] * distance
lam0g[*,*,10]=(lam0g[*,*,9]+lam0g[*,*,11])/2.d0 ; to correct a shitty dopplergram
maxi=0.
mini=100000.
FOR i=0,Nseq1+Nseq2-1 DO BEGIN
    temp = REFORM(lam0g[*,*,i])
    lam0g[*,*,i] = lam0g[*,*,i] - MEAN(temp[a])
    temp = REFORM(lam0g[*,*,i])
    IF (min(temp[a]) LE mini) THEN mini=MIN(temp[a])
    IF(max(temp[a]) GE maxi) THEN maxi=MAX(temp[a])
ENDFOR
LOADCT,3

PRINT,'PRODUCE A MPEG MOVIE'
; CODE FROM JESPER TO PRODUCE MPEGS
; 55 makes file twice as large, but it looks a lot better.
;nsmax=100
;nsc=((nsmax-1)/8+1)*8
;nbig=nx+2*nsc
;nbig = 360.
;mpegID = MPEG_OPEN([nbig,nbig],quality=100)
;TVLCT, r, g, b, /Get
;temp = BYTARR(3,nx,nx)
;FOR i=0,Nseq1+Nseq2-1 DO BEGIN
;    temp[0,*,*] = r(REFORM(lam0g[*,*,i]))
;    temp[1,*,*] = g(REFORM(lam0g[*,*,i]))
;    temp[2,*,*] = b(REFORM(lam0g[*,*,i]))
;    MPEG_PUT,mpegID,image=temp,frame=i,/color
;ENDFOR
;MPEG_SAVE, mpegID, FILENAME='dopplergram.mpg'
;MPEG_CLOSE, mpegID


SET_PLOT,'ps'

FOR i=0,Nseq1+Nseq2-1 DO BEGIN 
    IF i LE 9 THEN DEVICE,file='yo0'+STRTRIM(STRING(i),1)+'.ps',xoffset=0,yoffset=0,xsize=40,ysize=40,/color,bits=24 ELSE DEVICE,file='yo'+STRTRIM(STRING(i),1)+'.ps',xoffset=0,yoffset=0,xsize=40,ysize=40,/color,bits=24
    TVIM,lam0g[*,*,i],range=[mini,maxi]
    DEVICE,/close
ENDFOR


READ,pause


;SET_PLOT,'x'
!P.MULTI=0
;WINDOW,0,RETAIN=2,xsize=800,ysize=800
SET_PLOT,'ps'
DEVICE,FILE='yo.ps',BITS=24,xsize=22,ysize=27,xoffset=0,yoffset=0,/color
LOADCT,3
TVIM,lam0g,/scale,stit='Doppler Velocity (m/s)'
DEVICE,/CLOSE
SET_PLOT,'x'

;SAVE,lam0g,FILE='RESULTS/RESULTS_SCAN_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

READ,PAUSE



END
