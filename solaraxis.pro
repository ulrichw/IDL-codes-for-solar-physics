; THIS PROGRAM DETERMINES THE SOLAR ROTATION AXIS FROM HMI DOPPLERGRAMS


PRO solaraxis

nangle   = 1000
Residual = FLTARR(nangle)
nx       = 128;512

RESTORE,'dopplerMDI.bin' ; MDI Dopplergram 1024*1024 taken on 2006/02/24 at 22:52 UT, from /tmp80/rick/sebastian
dopplerMDI = REBIN(dopplerMDI,nx,nx)
dopplerMDI = REVERSE(dopplerMDI,1) ; comparison HMI/MDI magnetogram shows a need to flip the image


; DOPPLERGRAMS FROM DETUNE27 SEQUENCES
;list = 'listSun060224_224805'
;list    = 'listSun060222_174615'
;list    = 'listSun060222_204505'
 list    = 'listSun060222_002209'
;list    = 'listSun060222_225521'



RESTORE,'RESULTS/RESULTS2_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'  ; Dopplergram obtained with HMI during sunlight test
;lam0g = lam0g - 621.36914 ; l.o.s. velocity
;lam0g = lam0g - 462.2428    
;lam0g = lam0g - 628.13942
 lam0g = lam0g -1025.3233
;lam0g = lam0g - 795.02078



;a=WHERE(lam0g EQ MAX(lam0g))
;lam0g[a] = 0.d0

; OTHER DOPPLERGRAMS FROM LINE SCANS
;list = 'listScan_060224_172051'
;list = 'listScan_060225_171619'
;RESTORE,'RESULTS/RESULTS_SCAN_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;lam0g=SHIFT(lam0g,-3,-3) ; for listScan_060224_172051
;lam0g=lam0g+770.09573   ; for listScan_060224_172051
;lam0g=SHIFT(lam0g,-5,0) ; for listScan_060225_171619
;lam0g=lam0g+0.25678151   ; for listScan_060225_171619

;lam0g = lam0g*2933.1359 ; FROM RICHARD WACHTER

angle = FINDGEN(nangle)*360.0/(nangle-1.)-180.0

distance = SHIFT(DIST(nx,nx),nx/2,nx/2)*0.5d0*4096.d0/nx
a=WHERE(distance LT 900.d0,COMPLEMENT=b) ;920;760.
distance[a]=1.d0
distance[b]=0.d0
lam0g[*,*] = lam0g[*,*]*distance
dopplerMDI = dopplerMDI*distance

FOR i=0,nangle-1 DO BEGIN
    PRINT,i
    d = ROT(lam0g,angle[i])
    residual[i] = TOTAL((d[a]-dopplerMDI[a])^2.0)
ENDFOR

a = WHERE(residual EQ MIN(residual))
PRINT,'ANGLE OF THE ROTATION AXIS = ',90.d0+angle[a]


SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xsize=22,ysize=27,xoffset=-1,yoffset=0,bits=24
!P.MULTI=[0,2,2]
TVIM,dopplerMDI,/scale,tit='MDI Dopplergram'
TVIM,lam0g,/scale,tit='HMI Dopplergram 02/24/2006 22:48 UT'
TVIM,ROT(lam0g,angle[a[0]]),tit='HMI Dopplergram rotated'
PLOT,angle,residual,yst=1,xtit='angle',ytit='residual'
DEVICE,/CLOSE


year = 2006
month= 2
dom  = 25;24;25;24
heure= 17;17;17;22.
minute=16;20;16;48.
seconde=19;51;19;5.
time = heure*3600.+minute*60.+seconde ; time when the detune sequence was taken in UT
;ra = rotation(year,month,dom,time) ; returns the relative angle to the n-s solar axis

;PRINT,'ROTATION ANGLE FOR THE LOCKHEED HELIOSTAT',ra

command = '/auto/home0/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(dom),1)+'.'+STRTRIM(STRING(heure),1)+'.'+STRTRIM(STRING(minute),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -17 | head -1 | cut -b19-32'
SPAWN,command,res
pangle = DOUBLE(res)

PRINT,'P-ANGLE',pangle ; returns the angle between the N-S solar axis and the solar rotation axis

;PRINT,'HELIOSTAT CONSTANT =',(90.d0+angle[a[0]])-pangle+ra

d=ROT(lam0g,angle[a[0]])
device,file='yo2.ps',bits=24,xsize=34,ysize=29,xoffset=-8,yoffset=-0.5,/color
tvim,d[255:511,*],rmarg=-4.5,range=[-450,1650] & tvim,dopplermdi[255:511,*],rmarg=4,range=[-450,1650] & tvim,d[0:255,*],rmarg=-4.5,range=[-2300,600] & tvim,dopplermdi[0:255,*],rmarg=4,range=[-2300,600] & device,/close

; NB: THE CONVENTION IS POSITIVE ANGLE FOR CLOCKWISE FOR ANGLE[A[0]]

SET_PLOT,'x'
!p.multi=0
READ,pause


END
