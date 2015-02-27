; THIS PROGRAM SHOWS THE LOCATION OF THE HMI FILTER
; AS A FUNCTION OF THE TUNABLE-ELEMENT PHASE MAPS

PRO cotune

SET_PLOT,'x'
!p.multi = 0

dpi         = 2.d0*!dpi
lamref      = 6173.3433        ; Fe I 6173 line central wavelength in air
nseq        = 6.;20               ; number of positions in the co-tune sequence
nlam        = 60000            ; number of wavelengths
dlam        = 1.d0/4.d3       ; resolution in wavelength
lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam

; CO-TUNE SEQUENCE FOR THE SUNLIGHT TESTS
tuningJS     = DBLARR(nseq,3)
;tuningJS[0,*]=[  111 , 100 ,   67]
;tuningJS[1,*]=[  117 ,  88 ,   43]
;tuningJS[2,*]=[  123 ,  76 ,   79]
;tuningJS[3,*]=[  129 , 124 ,   55]
;tuningJS[4,*]=[  135 , 112 ,   91]
;tuningJS[5,*]=[   81 , 100 ,   67]
;tuningJS[6,*]=[   87 ,  88 ,   43]
;tuningJS[7,*]=[   93,   76 ,   79]
;tuningJS[8,*]=[   99 , 124 ,   55]
;tuningJS[9,*]=[  105 , 112 ,   91]


; FOR LINE SCANS
;tuningJS[0,*] =[       126 ,         70   ,      67]
;tuningJS[1,*] =[       129 ,        124   ,      55]
;tuningJS[2,*] =[       132 ,        118   ,      43]
;tuningJS[3,*] =[       135 ,        112   ,      91]
;tuningJS[4,*] =[       138 ,        106   ,      79]
;tuningJS[5,*] =[        81 ,        100   ,      67]
;tuningJS[6,*] =[        84 ,         94   ,      55]
;tuningJS[7,*] =[        87 ,         88   ,      43]
;tuningJS[8,*] =[        90 ,         82   ,      91]
;tuningJS[9,*] =[        93 ,         76   ,      79]
;tuningJS[10,*]=[        96 ,         70   ,      67]
;
;tuningJS=REVERSE(tuningJS,2)
;FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]

; FOR DETUNE27
;tuningJS[ 0,*]  = [         0.d0,          0.d0 ,         0.d0]
;tuningJS[ 1,*]  = [        80.d0,          0.d0 ,         0.d0]
;tuningJS[ 2,*]  = [       160.d0,          0.d0 ,         0.d0]
;tuningJS[ 3,*]  = [         0.d0,         80.d0,          0.d0]
;tuningJS[ 4,*]  = [        80.d0,         80.d0,          0.d0]
;tuningJS[ 5,*]  = [       160.d0,         80.d0,          0.d0]
;tuningJS[ 6,*]  = [         0.d0,        160.d0,          0.d0]
;tuningJS[ 7,*]  = [        80.d0,        160.d0,          0.d0]
;tuningJS[ 8,*]  = [       160.d0,        160.d0,          0.d0]
;tuningJS[ 9,*]  = [         0.d0,          0.d0,         80.d0]
;tuningJS[10,*]  = [        80.d0,          0.d0,         80.d0]
;tuningJS[11,*]  = [       160.d0,          0.d0,         80.d0]
;tuningJS[12,*]  = [         0.d0,         80.d0,         80.d0] 
;tuningJS[13,*]  = [        80.d0,         80.d0,         80.d0]
;tuningJS[14,*]  = [       160.d0,         80.d0,         80.d0]
;tuningJS[15,*]  = [         0.d0,        160.d0,         80.d0]
;tuningJS[16,*]  = [        80.d0,        160.d0,         80.d0]
;tuningJS[17,*]  = [       160.d0,        160.d0,         80.d0]
;tuningJS[18,*]  = [         0.d0,          0.d0,        160.d0]
;tuningJS[19,*]  = [        80.d0,          0.d0,        160.d0]
;tuningJS[20,*]  = [       160.d0,          0.d0,        160.d0]
;tuningJS[21,*]  = [         0.d0,         80.d0,        160.d0]
;tuningJS[22,*]  = [        80.d0,         80.d0,        160.d0]
;tuningJS[23,*]  = [       160.d0,         80.d0,        160.d0]
;tuningJS[24,*]  = [         0.d0,        160.d0,        160.d0]
;tuningJS[25,*]  = [        80.d0,        160.d0,        160.d0]
;tuningJS[26,*]  = [       160.d0,        160.d0,        160.d0]
;
;FOR i=0,nseq-1 DO tuningJS[i,*]  = tuningJS[i,*] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0

; FOR COTUNE-22 (JUNE 2006)

;tuningJS[ 0,*]  = [ 53,72,44]
;tuningJS[ 1,*]  = [ 41,66,47]
;tuningJS[ 2,*]  = [ 29,60,50]
;tuningJS[ 3,*]  = [ 77,54,53]
;tuningJS[ 4,*]  = [ 65,48,56]
;tuningJS[ 5,*]  = [ 53,42,59]
;tuningJS[ 6,*]  = [ 41,96,62]
;tuningJS[ 7,*]  = [ 29,90,65]
;tuningJS[ 8,*]  = [ 77,84,68]
;tuningJS[ 9,*]  = [ 65,78,71]
;tuningJS[10,*]  = [ 53,72,74]
;tuningJS[11,*]  = [ 41,66,77]
;tuningJS[12,*]  = [ 29,60,80]
;tuningJS[13,*]  = [ 77,54,83]
;tuningJS[14,*]  = [ 65,48,86]
;tuningJS[15,*]  = [ 53,42,89]
;tuningJS[16,*]  = [ 41,96,92]
;tuningJS[17,*]  = [ 29,90,95]
;tuningJS[18,*]  = [ 77,84,98]
;tuningJS[19,*]  = [ 65,78,101]
;tuningJS[ 0,*]  = [82,29,22]
;tuningJS[ 1,*]  = [58,77,28]
;tuningJS[ 2,*]  = [94,65,34]
;tuningJS[ 3,*]  = [70,53,40]
;tuningJS[ 4,*]  = [106,41,46]
;tuningJS[ 5,*]  = [82,29,52]
tuningJS[ 0,*]  = [0,-30,-15]
tuningJS[ 1,*]  = [-24,18,-9]
tuningJS[ 2,*]  = [12,6,-3]
tuningJS[ 3,*]  = [-12,-6,3]
tuningJS[ 4,*]  = [24,-18,9]
tuningJS[ 5,*]  = [0,-30,15]

FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]


; PHASE MAPS
phase    = DBLARR(3)
phase[0] = 0.;-117.922;-53;43.;-55.75    ; NB Michelson
phase[1] = 0.;8.7036;-72;-71.;130.80    ; WB Michelson
phase[2] = 0.;-137.9398;74;83.;106.14    ; E1 
phase    = phase*!dpi/180.d0

; COMPUTATION OF THE FILTER PROFILES

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172457         ; for the narrow-band Michelson
FSR[1]      = 0.344242         ; for the broad-band  Michelson
FSR[2]      = 0.7039;0.702     ; for E1
FSR[3]      = 1.405            ; for E2
FSR[4]      = 2.779            ; for E3
FSR[5]      = 5.682            ; for E4
FSR[6]      = 11.354           ; for E5

RESTORE,'frontwindow.bin' 
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam)
q           = READFITS('blocker11.fits')
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lamref,lam)

;profilef = DBLARR(nlam)+1.d0 ; no blocker filter
profilef  = blocker0
FOR i=0,3 DO profilef  = profilef * 0.5d0 * (1.d0+COS(dpi/FSR[i+3]*lam))
profile  = DBLARR(nseq,nlam)

FOR j=0,nseq-1 DO BEGIN
    profile[j,*] = profilef
    FOR i=0,2 DO profile[j,*] = profile[j,*] * 0.5d0 * (1.d0+COS(dpi/FSR[i]*lam+phase[i]+tuningJS[j,i]))
ENDFOR

WINDOW,0,RETAIN=2,xsize=1200,ysize=900
!p.multi=0
;!p.multi=[0,3,2];[0,5,6]
;FOR i=0,nseq-1 DO 
plot,lam,profile[0,*],thick=3,xst=1,tit='!17',charsize=1.5,xrange=[-.6,.6],yrange=[0,.7]
oplot,lam,profile[1,*],linestyle=2
oplot,lam,profile[2,*],linestyle=3
oplot,lam,profile[3,*],linestyle=4
oplot,lam,profile[4,*],linestyle=2,thick=2
oplot,lam,profile[5,*],linestyle=3,thick=2,col=180
;oplot,lam,profile[6,*],linestyle=4,thick=2
;oplot,lam,profile[7,*],linestyle=2,thick=3
;oplot,lam,profile[8,*],linestyle=3,thick=3
;oplot,lam,profile[9,*],linestyle=4,thick=3
;oplot,lam,profile[10,*],linestyle=4,thick=3
;
;oplot,lam,blocker0


;SAVE,profile,FILE='richard_0.9_tun.bin'

maxi = DBLARR(nseq)
FOR i=0,nseq-1 DO BEGIN
    a = WHERE(profile[i,*] EQ MAX(profile[i,*]))
    PRINT,'POSITION NUMBER:',i,lam[a[0]]
ENDFOR

; FOR THE WAVELENGTH AND SPATIAL DEPENDENCE TEST: WITH 27 COTUNE
; SEQUENCES, AT WHICH POSITION SHOULD THE AMPLITUDES BE MAXIMUM ?

;nl0=27
;lam0    = DBLARR(nl0)
;lam0[0] = 6173.463d0
;lam0[1] = 6173.647d0
;lam0[2] = 6173.277d0
;lam0[3] = 6173.904d0
;lam0[4] = 6174.150d0
;lam0[5] = 6174.400d0
;lam0[6] = 6174.552d0
;lam0[7] = 6174.913d0
;lam0[8] = 6175.399d0
;lam0[9] = 6175.900d0
;lam0[10]= 6176.906d0
;lam0[11]= 6177.398d0
;lam0[12]= 6177.901d0
;lam0[13]= 6173.400d0
;lam0[14]= 6173.147d0
;lam0[15]= 6172.900d0
;lam0[16]= 6172.668d0
;lam0[17]= 6172.276d0
;lam0[18]= 6172.150d0
;lam0[19]= 6171.902d0
;lam0[20]= 6171.399d0
;lam0[21]= 6170.900d0
;lam0[22]= 6170.400d0
;lam0[23]= 6169.896d0
;lam0[24]= 6169.399d0
;lam0[25]= 6168.899d0
;lam0[26]= 6168.400d0
;lam0    = lam0-lamref
;
;maxi = LONARR(nl0)
;FOR i=0,nl0-1 DO BEGIN
;PRINT,i
;    maximum = 0
;    FOR j=0,nseq-1 DO BEGIN
;        temp = INTERPOL(profile[j,*],lam,lam0[i])
;        PRINT,temp
;        IF(temp GT maximum) THEN BEGIN
;            maxi[i]=j 
;            maximum = temp
;        ENDIF
;    ENDFOR
;READ,PAUSE
;ENDFOR
;FOR i=0,nl0-1 DO PRINT,lam0[i]+lamref,maxi[i]



READ,pause

END
