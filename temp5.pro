FSR             = DBLARR(7)             ; FSR in Angstrom
FSR[0]          = 0.172457d0            ; for the narrow-band Michelson
FSR[1]          = 0.344242d0            ; for the broad-band  Michelson
FSR[2]          = 0.7039;0.702d0               ; for E1
FSR[3]          = 1.405d0               ; for E2
FSR[4]          = 2.779d0               ; for E3
FSR[5]          = 5.682d0               ; for E4
FSR[6]          = 11.354d0              ; for E5


nseq            = 20
tuningJS        = DBLARR(nseq,3)
tuningJS[ 0,*]  = [ 53,72,44]
tuningJS[ 1,*]  = [ 41,66,47]
tuningJS[ 2,*]  = [ 29,60,50]
tuningJS[ 3,*]  = [ 77,54,53]
tuningJS[ 4,*]  = [ 65,48,56]
tuningJS[ 5,*]  = [ 53,42,59]
tuningJS[ 6,*]  = [ 41,96,62]
tuningJS[ 7,*]  = [ 29,90,65]
tuningJS[ 8,*]  = [ 77,84,68]
tuningJS[ 9,*]  = [ 65,78,71]
tuningJS[10,*]  = [ 53,72,74]
tuningJS[11,*]  = [ 41,66,77]
tuningJS[12,*]  = [ 29,60,80]
tuningJS[13,*]  = [ 77,54,83]
tuningJS[14,*]  = [ 65,48,86]
tuningJS[15,*]  = [ 53,42,89]
tuningJS[16,*]  = [ 41,96,92]
tuningJS[17,*]  = [ 29,90,95]
tuningJS[18,*]  = [ 77,84,98]
tuningJS[19,*]  = [ 65,78,101]
FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]


; PHASE MAPS
phase    = DBLARR(3)
phase[0] =  56.;-55.75    ; NB Michelson
phase[1] = -58.;130.80    ; WB Michelson
phase[2] =  62.;106.14    ; E1 
phase    = phase*!dpi/180.d0

dpi     = 2.d0*!dpi
nlam    = 2250                ; number of wavelengths
dlam    = 3.6d0/1.d3          ; resolution in wavelength
lam     =(DINDGEN(nlam)-(nlam-1)/2.)*dlam
lam0    = 0.0

; GAUSSIAN LINE PROFILE
;fwidthg = 0.07
;fdepthg = 0.51
;Icg     = 1
;templine= EXP(-(lam-lam0)^2.d0/fwidthg^2.d0)
;line    = Icg -fdepthg*templine
inten   = fltarr(nseq)

; ATLAS PROFILE: OBSERVED SOLAR LINE
RESTORE,'atlas_profile.sav'
lb1=lb1-6175.025  ; vacuum solar line
lb1=[lam[0],lb1,lam[nlam-1]]
pr1=[1.0,pr1,1.0]
line=INTERPOL(pr1,lb1,lam)

RESTORE,'frontwindow.bin' 
lamref       = 6173.3433
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam)
q            = READFITS('blocker11.fits')
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam)
profilef     = blocker0
contrasts    = FLTARR(7)+1.0;[0.9,0.98,0.95,0.9,0.94,0.95,0.96]
phases       = FLTARR(7);[20.,-12.,14.,-17.]*!dpi/180.d0

FOR i=0,3 DO profilef  = profilef * 0.5d0 * (1.d0+contrasts[i+2]*COS(dpi/FSR[i+3]*lam+phases[i]))

profile  = DBLARR(nseq,nlam)
FOR j=0,nseq-1 DO BEGIN
    profile[j,*] = profilef
    FOR i=0,2 DO profile[j,*] = profile[j,*] * 0.5d0 * (1.d0+contrasts[i]*COS(dpi/FSR[i]*lam+phase[i]+tuningJS[j,i]))
    inten[j] = TOTAL(profile[j,*]*line)*dlam
ENDFOR

;window,0,retain=2
set_plot,'ps'
device,file='yo.ps'
plot,inten,psym=4
oplot,[10,10],[0,1]
restore,'temp.bin'
oplot,inten0,linestyle=2
device,/close

END
