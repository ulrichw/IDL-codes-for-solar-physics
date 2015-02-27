; to reconstruct the spatial-average intensity
; variation across a detune sequence

PRO HMI_reconsaverage

nlam        = 1125                ; number of wavelengths
dlam        = 6.d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam
contrasts   = 0.97
fwidthg     = 0.068510673
fdepthg     = 0.463
dpi         = 2.d0*!dpi
; PARAMETERS FOR THE I-RIPPLE ON THE NB MICHELSON
amplitudeNB=0.058936515d0
phaseNB  = -102.72668d0*!pi/180.
amplitudeWB=0.0096328533
phaseWB  = 130.88739d0*!pi/180.

;list ='listSun060224_225806'
;list ='listSun060621_170559' ; in Calmode
;list ='listSun060703_205709' ; in Calmode
;list ='listSun060714_221512' ; wobble
 list='listSun060714_231213' ; wobble
;list='listLamp060615_170702' ; lamp data
;list = 'list_060304_031619'

nx   = 128
nseq = 36;27
RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;Inten = REFORM(REBIN(imx[*,*,1:27],1,1,nseq))
 Inten = REBIN(imx,1,1,nseq)
;lam0=0.0018532850 ; for listSun060714_221512
lam0=0.0032872231 ; for listSun060714_231213
;lam0  = lam0*6173.3433/3.d8
Inten = Inten/MEAN(Inten)

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172457d0            ; for the narrow-band Michelson
FSR[1]      = 0.344242d0            ; for the broad-band  Michelson
FSR[2]      = 0.7039;0.702d0               ; for E1
FSR[3]      = 1.405d0               ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5

GOTO,yop
tuning       = FLTARR(3,nseq)
tuning[*,0]  = [         0.d0,          0.d0 ,         0.d0]
tuning[*,1]  = [        80.d0,          0.d0 ,         0.d0]
tuning[*,2]  = [       160.d0,          0.d0 ,         0.d0]
tuning[*,3]  = [         0.d0,         80.d0,          0.d0]
tuning[*,4]  = [        80.d0,         80.d0,          0.d0]
tuning[*,5]  = [       160.d0,         80.d0,          0.d0]
tuning[*,6]  = [         0.d0,        160.d0,          0.d0]
tuning[*,7]  = [        80.d0,        160.d0,          0.d0]
tuning[*,8]  = [       160.d0,        160.d0,          0.d0]
tuning[*,9]  = [         0.d0,          0.d0,         80.d0]
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
yop:

; FOR WOBBLE SEQUENCE 1
tuning       = FLTARR(3,nseq)
tuning[*,0]  = [   0,   0,   0]
tuning[*,1]  = [  10,   0,   0]
tuning[*,2]  = [  20,   0,   0]
tuning[*,3]  = [  30,   0,   0]
tuning[*,4]  = [  40,   0,   0]
tuning[*,5]  = [  50,   0,   0]
tuning[*,6]  = [  60,   0,   0]
tuning[*,7]  = [  70,   0,   0]
tuning[*,8]  = [  80,   0,   0]
tuning[*,9]  = [  90,   0,   0]
tuning[*,10] = [ 100,   0,   0]
tuning[*,11] = [ 110,   0,   0]
tuning[*,12] = [  0,   0,   0]
tuning[*,13] = [  0,  10,   0]
tuning[*,14] = [  0,  20,   0]
tuning[*,15] = [  0,  30,   0]
tuning[*,16] = [  0,  40,   0]
tuning[*,17] = [  0,  50,   0]
tuning[*,18] = [  0,  60,   0]
tuning[*,19] = [  0,  70,   0]
tuning[*,20] = [  0,  80,   0]
tuning[*,21] = [  0,  90,   0]
tuning[*,22] = [  0, 100,   0]
tuning[*,23] = [  0, 110,   0]
tuning[*,24] = [  0,   0,   0]
tuning[*,25] = [  0,   0,  10]
tuning[*,26] = [  0,   0,  20]
tuning[*,27] = [  0,   0,  30]
tuning[*,28] = [  0,   0,  40]
tuning[*,29] = [  0,   0,  50]
tuning[*,30] = [  0,   0,  60]
tuning[*,31] = [  0,   0,  70]
tuning[*,32] = [  0,   0,  80]
tuning[*,33] = [  0,   0,  90]
tuning[*,34] = [  0,   0, 100]
tuning[*,35] = [  0,   0, 110]

tuning       = REVERSE(tuning,1)
FOR i=0,nseq-1 DO tuning[*,i]=REFORM(tuning[*,i])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]



profilef = FLTARR(nlam)+1.0
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(2.*!pi/FSR[i]*lam))/2.d0

; THEORETICAL GAUSSIAN PROFILE
templine= EXP(-(lam-lam0)^2.d0/fwidthg^2.d0)
line    = 1.0 -fdepthg*templine

; ATLAS PROFILE: OBSERVED SOLAR LINE
;RESTORE,'atlas_profile.sav'
;lb1=lb1-6175.025  ; vacuum solar line
;lb1=[lam[0],lb1,lam[nlam-1]]
;pr1=[1.0,pr1,1.0]
;line=INTERPOL(pr1,lb1,lam)


ntest   = 47.
Phig    = FLTARR(3,ntest)
Phig[0,*] = (FINDGEN(ntest)/ntest*280.-140)*!pi/180
Phig[1,*] = (FINDGEN(ntest)/ntest*280.-140)*!pi/180
Phig[2,*] = (FINDGEN(ntest)/ntest*280.-140)*!pi/180

fit     = FLTARR(nseq)

residu  = FLTARR(ntest,ntest,ntest)


mini = 10000.
FOR a=0,ntest-1 DO BEGIN
PRINT,a
FOR b=0,ntest-1 DO BEGIN
FOR c=0,ntest-1 DO BEGIN

FOR j=0,nseq-1 DO BEGIN
    profileg      = 0.125d0 * profilef * (1.d0+COS(2.*!pi/FSR[0]*lam+Phig[0,a]+tuning[0,j])) * (1.d0+COS(2.*!pi/FSR[1]*lam+Phig[1,b]+tuning[1,j])) * (1.d0+COS(2.*!pi/FSR[2]*lam+Phig[2,c]+tuning[2,j]))
   ;fit[j] = TOTAL(line*profileg)*dlam ; solar line
    fit[j] = TOTAL(profileg)*dlam*(1.d0+amplitudeNB*COS(tuning[0,j]+phaseNB))*(1.d0+amplitudeWB*COS(tuning[1,j]+phaseWB)) ; lamp data
ENDFOR
fit = fit/MEAN(fit)
residu[a,b,c]=TOTAL((fit-inten)^2.0)
if residu[a,b,c] lt mini then begin
    mini = residu[a,b,c]
    fitmini = fit
endif
ENDFOR
ENDFOR
ENDFOR

!P.MULTI=0
SET_PLOT,'ps'
;WINDOW,0,retain=2
device,file='yo.ps',xsize=15,ysize=12,xoffset=0,yoffset=0
plot,inten,xst=1
oplot,fitmini,linestyle=2
device,/close
set_plot,'x'

READ,pause

END
