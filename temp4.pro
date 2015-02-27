dpi = 2.d0*!dpi
nseq= 27
lamref=6173.3433d0

depth=0.5165;0.737857;0.5165
width=0.0843;0.10
contrasts=0.85;0.98
phases=20.*!dpi/180.d0
continuum=75000.

lam0=0.0
nlam        = 2250                ; number of wavelengths
dlam        = 3.6d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam 
tuning      = fltarr(3,nseq)

line=continuum-depth*continuum*exp(-(lam-lam0)^2.0/width^2.)

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172457d0            ; for the narrow-band Michelson
FSR[1]      = 0.344242d0            ; for the broad-band  Michelson
FSR[2]      = 0.7039;0.702d0               ; for E1
FSR[3]      = 1.405d0               ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5

FSR         = dpi/FSR

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


RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

q            = READFITS('blocker11.fits') ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
                                          ; this profile, obtained with a Cary, is centered on 6169.8 A, but measures made by Andover
                                          ; with a spectrograph show the center on 6172 A instead of 6173.3433. Also
                                          ; we shift the profile by 1.3433 A

blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam);,/LSQUADRATIC)

profilef= blocker0  ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam+phases))/2.d0
 
In = FLTARR(nseq)
FOR j=0,nseq-1 DO BEGIN 
    profileg = 0.125d0 * profilef *(1.d0+contrasts*cos(FSR[0]*lam+phases+tuning[0,j]))* (1.d0+contrasts*cos(FSR[1]*lam+phases+tuning[1,j])) * (1.d0+contrasts*cos(FSR[2]*lam+phases+tuning[2,j])) 
    in[j]= TOTAL(line*profileg)*dlam 
endfor

;GOTO,fuck
seed=1l
contrast = RANDOMU(seed,64,64)*0.15+0.8
phase    = RANDOMU(seed,64,64,3)*dpi-!dpi

Inten=fltarr(64,64,27)
FOR i=0,63 do begin
    PRINT,i
    for j=0,63 do begin
        for k=0,26 do begin
            profileg = 0.125d0 * profilef *(1.d0+contrast[i,j]*cos(FSR[0]*lam+phase[i,j,0]+tuning[0,k]))* (1.d0+contrast[i,j]*cos(FSR[1]*lam+phase[i,j,1]+tuning[1,k])) * (1.d0+contrast[i,j]*cos(FSR[2]*lam+phase[i,j,2]+tuning[2,k])) 
            inten[i,j,k]= TOTAL(line*profileg)*dlam 
            inten[i,j,k]= inten[i,j,k] + SQRT(inten[i,j,k])*RANDOMN(seed)
        endfor
    endfor
endfor

save,inten,phase,contrast,file='temp.bin'

read,pause
fuck:


mat      = DBLARR(3,2,2)
vec      = DBLARR(1,2)
FOR i=0,2 DO BEGIN
    mat[i,0,0] =  COS(dpi/FSR[i]*lam0)
    mat[i,1,0] = -SIN(dpi/FSR[i]*lam0)
    mat[i,0,1] =  COS(dpi/FSR[i]*lam0 + dpi/3.d0)
    mat[i,1,1] = -SIN(dpi/FSR[i]*lam0 + dpi/3.d0)
    mat[i,*,*] =  LA_INVERT(REFORM(mat[i,*,*]),/DOUBLE)
ENDFOR

FOR i=0,2 DO BEGIN
    IF(i EQ 0) THEN BEGIN
        I00= (In[0]+In[3]+In[6]+In[9]+In[12]+In[15]+In[18]+In[21]+In[24])
        I1 = (In[1]+In[4]+In[7]+In[10]+In[13]+In[16]+In[19]+In[22]+In[25])
        I2 = (In[2]+In[5]+In[8]+In[11]+In[14]+In[17]+In[20]+In[23]+In[26])
    ENDIF
    IF(i EQ 1) THEN BEGIN
        I00= (In[0]+In[1]+In[2]+In[9]+In[10]+In[11]+In[18]+In[19]+In[20])
        I1 = (In[3]+In[4]+In[5]+In[12]+In[13]+In[14]+In[21]+In[22]+In[23])
        I2 = (In[6]+In[7]+In[8]+In[15]+In[16]+In[17]+In[24]+In[25]+In[26])
    ENDIF
    IF(i EQ 2) THEN BEGIN
        I00= (In[0]+In[1]+In[2]+In[3]+In[4]+In[5]+In[6]+In[7]+In[8])
        I2 = (In[9]+In[10]+In[11]+In[12]+In[13]+In[14]+In[15]+In[16]+In[17])
        I1 = (In[18]+In[19]+In[20]+In[21]+In[22]+In[23]+In[24]+In[25]+In[26])
    ENDIF

    tot      = I00+I1+I2
    vec[0,0] = I00*3.d0/tot-1.d0
    vec[0,1] = I1 *3.d0/tot-1.d0
    vec      = REFORM(mat[i,*,*])##vec
    print,SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)
    print,ATAN(vec[0,1],vec[0,0])
ENDFOR    

READ,pause




END
