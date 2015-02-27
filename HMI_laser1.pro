; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE DE-TUNE SEQUENCE
; OBJECTIVE: TO LEARN HOW TO CO-TUNE THE ELEMENTS
; AND THE PHASES AND CONTRASTS OF THE TUNABLE ELEMENTS
;
; THIS PROGRAMS MENTIONS THE PRESENCE OF BAD EXPOSURES AND OVERSCANS
;
; ver 1.3 June 12, 2006
;
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser1,draw;,list,lam0

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0 ; solar Fe I 6173 line central wavelength in air
; when phases are at 0, we are centered on this wavelength

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 128;256              ; number of rows
ny          = 128;256              ; number of columns

Bg0         = DBLARR(nx,ny,3)  ; contrasts of the tunable elements
Phig0       = DBLARR(nx,ny,3)  ; phases of the tunable elements
nseq        = 27               ; number of positions in the de-tune sequence
Inten       = DBLARR(nx,ny,nseq);measured output intensities (measured on a HMI CCD)
tuning      = DBLARR(3,nseq)   ; tuning positions for the detune sequence
anglim      = 930.;960.

xcenter     = 66;nx/2; depends on the image !
ycenter     = 62;nx/2;

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;--------------------------------------------

FSR         = DBLARR(7)        ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.170942d0;0.172d0;0.1709709199186409401d0;0.1720025438132122275d0;0.172457d0 ; for the narrow-band Michelson
FSR[1]      = 0.341923d0;0.342d0;0.3423194217700681330d0;0.3419945321588832021d0;;0.344242 ; for the broad-band  Michelson
FSR[2]      = 0.693483d0;0.693d0;0.6930111472833409003d0;lamref/2.d0/!dpi/1424.d0;0.690d0;0.7039;0.702     ; for E1
FSR[3]      = 1.407d0;1.404948406918525405d0            ; for E2
FSR[4]      = 2.779d0;2.779533228275551604d0            ; for E3
FSR[5]      = 5.682d0;5.684478176795580318d0            ; for E4
FSR[6]      = 11.354d0;11.34805753676470630d0            ; for E5

; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
;                NB MICHELSON   WB MICHELSON         E1 LYOT
;-----------------------------------------------------------


tuning[*,0]  = [         0.d0,          0.d0,          0.d0]
tuning[*,1]  = [        80.d0,          0.d0,          0.d0]
tuning[*,2]  = [       160.d0,          0.d0,          0.d0]
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


; INFORMATION TO BE ADDED BY THE USER
;----------------------------------------

;list         = 'listLaser060131_014331' ; in Calmode
;list         = 'listLaser060210_230522' ; in Calmode
;list         = 'listLaser060210_225517' ; in Obsmode
;list         = 'listLaser060210_221457' ; in Calmode
;list         = 'listLaser060303_181526' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_182759' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_184032' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_185253' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_190514' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_191737' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060620_232532' ; in Calmode
;list         = 'listLaser060620_000903' ; in Calmode LOT OF PROBLEMS
;list         = 'listLaser060620_234618' ; in Calmode
;list         = 'listLaser060622_195601' ; in Calmode
;list         = 'listLaser060622_182748' ; in Obsmode
;list         = 'listLaser060622_194336' ; in Obsmode
list         = 'listLaser060622_201041' ; in Obsmode with diffuser
;list         = 'listLaser060622_202816' ; in Obsmode with diffuser
;list         = 'listLaser060622_221251' ; in Obsmode 
;list         = 'listLaser060622_214137' ; in Obsmode with diffuser
;list         = 'listLaser060622_223947' ; in Obsmode
;list         = 'listLaser060622_233055' ; in Obsmode
;list         = 'listLaser060720_205919' ; in Obsmode
;list         = 'listLaser060720_173252' ; in Obsmode
;list         = 'listLaser060720_170649' ; in Obsmode
;list         = 'listLaser060721_184158' ; in Obsmode
;list ='listSun060616_211450'
;list ='listSun060224_225806'
;list ='listSun060224_224805'
;list='listSun060623_230811'
;list='listLaser070108_40811' ; in Obsmode ; problem: random dots target forgotten USED datatreat.pro
;list='listLaser070108_40843'
;list='listLaser070108_40875'
;list='listLaser070108_40907'
;list='listLaser070108_40939'
;list ='listLaser070220_74462' ; OBSMODE in 512*512 and 256*256
;list ='listLaser070220_74492'
;list='listLaser070220_74548' ; OBSMODE
;list='listLaser070220_74578'
;list='listLaser070220_74998' ; OBSMODE
;list='listLaser070220_75028'
;list='listLaser070223_79336'
;list='listLaser070223_79886' ; OBSMODE
;list='listLaser070222_79306'
;list='listLaser070222_79916'

;list='listLaser070225_85408'
;list='listLaser070225_85439'
;list='listLaser070225_85470'
;list='listLaser070225_85532'
;list='listLaser070225_85593'
;list='listLaser070225_85624'
;list='listLaser070225_85655'
;list='listLaser070225_85686'
;list='listLaser070225_85500'

;list='listLaser070225_86978'
;list='listLaser070225_87009'
;list='listLaser070225_87040'
;list='listLaser070225_87071'
;list='listLaser070225_87102'
;;list='listLaser070225_87133' ; fisrt darks missing, and then shitty data...
;list='listLaser070225_87164'
;list='listLaser070225_87195'
;list='listLaser070225_87226'

;list='listLaser070225_87342'
;list='listLaser070225_87373'
;list='listLaser070225_87404'
;list='listLaser070225_87435'
;list='listLaser070225_87466'
;list='listLaser070225_87497'
;list='listLaser070225_87528'
;list='listLaser070225_87559'
;list='listLaser070225_87590'
;list='listLaser070225_85780' ; OBSMODE 1 image is missing (29/30)
;list='listLaser070225_85814'
;list='listLaser070225_85882'
;list='listLaser070225_85916'
;list='listLaser070225_85950'
;list='listLaser070225_85984'
;list='listLaser070225_86018'
;list='listLaser070225_86052'

;list='listLaser070228_98858'
;list='listLaser070228_98888'
;list='listLaser070228_98918'
;list='listLaser070228_98948'
;list='listLaser070228_98978'
;list='listLaser070228_99008'
;list='listLaser070228_99038'
;list='listLaser070228_99068'
;list='listLaser070228_99098'
;list='listLaser070228_99128'
;list='listLaser070228_99158'
;list='listLaser070228_99188'
;list='listLaser070228_99218'
;list='listLaser070228_99248'
;list='listLaser070228_99278'

;list='listLaser070227_93773'
;list='listLaser070227_93803'
;list='listLaser070227_93833' ; one image missing (10/30), small field stop
;list='listLaser070228_98035'
;list='listLaser070228_98156' ; one image missing (8/30), small field stop
;list='listLaser070228_98065'

;list='listLaser070227_94221' ; OBSMODE, thermal test with laser (from 28 to 35 degrees C)
;list='listLaser070227_94251'
;list='listLaser070227_94281'
;list='listLaser070227_94311'
;list='listLaser070227_94341'
;list='listLaser070227_94371'
;list='listLaser070227_94401'
;list='listLaser070227_94431'
;list='listLaser070227_94461'
;list='listLaser070227_94491'
;list='listLaser070227_94521'
;list='listLaser070227_94551'

;list='listLaser070225_86610' ; OBSMODE
;list='listLaser070225_85746' ; OBSMODE
;list='listLaser070225_86436' ; OBSMODE
;list='listLaser070225_87624' ; OBSMODE
;list='listLaser070222_78530' ; CALMODE
;list='listLaser070222_78560' ; OBSMODE
;list='listLaser070225_87311' ; CALMODE

;list='listSun070227_89915' ; SUN; CALMODE; TO CHECK THE PHASES

;list='listLaser070227_90276' ; OBSMODE
;list='listLaser070227_90726' ; OBSMODE

;wgood        = 0.0
;READIMAGES,list,images,headers,time,29,wgood,nbin=nx,iover=iover
;CLEAN,images,imx,headers
;SAVE,imx,time,wgood,FILE='SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;Inten = DOUBLE(imx[*,*,2:28]) ; after February 2007
Inten = imx[*,*,1:27]
;Inten = images[*,*,1:27]
;wgood=FINDGEN(29) ; meaning everything is good !

; WE GET RID OF THE OVERSCANS AND BAD EXPOSURES

FOR i=0,nseq-1 DO BEGIN
    a=WHERE( wgood EQ (i+1))
    IF(a[0] EQ -1) THEN BEGIN
        Inten[*,*,i] = -1.0
        PRINT,'SOME BAD IMAGES !!!!!!'
    ENDIF
ENDFOR

;Inten[*,*,26]=-1.0 ;for 'listLaser070225_85780'
;Inten[*,*,5]=-1.0 ;for 'listLaser070228_98156'
;Inten[*,*,7]=-1.0 ;for 'listLaser070227_93833'

; IF BAD DARK FRAME REMOVAL, INTEN HAS NEGATIVE VALUES
; NO NEED TO REMOVE APPARENTLY
;a  = WHERE(Inten lt 0.0)
;Inten = Inten + MIN(Inten)
;IF(a[0] NE -1) THEN Inten[a] = 0.0


; LASER CENTRAL WAVELENGTH

;GOTO,jose
CASE list OF
    'listLaser060131_014331': lam0          = 6173.538d0-lamref 
    'listLaser060210_221457': lam0          = 6173.585d0-lamref 
    'listLaser060210_225517': lam0          = 6173.463d0-lamref 
    'listLaser060210_230522': lam0          = 6173.463d0-lamref
    'listLaser060303_181526': lam0          = 6173.467d0-lamref 
    'listLaser060303_182759': lam0          = 6173.467d0-lamref
    'listLaser060303_184032': lam0          = 6173.467d0-lamref
    'listLaser060303_185253': lam0          = 6173.467d0-lamref
    'listLaser060303_190514': lam0          = 6173.467d0-lamref
    'listLaser060303_191737': lam0          = 6173.467d0-lamref
    'listLaser060620_232532': lam0          = 6173.227d0-lamref
    'listLaser060620_000903': lam0          = 6173.227d0-lamref
    'listLaser060620_234618': lam0          = 6173.227d0-lamref
    'listLaser060622_195601': lam0          = 6173.397d0-lamref
    'listLaser060622_194336': lam0          = 6173.397d0-lamref
    'listLaser060622_182748': lam0          = 6173.450d0-lamref
    'listLaser060622_201041': lam0          = 6173.397d0-lamref
    'listLaser060622_202816': lam0          = 6173.243d0-lamref
    'listLaser060622_221251': lam0          = 6173.512d0-lamref
    'listLaser060622_214137': lam0          = 6173.512d0-lamref
    'listLaser060622_223947': lam0          = 6173.640d0-lamref
    'listLaser060622_233055': lam0          = 6173.080d0-lamref
    'listLaser060720_205919': lam0          = 6328.0d0  -lamref ; Lyot had been removed
    'listLaser060720_173252': lam0          = 6328.0d0  -lamref ; Lyot had been removed
    'listLaser060720_170649': lam0          = 6328.0d0  -lamref ; Lyot had been removed
    'listLaser060721_184158': lam0          = 6120.0d0  -lamref
    'listSun060616_211450'  : lam0          = 0.005149
    'listSun060224_225806'  : lam0          = 0.013049383
    'listSun060224_224805'  : lam0          = 0.012732526
    'listSun060623_230811'  : lam0          = 0.0067533583
    'listLaser070108_40811' : lam0          = 6173.325d0-lamref
    'listLaser070108_40843' : lam0          = 6173.484d0-lamref
    'listLaser070108_40875' : lam0          = 6173.238d0-lamref
    'listLaser070108_40907' : lam0          = 6173.111d0-lamref
    'listLaser070108_40939' : lam0          = 6173.600d0-lamref
    'listLaser070220_74462' : lam0          = 6173.3255d0-lamref;6173.330d0-lamref
    'listLaser070220_74492' : lam0          = 6173.3265d0-lamref
    'listLaser070220_74548' : lam0          = 6173.3285d0-lamref;6173.330d0-lamref
    'listLaser070220_74578' : lam0          = 6173.3295d0-lamref
    'listLaser070220_74998' : lam0          = 6173.400d0-lamref
    'listLaser070220_75028' : lam0          = 6173.257d0-lamref;6173.400d0-lamref
    'listLaser070223_79336' : lam0          = 6173.340d0-lamref
    'listLaser070223_79886' : lam0          = 6173.340d0-lamref
    'listLaser070222_79306' : lam0          = 6173.340d0-lamref
    'listLaser070222_79916' : lam0          = 6173.340d0-lamref
    'listLaser070225_85408' : lam0          = 6173.340d0-lamref
    'listLaser070225_85439' : lam0          = 6173.341d0-lamref
    'listLaser070225_85470' : lam0          = 6173.341d0-lamref
    'listLaser070225_85532' : lam0          = 6173.3405d0-lamref
    'listLaser070225_85593' : lam0          = 6173.340d0-lamref
    'listLaser070225_85624' : lam0          = 6173.3395d0-lamref
    'listLaser070225_85655' : lam0          = 6173.339d0-lamref
    'listLaser070225_85686' : lam0          = 6173.3395d0-lamref
    'listLaser070225_85500' : lam0          = 6173.341d0-lamref
    'listLaser070225_86978' : lam0          = 6172.327d0-lamref
    'listLaser070225_87009' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87040' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87071' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87102' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87164' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87195' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87226' : lam0          = 6172.331d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87342' : lam0          = 6174.327d0-lamref
    'listLaser070225_87373' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87404' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87435' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87466' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87497' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87528' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87559' : lam0          = 6174.326d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_87590' : lam0          = 6174.325d0-lamref ; WARNING: CHANGE VALUE !!!!
    'listLaser070225_85780' : lam0          = 6173.348d0-lamref ; WAVELENGTH DRIFTS FROM .344 TO .348 !
    'listLaser070225_85814' : lam0          = 6173.349d0-lamref
    'listLaser070225_85882' : lam0          = 6173.349d0-lamref
    'listLaser070225_85916' : lam0          = 6173.349d0-lamref
    'listLaser070225_85950' : lam0          = 6173.349d0-lamref
    'listLaser070225_85984' : lam0          = 6173.349d0-lamref
    'listLaser070225_86018' : lam0          = 6173.348d0-lamref
    'listLaser070225_86052' : lam0          = 6173.348d0-lamref
    'listLaser070228_98858' : lam0          = 6173.452d0-lamref
    'listLaser070228_98888' : lam0          = 6173.452d0-lamref
    'listLaser070228_98918' : lam0          = 6173.452d0-lamref
    'listLaser070228_98948' : lam0          = 6173.452d0-lamref
    'listLaser070228_98978' : lam0          = 6173.452d0-lamref
    'listLaser070228_99008' : lam0          = 6173.452d0-lamref
    'listLaser070228_99038' : lam0          = 6173.452d0-lamref
    'listLaser070228_99068' : lam0          = 6173.452d0-lamref
    'listLaser070228_99098' : lam0          = 6173.452d0-lamref
    'listLaser070228_99128' : lam0          = 6173.452d0-lamref
    'listLaser070228_99158' : lam0          = 6173.452d0-lamref
    'listLaser070228_99188' : lam0          = 6173.452d0-lamref
    'listLaser070228_99218' : lam0          = 6173.452d0-lamref
    'listLaser070228_99248' : lam0          = 6173.452d0-lamref
    'listLaser070228_99278' : lam0          = 6173.452d0-lamref
    'listLaser070227_93773' : lam0          = 6173.452d0-lamref
    'listLaser070227_93803' : lam0          = 6173.452d0-lamref
    'listLaser070227_93833' : lam0          = 6173.424d0-lamref
    'listLaser070228_98035' : lam0          = 6173.452d0-lamref
    'listLaser070228_98156' : lam0          = 6173.452d0-lamref
    'listLaser070228_98065' : lam0          = 6173.452d0-lamref
    'listLaser070227_94221' : lam0          = 6173.173d0-lamref ; estimate from log, see wavemeter file to confirm
    'listLaser070227_94251' : lam0          = 6173.173d0-lamref
    'listLaser070227_94281' : lam0          = 6173.173d0-lamref
    'listLaser070227_94311' : lam0          = 6173.173d0-lamref
    'listLaser070227_94341' : lam0          = 6173.173d0-lamref
    'listLaser070227_94371' : lam0          = 6173.173d0-lamref
    'listLaser070227_94401' : lam0          = 6173.173d0-lamref
    'listLaser070227_94431' : lam0          = 6173.173d0-lamref
    'listLaser070227_94461' : lam0          = 6173.173d0-lamref
    'listLaser070227_94491' : lam0          = 6173.173d0-lamref
    'listLaser070227_94521' : lam0          = 6173.173d0-lamref
    'listLaser070227_94551' : lam0          = 6173.173d0-lamref
    'listLaser070225_86610' : lam0          = 6172.331d0-lamref ; wavemeter file
    'listLaser070225_85746' : lam0          = 6173.344d0-lamref ; wavemeter file ; OBSMODE
    'listLaser070225_86436' : lam0          = 6172.3315d0-lamref ; wavemeter file ; OBSMODE
    'listLaser070225_87624' : lam0          = 6174.325d0-lamref ; wavemeter file ; OBSMODE
    'listLaser070222_78530' : lam0          = 6173.340d0-lamref
    'listLaser070222_78560' : lam0          = 6173.340d0-lamref
    'listLaser070225_87311' : lam0          = 6174.3275d0-lamref
    'listSun070227_89915' : lam0 = 211.24707*6173.3433d0/2.99792458d8
    'listLaser070227_90276' : lam0          = 6173.347d0-lamref
    'listLaser070227_90726' : lam0          = 6173.347d0-lamref
ENDCASE
;jose:

;RESTORE,'temp.bin'
;Inten0=fltarr(1,1,27)
;Inten0[0,0,*]=inten
;inten=inten0
;lam0=0.0617d0
;nx=1
;ny=1
;anglim=2048.
;xcenter=0
;ycenter=0
;distance=0.
;Bg0         = DBLARR(1,1,3)  ; contrasts of the tunable elements
;Phig0       = DBLARR(1,1,3)  ; phases of the tunable elements

distance=FLTARR(nx,nx)
for i=0,nx-1 do for j=0,nx-1 do distance[i,j]=SQRT((i-xcenter)^2.d0 + (j-ycenter)^2.d0)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center

IF(draw EQ 1) THEN GOTO,draw

; BEGINNING OF THE ITERATIONS
;----------------------------------------

mat      = DBLARR(3,2,2)
vec      = DBLARR(1,2)
FOR i=0,2 DO BEGIN
    mat[i,0,0] =  COS(dpi/FSR[i]*lam0)
    mat[i,1,0] = -SIN(dpi/FSR[i]*lam0)
    mat[i,0,1] =  COS(dpi/FSR[i]*lam0 + dpi/3.d0)
    mat[i,1,1] = -SIN(dpi/FSR[i]*lam0 + dpi/3.d0)
    mat[i,*,*] =  LA_INVERT(REFORM(mat[i,*,*]),/DOUBLE)
ENDFOR


; FOR RECONSTRUCTION OF THE DETUNE SEQUENCE
Intenrecons = DBLARR(nx,ny,nseq)
RESTORE,'frontwindow.bin'
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam0)
q           = READFITS('blocker11.fits')
;blocker0   = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam0);,/LSQUADRATIC)
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+2.6-lamref,lam0);,/LSQUADRATIC)


; BEGINNING OF COMPUTATIONS
In = DBLARR(nseq) 
FOR jjj=0,ny-1 DO BEGIN
    FOR iii=0,nx-1 DO BEGIN
  ; GUESS OF THE PHASES AND CONTRASTS OF TUNABLE ELEMENTS 
  ; ASSUMING THE LASER LINE PROFILE IS A DELTA FUNCTION
  ; we use the fact that cos(a)+cos(a+2pi/3)+cos(a+4pi/3) = 0
  ; I0 = [ cos(2\pi l0/FSR0)         -sin(2\pi l0/FSR0)         ] * B cos(\phi) * T/2 + T/2
  ; I1   [ cos(2\pi l0/FSR0+2\pi/3)  -sin(2\pi l0/FSR0+2\pi/3)  ]   B sin(\phi) * T/2 + T/2
  ; with 3 T/2 = I0+I1+I2
  ; USING THE MEASURED INTENSITIES
  ; THE FOLLOWING PROCEDURE IS WRITTEN FOR [1.5,1.5,-1.5]
  ; (WHEN ONE OF THE COEFFICIENT IS -1.5, THEN I1=I2 AND I2=I1) 
  ;----------------------------------------------------------

        IF(distance[iii,jjj] LE anglim) THEN BEGIN
           FOR i=0,2 DO BEGIN
           CASE i OF
               0: BEGIN
                   In[*] = Inten[iii,jjj,*]
                 ; TO TAKE CARE OF THE DUBIOUS IMAGES
                   IF(In[0] EQ -1.0 OR In[1] EQ -1.0 OR In[2] EQ -1.0) THEN BEGIN
                       In[0] = 0.0
                       In[1] = 0.0
                       In[2] = 0.0
                   ENDIF
                   IF(In[3] EQ -1.0 OR In[4] EQ -1.0 OR In[5] EQ -1.0) THEN BEGIN
                       In[3] = 0.0
                       In[4] = 0.0
                       In[5] = 0.0
                   ENDIF
                   IF(In[6] EQ -1.0 OR In[7] EQ -1.0 OR In[8] EQ -1.0) THEN BEGIN
                       In[6] = 0.0
                       In[7] = 0.0
                       In[8] = 0.0
                   ENDIF
                   IF(In[9] EQ -1.0 OR In[10] EQ -1.0 OR In[11] EQ -1.0) THEN BEGIN
                       In[9] = 0.0
                       In[10] = 0.0
                       In[11] = 0.0
                   ENDIF
                   IF(In[12] EQ -1.0 OR In[13] EQ -1.0 OR In[14] EQ -1.0) THEN BEGIN
                       In[12] = 0.0
                       In[13] = 0.0
                       In[14] = 0.0
                   ENDIF
                   IF(In[15] EQ -1.0 OR In[16] EQ -1.0 OR In[17] EQ -1.0) THEN BEGIN
                       In[15] = 0.0
                       In[16] = 0.0
                       In[17] = 0.0
                   ENDIF
                   IF(In[18] EQ -1.0 OR In[19] EQ -1.0 OR In[20] EQ -1.0) THEN BEGIN
                       In[18] = 0.0
                       In[19] = 0.0
                       In[20] = 0.0
                   ENDIF
                   IF(In[21] EQ -1.0 OR In[22] EQ -1.0 OR In[23] EQ -1.0) THEN BEGIN
                       In[21] = 0.0
                       In[22] = 0.0
                       In[23] = 0.0
                   ENDIF
                   IF(In[24] EQ -1.0 OR In[25] EQ -1.0 OR In[26] EQ -1.0) THEN BEGIN
                       In[24] = 0.0
                       In[25] = 0.0
                       In[26] = 0.0
                   ENDIF
                   I00= (In[0]+In[3]+In[6]+In[9]+In[12]+In[15]+In[18]+In[21]+In[24])
                   I1 = (In[1]+In[4]+In[7]+In[10]+In[13]+In[16]+In[19]+In[22]+In[25])
                   I2 = (In[2]+In[5]+In[8]+In[11]+In[14]+In[17]+In[20]+In[23]+In[26])
               END
               1: BEGIN
                   In[*] = Inten[iii,jjj,*]
                   IF(In[0] EQ -1.0 OR In[3] EQ -1.0 OR In[6] EQ -1.0) THEN BEGIN
                       In[0] = 0.0
                       In[3] = 0.0
                       In[6] = 0.0
                   ENDIF
                   IF(In[1] EQ -1.0 OR In[4] EQ -1.0 OR In[7] EQ -1.0) THEN BEGIN
                       In[1] = 0.0
                       In[4] = 0.0
                       In[7] = 0.0
                   ENDIF
                   IF(In[2] EQ -1.0 OR In[5] EQ -1.0 OR In[8] EQ -1.0) THEN BEGIN
                       In[2] = 0.0
                       In[5] = 0.0
                       In[8] = 0.0
                   ENDIF
                   IF(In[9] EQ -1.0 OR In[12] EQ -1.0 OR In[15] EQ -1.0) THEN BEGIN
                       In[9] = 0.0
                       In[12] = 0.0
                       In[15] = 0.0
                   ENDIF
                   IF(In[10] EQ -1.0 OR In[13] EQ -1.0 OR In[16] EQ -1.0) THEN BEGIN
                       In[10] = 0.0
                       In[13] = 0.0
                       In[16] = 0.0
                   ENDIF
                   IF(In[11] EQ -1.0 OR In[14] EQ -1.0 OR In[17] EQ -1.0) THEN BEGIN
                       In[11] = 0.0
                       In[14] = 0.0
                       In[17] = 0.0
                   ENDIF
                   IF(In[18] EQ -1.0 OR In[21] EQ -1.0 OR In[24] EQ -1.0) THEN BEGIN
                       In[18] = 0.0
                       In[21] = 0.0
                       In[24] = 0.0
                   ENDIF
                   IF(In[19] EQ -1.0 OR In[22] EQ -1.0 OR In[25] EQ -1.0) THEN BEGIN
                       In[19] = 0.0
                       In[22] = 0.0
                       In[25] = 0.0
                   ENDIF
                   IF(In[20] EQ -1.0 OR In[23] EQ -1.0 OR In[26] EQ -1.0) THEN BEGIN
                       In[20] = 0.0
                       In[23] = 0.0
                       In[26] = 0.0
                   ENDIF
                   I00= (In[0]+In[1]+In[2]+In[9]+In[10]+In[11]+In[18]+In[19]+In[20])
                   I1 = (In[3]+In[4]+In[5]+In[12]+In[13]+In[14]+In[21]+In[22]+In[23])
                   I2 = (In[6]+In[7]+In[8]+In[15]+In[16]+In[17]+In[24]+In[25]+In[26])
               END
               2: BEGIN
                   In[*] = Inten[iii,jjj,*]
                   IF(In[0] EQ -1.0 OR In[9] EQ -1.0 OR In[18] EQ -1.0) THEN BEGIN
                       In[0] = 0.0
                       In[9] = 0.0
                       In[18] = 0.0
                   ENDIF
                   IF(In[1] EQ -1.0 OR In[10] EQ -1.0 OR In[19] EQ -1.0) THEN BEGIN
                       In[1] = 0.0
                       In[10] = 0.0
                       In[19] = 0.0
                   ENDIF
                   IF(In[2] EQ -1.0 OR In[11] EQ -1.0 OR In[20] EQ -1.0) THEN BEGIN
                       In[2] = 0.0
                       In[11] = 0.0
                       In[20] = 0.0
                   ENDIF
                   IF(In[3] EQ -1.0 OR In[12] EQ -1.0 OR In[21] EQ -1.0) THEN BEGIN
                       In[3] = 0.0
                       In[12] = 0.0
                       In[21] = 0.0
                   ENDIF
                   IF(In[4] EQ -1.0 OR In[13] EQ -1.0 OR In[22] EQ -1.0) THEN BEGIN
                       In[4] = 0.0
                       In[13] = 0.0
                       In[22] = 0.0
                   ENDIF
                   IF(In[5] EQ -1.0 OR In[14] EQ -1.0 OR In[23] EQ -1.0) THEN BEGIN
                       In[5] = 0.0
                       In[14] = 0.0
                       In[23] = 0.0
                   ENDIF
                   IF(In[6] EQ -1.0 OR In[15] EQ -1.0 OR In[24] EQ -1.0) THEN BEGIN
                       In[6] = 0.0
                       In[15] = 0.0
                       In[24] = 0.0
                   ENDIF
                   IF(In[7] EQ -1.0 OR In[16] EQ -1.0 OR In[25] EQ -1.0) THEN BEGIN
                       In[7] = 0.0
                       In[16] = 0.0
                       In[25] = 0.0
                   ENDIF
                   IF(In[8] EQ -1.0 OR In[17] EQ -1.0 OR In[26] EQ -1.0) THEN BEGIN
                       In[8] = 0.0
                       In[17] = 0.0
                       In[26] = 0.0
                   ENDIF
                   I00= (In[0]+In[1]+In[2]+In[3]+In[4]+In[5]+In[6]+In[7]+In[8])
                   I2 = (In[9]+In[10]+In[11]+In[12]+In[13]+In[14]+In[15]+In[16]+In[17])
                   I1 = (In[18]+In[19]+In[20]+In[21]+In[22]+In[23]+In[24]+In[25]+In[26])
               END
           ENDCASE

                tot      = I00+I1+I2
                vec[0,0] = I00*3.d0/tot-1.d0
                vec[0,1] = I1 *3.d0/tot-1.d0
                vec      = REFORM(mat[i,*,*])##vec
                Bg0[iii,jjj,i]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)

               ;WE CLEAN THE SOLUTION
                IF(Bg0[iii,jjj,i] GT 1.6d0) THEN Bg0[iii,jjj,i] = 1.6d0
                Phig0[iii,jjj,i] = ATAN(vec[0,1],vec[0,0])
            ENDFOR

      ; TO RECONSTRUCT THE INTENSITIES AT THE WAVELENGTH lam2
            FOR j=0,nseq-1 DO BEGIN
                profile = 1.0;blocker0
                FOR i=0,2 DO profile = profile * (1.d0+Bg0[iii,jjj,i]*COS(2.d0*!dpi/FSR[i]*lam0+Phig0[iii,jjj,i]+tuning[i,j]))/2.d0
                profile = profile * TOTAL(Inten[iii,jjj,*])*8.d0/FLOAT(nseq) ; laser intensity * non-tunable transmission
               ;FOR i=3,6 DO profile = profile * (1.d0+0.98d0*COS(2.d0*!dpi/FSR[i]*lam0))/2.d0
                Intenrecons[iii,jjj,j] = profile
                IF(FINITE(Intenrecons[iii,jjj,j]) EQ 0) THEN Intenrecons[iii,jjj,j] = 0.0
            ENDFOR
            
        ENDIF ELSE BEGIN

            Bg0[iii,jjj,*]  = 0.d0
            Phig0[iii,jjj,*]= 0.d0

        ENDELSE
           
    ENDFOR
ENDFOR

SAVE,Bg0,Phig0,nx,FILE='RESULTS/RESULTS_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'TOTAL TIME ELAPSED',SYSTIME(1)-time0

draw:
RESTORE,'RESULTS/RESULTS_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'


; WE APPLY A MASK TO KEEP ONLY THE PHASES AND CONTRASTS ON THE ELEMENT
; APERTURE
a        = WHERE(distance GT anglim,COMPLEMENT=ba)

FOR i=0,2 DO BEGIN
    temp=REFORM(Phig0[*,*,i])
    temp[a]=-10000.d0
    phig0[*,*,i]=temp[*,*]
    temp=REFORM(Bg0[*,*,i])
    temp[a]=-10000.d0
    Bg0[*,*,i]=temp[*,*]
   ; TO GET RID OF THE STRIPES (LEFT AND UPPER EDGES)
   ;Phig0[0:5,*,i]=-10000.d0
   ;Phig0[*,122:127,i]=-10000.d0
   ;Bg0[0:5,*,i]=-10000.d0
   ;Bg0[*,122:127,i]=-10000.d0
ENDFOR

a = WHERE(FINITE(Phig0) EQ 0)
IF(a[0] NE -1) THEN Phig0[a]= -10000.0
a = WHERE(FINITE(Bg0) EQ 0)
IF(a[0] NE -1) THEN Bg0[a]  = -10000.0

a        = WHERE(REFORM(Bg0[*,*,0]) EQ -10000.d0,COMPLEMENT=ba)

; ADD 180 DEGREES TO PHASES < -180
FOR i=0,2 DO BEGIN
temp = REFORM(Phig0[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi
        Phig0[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+dpi
        Phig0[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(3)
moyb=moy
mini=DBLARR(3)
maxi=DBLARR(3)
mini2=DBLARR(3)
maxi2=DBLARR(3)
FOR i=0,2 DO BEGIN
    temp=REFORM(Phig0[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    moy[i]=MEAN(temp[a])
    mini[i]=MIN(temp[a])*180./!dpi
    maxi[i]=MAX(temp[a])*180./!dpi
    temp[b]=moy[i]
    Phig0[*,*,i]=temp
    temp=REFORM(Bg0[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    moyb[i]=MEAN(temp[a])
    temp[b]=moyb[i]
    mini2[i]=MIN(temp[a])
    maxi2[i]=MAX(temp[a])
    Bg0[*,*,i]=temp
ENDFOR


;GOTO,jose2
; WE PLOT THE RESULT
SET_PLOT,'ps'
!p.multi=[0,2,3]
device,file='yo.ps',bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27,/color
loadct,4;3

tvim,phig0[*,*,0]*180.d0/!dpi,/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[0],maxi[0]],pcharsize=1.5

tvim,phig0[*,*,1]*180.d0/!dpi,/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[1],maxi[1]],pcharsize=1.5

tvim,phig0[*,*,2]*180.d0/!dpi,/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[2],maxi[2]],pcharsize=1.5

tvim,Bg0[*,*,0],/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[0],maxi2[0]],pcharsize=1.5
temp = Bg0[*,*,0]
hist = histogram(temp[ba],binsize=0.001,min=0.7,max=1.05)
plot,FINDGEN(N_ELEMENTS(hist))*0.001+0.7,hist/FLOAT(N_ELEMENTS(ba)),psym=10,tit='!7r!17='+STRING(SIGMA(temp[ba]))+'/!7l!17='+STRING(MEAN(temp[ba])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5

tvim,Bg0[*,*,1],/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[1],maxi2[1]],pcharsize=1.5
temp = Bg0[*,*,1]
hist = histogram(temp[ba],binsize=0.001,min=0.7,max=1.05)
plot,FINDGEN(N_ELEMENTS(hist))*0.001+0.7,hist/FLOAT(N_ELEMENTS(ba)),psym=10,tit='!7r!17='+STRING(SIGMA(temp[ba]))+'/!7l!17='+STRING(MEAN(temp[ba])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5

tvim,Bg0[*,*,2],/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[2],maxi2[2]],pcharsize=1.5
temp = Bg0[*,*,2]
hist = histogram(temp[ba],binsize=0.001,min=0.7,max=1.05)
plot,FINDGEN(N_ELEMENTS(hist))*0.001+0.7,hist/FLOAT(N_ELEMENTS(ba)),psym=10,tit='!7r!17='+STRING(SIGMA(temp[ba]))+'/!7l!17='+STRING(MEAN(temp[ba])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5

TVIM,TOTAL(Inten[*,*,*],3),/scale,tit='!17Laser Intensity * n.t. transmittance',xtit='!17pixels',ytit='!17pixels'

DEVICE,/close
PRINT,'AVERAGES'

FOR i=0,2 DO PRINT,moy[i]*180.d0/!dpi,mini[i],maxi[i]
FOR i=0,2 DO PRINT,moyb[i]
jose2:


;SET_PLOT,'ps'
;!p.multi=0
;device,file=list+'_NB.ps',bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27,/color
;loadct,4;3
;tvim,phig0[*,*,0]*180.d0/!dpi-moy[0]*180.d0/!dpi,/scale,tit='!17Narrow-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[0],maxi[0]]-moy[0]*180.d0/!dpi,pcharsize=1.5
;device,/close
;
;device,file=list+'_WB.ps',bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27,/color
;tvim,phig0[*,*,1]*180.d0/!dpi-moy[1]*180.d0/!dpi,/scale,tit='!17Broad-Band Michelson',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[1],maxi[1]]-moy[1]*180.d0/!dpi,pcharsize=1.5
;device,/close
;
;device,file=list+'_E1.ps',bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27,/color
;tvim,phig0[*,*,2]*180.d0/!dpi-moy[2]*180.d0/!dpi,/scale,tit='!17Lyot Element E1',xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[2],maxi[2]]-moy[2]*180.d0/!dpi,pcharsize=1.5
;device,/close

;RESTORE,'TMPREADFILE2_BIS_FIV_512.BIN' ; 2nd detune with dark frame computed line by line and for each position size 512*512
;nx2=512
;RESTORE,'TMPREADFILE3_FIV_256.BIN'
;nx2=nx
;Inten2        = imx
;moncul=DBLARR(27)
;for i=0,26 do moncul[i]=MEAN(intenrecons[1200/(4096./nx):1760/(4096./nx),2400/(4096./nx):3040/(4096./nx),i])
;for i=0,26 do moncul[i]=MEAN(intenrecons[25:150,75:150,i])
;!p.multi=[0,1,2]
;device,file='yo2.ps',bits=24,xoffset=0,yoffset=0,xsize=15,ysize=22
;plot,moncul/MAX(moncul),tit='Reconstructed Intensities',xst=1
;for i=0,26 do moncul[i]=MEAN(inten2[1200/(4096./nx2):1760/(4096./nx2),2400/(4096./nx2):3040/(4096./nx2),i])
;for i=0,26 do moncul[i]=MEAN(inten2[25:150,75:150,i])
;plot,moncul/MAX(moncul),tit='Measured Intensities',xst=1
;device,/close


; FIND THE MOTOR POSITIONS FOR THE CO-TUNE SEQUENCE

; DEFINITION OF A CO-TUNE IN RADIANS
;!p.multi=0
;nseq2=20 ; NUMBER OF POSITION IN THE CO-TUNE
;tuning2=DBLARR(3,nseq2)
;;FOR i=0,nseq2-1 DO tuning2[*,i] = [ (DOUBLE(i)/DOUBLE(nseq2)*2.d0*!dpi-!dpi)*FSR[2]/FSR[0],(DOUBLE(i)/DOUBLE(nseq2)*2.d0*!dpi-!dpi)*FSR[2]/FSR[1],(DOUBLE(i)/DOUBLE(nseq2)*2.d0*!dpi-!dpi)*FSR[2]/FSR[2]]
;
;tuning2[*,0] =[ -30,   0,   0]
;tuning2[*,1] =[ -27,  -6, -12]
;tuning2[*,2] =[ -24, -12, -24]
;tuning2[*,3] =[ -21, -18,  24]
;tuning2[*,4] =[ -18, -24,  12]
;tuning2[*,5] =[ -15, -30,   0]
;tuning2[*,6] =[ -12,  24, -12]
;tuning2[*,7] =[  -9,  18, -24]
;tuning2[*,8] =[  -6,  12,  24]
;tuning2[*,9] =[  -3,   6,  12]
;tuning2[*,10]=[  0 ,  0 ,   0]
;tuning2[*,11]=[  3 , -6 , -12]
;tuning2[*,12]=[  6 ,-12 , -24]
;tuning2[*,13]=[  9 ,-18 ,  24]
;tuning2[*,14]=[ 12 ,-24 ,  12]
;tuning2[*,15]=[ 15 ,-30 ,   0]
;tuning2[*,16]=[ 18 , 24 , -12]
;tuning2[*,17]=[ 21 , 18 , -24]
;tuning2[*,18]=[ 24 , 12 ,  24]
;tuning2[*,19]=[ 27 ,  6 ,  12]
;tuning2      = tuning2 + [111,100,67]
;
;tuning2=REVERSE(tuning2,1)
;for i=0,nseq2-1 do tuning2[*,i]=REFORM(tuning2[*,i])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]
;
;
;nlam        = 5000    ; number of wavelengths            
;dlam        = 2.d0/1.d3        ; resolution in wavelength
;lam         = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
;SET_PLOT,'PS'
;DEVICE,file='yo.ps',bits=24,/color
;loadct,4
;FOR i=0,nseq2-1 DO BEGIN & profile=DBLARR(nlam)+1.d0 & FOR j=0,2 DO profile = profile * (1.d0+COS(2.d0*!dpi/FSR[j]*lam+tuning2[j,i]+moy[j]*!dpi/180.d0))/2.d0 & FOR j=3,6 DO profile = profile * (1.d0+COS(2.d0*!dpi/FSR[j]*lam))/2.d0 & IF i EQ 0 THEN plot,lam,profile,yrange=[0,1],yst=1,xrange=[-0.6,0.6],xst=1,col=i*254/(nseq2-1),thick=3 ELSE oplot,lam,profile,col=i*254/(nseq2-1) & ENDFOR
;DEVICE,/close

SET_PLOT,'x'
!p.multi=0

END
