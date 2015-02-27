; IDL VERSION OF THE lookup.c CODE
;
;
; PROGRAM TO GENERATE THE HMI FILTER PROFILES
; ASSUME 6 FILTERS WITH THE FOLLOWING NAMES AND CENTRAL LOCATIONS:
; I0 is centered at +170 mA                                                                      
; I1 is centered at +102 mA                                                                      
; I2 is centered at +34  mA                                                                      
; I3 is centered at -34  mA                                                                      
; I4 is centered at -102 mA                                                                      
; I5 is centered at -170 mA                                                                      
;--------------------------------------------------------------------------------------------------------------------------------


;-------------------------------------------------------------------------------------------------------------------------
;
; MAIN PROGRAM
;
; nx2 is the number of columns and rows for which the user
; wants to produce filter profiles.
; The actual calibration data used to produce the filter profiles are
; currently in 256x256 (the gradient across these data are relatively
; small, so there is no need for the full 4096x4096)
; So if nx2>256, the code will do a bilinear interpolation from a
; 256*256 grid to a nx2*nx2 one
;
;-------------------------------------------------------------------------------------------------------------------------


PRO lookup2,row,column,Num_lambda_filter

nx2=256

; Hollow Core Motors steps for the 3 tunable elements
; we assume the tuning polarizaer is fixed at 0
;----------------------------------------------------

HCME1=37
HCMWB=59
HCMNB=82


; Free Spectral Ranges of the 7 optical-filter elements
; these FSRs are ESTIMATES only, and more accurate values will be
; derived once SDO is launched
;-----------------------------------

FSR     = DBLARR(7)    ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]  = 0.172d0-0.0010576d0;0.1710098d0  ; for the narrow-band Michelson in Angstrom
FSR[1]  = 0.344d0-0.00207683d0;0.3421506d0  ; for the broad-band  Michelson
FSR[2]  = 0.693d0+0.000483467d0;0.6943613d0  ; for the Lyot element E1
FSR[3]  = 1.407d0      ; for E2
FSR[4]  = 2.779d0      ; for E3
FSR[5]  = 5.682d0      ; for E4
FSR[6]  = 11.354d0     ; for E5

lam0    = 6173.3433d0 ;target wavelength

ntune   = Num_lambda_filter           ; number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0   ; number of tuning position in each wing

contrast   = FLTARR(7)+1.0

dtune  = FSR[0]/inttune  ; interval between two tuning positions
tune   = DBLARR(3,ntune)
FOR  i = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune

; GRID WE WANT IN WAVELENGTH
dvtest        = 24.d0
dlamdv        = 2.059205672212074294d-5
nlam          = 16001   ; number of wavelength points
dlam          = dvtest*dlamdv ; sampling rate in Angstroms
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam
solarradiustable = 976.22d0
ntest         = 1599.;1399.;1399;1369;621
;vtest         = (DINDGEN(ntest)-310.)*dvtest
vtest         = (DINDGEN(ntest)-(ntest-1)/2)*dvtest


; TUNING POSITIONS OF THE NB AND WB MICHELSONS AND OF THE LYOT ELEMENT
; E1
;---------------------------------------------------------------------
HCME1phase = DBLARR(Num_lambda_filter)
HCMWBphase = DBLARR(Num_lambda_filter)
HCMNBphase = DBLARR(Num_lambda_filter)

IF(Num_lambda_filter EQ 6) THEN BEGIN
HCME1phase[0]= double( (HCME1+15)*6 MOD 360)*!pi/180.0 ; //I0
HCME1phase[1]= double( (HCME1+9 )*6 MOD 360)*!pi/180.0 ; //I1
HCME1phase[2]= double( (HCME1+3 )*6 MOD 360)*!pi/180.0 ; //I2
HCME1phase[3]= double( (HCME1-3 )*6 MOD 360)*!pi/180.0 ; //I3
HCME1phase[4]= double( (HCME1-9 )*6 MOD 360)*!pi/180.0 ; //I4
HCME1phase[5]= double( (HCME1-15)*6 MOD 360)*!pi/180.0 ; //I5
                                 
HCMWBphase[0]= double( (HCMWB-30)*6 MOD 360)*!pi/180.0 ;
HCMWBphase[1]= double( (HCMWB-18)*6 MOD 360)*!pi/180.0 ;
HCMWBphase[2]= double( (HCMWB-6 )*6 MOD 360)*!pi/180.0 ;
HCMWBphase[3]= double( (HCMWB+6 )*6 MOD 360)*!pi/180.0 ;
HCMWBphase[4]= double( (HCMWB+18)*6 MOD 360)*!pi/180.0 ;
HCMWBphase[5]= double( (HCMWB-30)*6 MOD 360)*!pi/180.0 ;
                                 
HCMNBphase[0]= double( (HCMNB+0 )*6 MOD 360)*!pi/180.0 ;
HCMNBphase[1]= double( (HCMNB+24)*6 MOD 360)*!pi/180.0 ;
HCMNBphase[2]= double( (HCMNB-12)*6 MOD 360)*!pi/180.0 ;
HCMNBphase[3]= double( (HCMNB+12)*6 MOD 360)*!pi/180.0 ;
HCMNBphase[4]= double( (HCMNB-24)*6 MOD 360)*!pi/180.0 ;
HCMNBphase[5]= double( (HCMNB+0 )*6 MOD 360)*!pi/180.0 ;
ENDIF
              
ContinuumE1  = double( (HCME1-30)*6 MOD 360)*!pi/180.0 ; //CONTINUUM
ContinuumWB  = double( (HCMWB+0 )*6 MOD 360)*!pi/180.0
ContinuumNB  = double( (HCMNB+0 )*6 MOD 360)*!pi/180.0

IF(Num_lambda_filter EQ 10) THEN BEGIN
    HCME1phase[0]= double( ((HCME1+27)*6 MOD 360))*!dpi/180.d0 ; //I9
    HCME1phase[1]= double( ((HCME1+21)*6 MOD 360))*!dpi/180.d0 ; //I7
    HCME1phase[2]= double( ((HCME1+15)*6 MOD 360))*!dpi/180.d0 ; //I0
    HCME1phase[3]= double( ((HCME1+9 )*6 MOD 360))*!dpi/180.d0 ; //I1
    HCME1phase[4]= double( ((HCME1+3 )*6 MOD 360))*!dpi/180.d0 ; //I2
    HCME1phase[5]= double( ((HCME1-3 )*6 MOD 360))*!dpi/180.d0 ; //I3
    HCME1phase[6]= double( ((HCME1-9 )*6 MOD 360))*!dpi/180.d0 ; //I4
    HCME1phase[7]= double( ((HCME1-15)*6 MOD 360))*!dpi/180.d0 ; //I5
    HCME1phase[8]= double( ((HCME1-21)*6 MOD 360))*!dpi/180.d0 ; //I6
    HCME1phase[9]= double( ((HCME1-27)*6 MOD 360))*!dpi/180.d0 ; //I8

    HCMWBphase[0]= double( ((HCMWB+6) *6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[1]= double( ((HCMWB+18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[2]= double( ((HCMWB-30)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[3]= double( ((HCMWB-18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[4]= double( ((HCMWB-6 )*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[5]= double( ((HCMWB+6 )*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[6]= double( ((HCMWB+18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[7]= double( ((HCMWB-30)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[8]= double( ((HCMWB-18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[9]= double( ((HCMWB-6) *6 MOD 360))*!dpi/180.d0 ;

    HCMNBphase[0]= double( ((HCMNB+12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[1]= double( ((HCMNB-24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[2]= double( ((HCMNB-0 )*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[3]= double( ((HCMNB+24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[4]= double( ((HCMNB-12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[5]= double( ((HCMNB+12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[6]= double( ((HCMNB-24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[7]= double( ((HCMNB+0 )*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[8]= double( ((HCMNB+24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[9]= double( ((HCMNB-12)*6 MOD 360))*!dpi/180.d0 ;
ENDIF
IF(Num_lambda_filter EQ 8) THEN BEGIN
    HCME1phase[0]= double( ((HCME1+21)*6 MOD 360))*!dpi/180.d0 ; //I7
    HCME1phase[1]= double( ((HCME1+15)*6 MOD 360))*!dpi/180.d0 ; //I0
    HCME1phase[2]= double( ((HCME1+9 )*6 MOD 360))*!dpi/180.d0 ; //I1
    HCME1phase[3]= double( ((HCME1+3 )*6 MOD 360))*!dpi/180.d0 ; //I2
    HCME1phase[4]= double( ((HCME1-3 )*6 MOD 360))*!dpi/180.d0 ; //I3
    HCME1phase[5]= double( ((HCME1-9 )*6 MOD 360))*!dpi/180.d0 ; //I4
    HCME1phase[6]= double( ((HCME1-15)*6 MOD 360))*!dpi/180.d0 ; //I5
    HCME1phase[7]= double( ((HCME1-21)*6 MOD 360))*!dpi/180.d0 ; //I6

    HCMWBphase[0]= double( ((HCMWB+18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[1]= double( ((HCMWB-30)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[2]= double( ((HCMWB-18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[3]= double( ((HCMWB-6 )*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[4]= double( ((HCMWB+6 )*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[5]= double( ((HCMWB+18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[6]= double( ((HCMWB-30)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[7]= double( ((HCMWB-18)*6 MOD 360))*!dpi/180.d0 ;

    HCMNBphase[0]= double( ((HCMNB-24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[1]= double( ((HCMNB-0 )*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[2]= double( ((HCMNB+24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[3]= double( ((HCMNB-12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[4]= double( ((HCMNB+12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[5]= double( ((HCMNB-24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[6]= double( ((HCMNB+0 )*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[7]= double( ((HCMNB+24)*6 MOD 360))*!dpi/180.d0 ;
ENDIF

; PHASE AND CONTRAST MAPS OF THE 4 NON-TUNABLE ELEMENTS (LYOT E2 TO
; E5), in 256x256
;------------------------------------------------------------------

nx=256l ; the phase and contrast maps are in 256x256 (ONLY TEMPORARY, WILL BE UPDATED LATER)
;RESTORE,'/auto/home0/couvidat/public_html/HMI/LEKA/RESULTS_June09_710660_CAL_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE.BIN' ; data from 2007/10/15 taken in vacuum with the dye laser (side camera)
; we want PHIG: the phases (in degrees), and BG: the contrasts
;I0g=0.0 ; to save some memory
;lateconv=0.0
;threshold=0.0

PHIG=DBLARR(nx,nx,7)
tempmerde=FLTARR(4l*nx*nx)
OPENR,1,"/home/jsoc/cvs/Development/JSOC/proj/lev1.5_hmi/apps/non_tunable_phases_710660_June09_cal_256_2.bin"
READU,1,tempmerde
CLOSE,1
FOR i=0l,nx-1l DO FOR j=0l,nx-1l DO FOR k=0l,3l DO PHIG[i,j,k+3]=double(tempmerde[i+j*nx+k*nx*nx])

BG=DBLARR(nx,nx,7)
tempmerde=FLTARR(4l*nx*nx)
OPENR,1,"/home/jsoc/cvs/Development/JSOC/proj/lev1.5_hmi/apps/non_tunable_contrasts_710660_June09_cal_256_2.bin"
READU,1,tempmerde
CLOSE,1
FOR i=0l,nx-1l DO FOR j=0l,nx-1l DO FOR k=0l,3l DO BG[i,j,k+3]=double(tempmerde[i+j*nx+k*nx*nx])


; PHASE MAPS OF THE TUNABLE ELEMENTS (MICHELSONS NB, WB, AND LYOT E1)
; in 256x256
;-------------------------------------------------------------------

;phases=READFITS("/auto/home0/couvidat/public_html/HMI/LEKA/PHASE_MAPS_707119.fits") ; data from 2007/10/14, taken in vacuum and sunlight
;FOR i=0,2 DO PHIG[*,*,i]=phases[i,*,*]
;phases=0.0

tempmerde=READFITS("/SUM0/D109538835/D46717947/S00000/phases.fits")
FOR i=0,2 DO PHIG[*,*,i]=double(tempmerde[i,*,*])


; CONTRAST MAPS OF THE TUNABLE ELEMENTS (MICHELSONS NB, WB, AND LYOT E1)
; in 256x256
;-------------------------------------------------------------------

;contrasts=READFITS("/auto/home0/couvidat/public_html/HMI/LEKA/CONTRAST_MAPS_710660.fits") ; data from 2007/10/15 taken in vacuum with the dye laser
;FOR i=0,2 DO BG[*,*,i]=contrasts[*,*,i]
;contrasts=0.0

tempmerde=FLTARR(3l*nx*nx)
OPENR,1,"/home/jsoc/cvs/Development/JSOC/proj/lev1.5_hmi/apps/tunable_contrasts_710660_June09_cal_256.bin"
READU,1,tempmerde
CLOSE,1
FOR i=0l,nx-1l DO FOR j=0l,nx-1l DO FOR k=0l,2l DO BG[i,j,k]=double(tempmerde[i+j*nx+k*nx*nx])



; BILINEAR INTERPOLATION OF THE PHASE AND CONTRAST MAPS
; FROM 256*256 to nx2*nx2
;------------------------------------------------------

PHIG=REBIN(PHIG,nx2,nx2,7)
BG  =REBIN(BG,nx2,nx2,7)


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;-----------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
; THIS IS THE FRONT WINDOW S/N 3?
; IN OCTOBER 2007 THE FRONT WINDOW WAS CHANGED FROM S/N 1 TO S/N 3:
transmission = DBLARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength   = REFORM(transmission[0,*])
transmission = REFORM(transmission[1,*])
blocker      = INTERPOL(transmission/100.d0,wavelength*10.d0-lam0,lam)
q            = READFITS('blocker11.fits')
blocker      = blocker * INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0-lam0,lam) ; I center the profile (VALUE OF CENTER AS OF JUNE 2009)


; SOLAR LINE
;-------------------------------------------------


; LINE FROM R. K. ULRICH
tempcul=FLTARR(7000)
OPENR,1,'/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceFeLine.bin'
READU,1,tempcul
CLOSE,1
templineref=double(tempcul)
tempcul=FLTARR(7000)
OPENR,1,'/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceWavelength.bin'
READU,1,tempcul
CLOSE,1
wavelengthref=double(tempcul)
wavelengthref=[-10.,wavelengthref,10.]
templineref=[1.0,templineref,1.0]

lineref = INTERPOL(templineref,wavelengthref,lam)

; GAUSSIAN LINE
;lineref = 1.d0-0.655957d0*exp(-lam^2.d0/0.00300482d0)

; KITT PEAK
goto,merde
kp=fltarr(3,5911)
openr,1,'KittPeakSolarAtlas.txt'
readf,1,kp
close,1
lineref = INTERPOL(kp[1,*],kp[0,*]*10.,lam+6173.3433)
merde:

IF(Num_lambda_filter EQ 6) THEN BEGIN
cosi=COS([+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0]*2d0*!dpi/6.d0)
sini=SIN([+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0]*2d0*!dpi/6.d0)
cos2i=COS(2.d0*[+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0]*2d0*!dpi/6.d0)
sin2i=SIN(2.d0*[+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0]*2d0*!dpi/6.d0)
ENDIF
IF(Num_lambda_filter EQ 10) THEN BEGIN
cosi=COS([+4.5d0,+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0,-4.5d0]*2d0*!dpi/10.d0)
sini=SIN([+4.5d0,+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0,-4.5d0]*2d0*!dpi/10.d0)
cos2i=COS(2.d0*[+4.5d0,+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0,-4.5d0]*2d0*!dpi/10.d0)
sin2i=SIN(2.d0*[+4.5d0,+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0,-4.5d0]*2d0*!dpi/10.d0)
ENDIF
IF(Num_lambda_filter EQ 8) THEN BEGIN
cosi=COS([+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0]*2d0*!dpi/8.d0)
sini=SIN([+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0]*2d0*!dpi/8.d0)
cos2i=COS(2.d0*[+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0]*2d0*!dpi/8.d0)
sin2i=SIN(2.d0*[+3.5d0,+2.5d0,+1.5d0,+0.5d0,-0.5d0,-1.5d0,-2.5d0,-3.5d0]*2d0*!dpi/8.d0)
ENDIF


minimumCoeffs=[0.41922611d0,0.24190794d0]
FWHMCoeffs=[151.34559d0,-58.521771d0]       
distance=SQRT((row-127.5d0)^2.d0+(column-127.5d0)^2.d0)*0.5d0*4096.d0/nx2
distance=cos(asin(distance/solarradiustable))
minimumI=[0.662964d0,0.585863d0,0.542768d0]
FWHMline=[93.878368d0,107.41855d0,123.57606d0]

FWHM=FWHMCoeffs[0]+FWHMCoeffs[1]*distance
minimum=minimumCoeffs[0]+minimumCoeffs[1]*distance
lineprofile2    = (1.d0-lineref)*minimum/minimumI[0] ; //scaling by the ratio of linedepths
lineprofile2    =  1.d0-lineprofile2 ;
wavelength2     = lam*FWHM/FWHMline[0] ;

lineprofile = INTERPOL(lineprofile2,wavelength2,lam)

filters     = FLTARR(nlam,ntune)
lyot        = FLTARR(nlam)

lyot = blocker

PRINT,'PHASES AND CONTRAST',PHIG[column,row,0]*!dpi/180.d0,PHIG[column,row,1]*!dpi/180.d0,PHIG[column,row,2]*!dpi/180.d0,PHIG[column,row,3]*!dpi/180.d0,PHIG[column,row,4]*!dpi/180.d0,PHIG[column,row,5]*!dpi/180.d0,PHIG[column,row,6]*!dpi/180.d0,BG[column,row,0],BG[column,row,1],BG[column,row,2],BG[column,row,3],BG[column,row,4],BG[column,row,5]

                                ;NON-TUNABLE ELEMENTS
FOR i = 0,3 DO BEGIN
    lyot = lyot*(1.d0+BG[column,row,i+3]*COS(2.d0*!dpi/FSR[i+3]*lam+PHIG[column,row,i+3]*!dpi/180.d0))/2.d0
ENDFOR
                                ;TUNABLE ELEMENTS
FOR itune = 0,ntune-1 DO BEGIN  
                              
    filters[*,itune] = lyot*(1.d0+BG[column,row,0]*COS(2.d0*!dpi/FSR[0]*lam[*]+(PHIG[column,row,0])*!dpi/180.d0+HCMNBphase[itune]))/2.d0*(1.d0+BG[column,row,1]*COS(2.d0*!dpi/FSR[1]*lam[*]+(PHIG[column,row,1])*!dpi/180.d0+HCMWBphase[itune]))/2.d0*(1.d0+BG[column,row,2]*COS(2.d0*!dpi/FSR[2]*lam[*]+((PHIG[column,row,2]))*!dpi/180.d0-HCME1phase[itune]))/2.d0
ENDFOR

filtersconti =  FLTARR(nlam) ;CONTINUUM
filtersconti[*] = lyot*(1.d0+BG[column,row,0]*COS(2.d0*!dpi/FSR[0]*lam[*]+(PHIG[column,row,0])*!dpi/180.d0+ContinuumNB))/2.d0*(1.d0+BG[column,row,1]*COS(2.d0*!dpi/FSR[1]*lam[*]+(PHIG[column,row,1])*!dpi/180.d0+ContinuumWB))/2.d0*(1.d0+BG[column,row,2]*COS(2.d0*!dpi/FSR[2]*lam[*]+((PHIG[column,row,2]))*!dpi/180.d0-ContinuumE1))/2.d0

inten=DBLARR(Num_lambda_filter)
pv1 = dtune/dlamdv*(Num_lambda_filter-1.d0) ; ERROR !!!!!! SHOULD BE Num_lambda_filter NOT Num_lambda_filter-1.d0
pv2 = pv1/2.d0
vela = DBLARR(ntest)
velb = vela

maxshiftlam=round(vtest[ntest-1]*dlamdv/dlam)
FOR i=0,ntest-1 DO BEGIN
    shiftlam=round(vtest[i]*dlamdv/dlam) ; //SIGN CONVENTION: POSITIVE VELOCITIES CORRESPOND TO REDSHIFT (MOVEMENTS AWAY FROM OBSERVER)
    ;print,vtest[i],shiftlam
    FOR j=0,Num_lambda_filter-1 DO BEGIN
        inten[j]=0.0
        FOR k=maxshiftlam,nlam-maxshiftlam-1 DO inten[j] = inten[j]+filters[k,j]*lineprofile[k-shiftlam]
    ENDFOR

    f1c=TOTAL(cosi*inten)
    f1s=TOTAL(sini*inten)
    f2c=TOTAL(cos2i*inten)
    f2s=TOTAL(sin2i*inten)

    vel1 = atan(-f1s,-f1c)*pv1/2.d0/!dpi
    vela[i]=vel1
    vel2 = atan(-f2s,-f2c)*pv2/2.d0/!dpi
    velb[i]=((vel2-vel1+10.5*pv2) MOD pv2-pv2/2.0+vel1)

    print,vtest[i],f1c,f1s,f2c,f2s

ENDFOR

; PLOT OF THE FILTER PROFILES WITH RESPECT TO THE Fe I LINE

!P.MULTI=0
SET_PLOT,'PS'
DEVICE,FILE='filterprofiles.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
LOADCT,4

PLOT,vtest,vela,xst=1,yst=1,thick=2
OPLOT,vtest,velb,LINESTYLE=2,thick=2

maxfil=MAX(filters)
PLOT,lam, filters[*,0]/maxfil,xrange=[-0.9,0.9],xst=1,charsize=1.5,tit='!17',xtit='Wavelength (A)',ytit='Normalized Intensities',thick=3,yrange=[0,1.005],yst=1
OPLOT,lam,filters[*,1]/maxfil,col=40,thick=3
OPLOT,lam,filters[*,2]/maxfil,col=80,thick=3
OPLOT,lam,filters[*,3]/maxfil,col=120,thick=3
OPLOT,lam,filters[*,4]/maxfil,col=160,thick=3
OPLOT,lam,filters[*,5]/maxfil,col=200,thick=3
IF(Num_lambda_filter EQ 10) THEN BEGIN
OPLOT,lam,filters[*,6]/maxfil,col=40,thick=3,linestyle=2
OPLOT,lam,filters[*,7]/maxfil,col=80,thick=3,linestyle=2
OPLOT,lam,filters[*,8]/maxfil,col=120,thick=3,linestyle=2
OPLOT,lam,filters[*,9]/maxfil,col=160,thick=3,linestyle=2
ENDIF
IF(Num_lambda_filter EQ 8) THEN BEGIN
OPLOT,lam,filters[*,6]/maxfil,col=40,thick=3,linestyle=2
OPLOT,lam,filters[*,7]/maxfil,col=80,thick=3,linestyle=2
ENDIF

OPLOT,lam,lineref,thick=3
OPLOT,lam,filtersconti,thick=3,linestyle=2
DEVICE,/CLOSE
SET_PLOT,'X'

READ,pause


; TO UNDERSTAND DIFFERENCES BETWEEN VELOCITY RETURNED BY FIRST AND BY
; SECOND FOURIER COEFFICIENT:
; we calculate from lev 0 images of an observable sequence of October
; 14, 2007 taken at 23:40 UT in sunlight and vacuum, the first and
; second fourier coefficients (after applying the HMI filters)
;----------------------------------------------------------------------------------------------------------------

filename=STRARR(24)
filename[0]="/SUM4/D21748040/D1988520/S00000/000706939.fits"
filename[1]="/SUM4/D21748040/D1988516/S00000/000706940.fits"
filename[2]="/SUM4/D21748040/D1988522/S00000/000706941.fits"
filename[3]="/SUM4/D21748040/D1988532/S00000/000706942.fits"
filename[4]="/SUM5/D21748041/D1988533/S00000/000706943.fits"
filename[5]="/SUM5/D21748041/D1988539/S00000/000706944.fits"
filename[6]="/SUM5/D21748041/D1988541/S00000/000706945.fits"
filename[7]="/SUM5/D21748041/D1988547/S00000/000706946.fits"
filename[8]="/SUM6/D21748042/D1988555/S00000/000706947.fits"
filename[9]="/SUM6/D21748042/D1988557/S00000/000706948.fits"
filename[10]="/SUM7/D21748043/D1988572/S00000/000706949.fits"
filename[11]="/SUM6/D21748042/D1988569/S00000/000706950.fits"
filename[12]="/SUM7/D21748043/D1988574/S00000/000706951.fits"
filename[13]="/SUM7/D21748043/D1988578/S00000/000706952.fits"
filename[14]="/SUM7/D21748043/D1988587/S00000/000706953.fits"
filename[15]="/SUM8/D21748044/D1988588/S00000/000706954.fits"
filename[16]="/SUM8/D21748044/D1988597/S00000/000706955.fits"
filename[17]="/SUM8/D21748044/D1988598/S00000/000706956.fits"
filename[18]="/SUM8/D21748044/D1988602/S00000/000706957.fits"
filename[19]="/SUM9/D21748045/D1988612/S00000/000706958.fits"
filename[20]="/SUM9/D21748045/D1988615/S00000/000706959.fits"
filename[21]="/SUM0/D21748046/D1988623/S00000/000706960.fits"
filename[22]="/SUM0/D21748046/D1988627/S00000/000706961.fits"
filename[23]="/SUM1/D21748047/D1988632/S00000/000706962.fits"

inten=fltarr(12,4096,4096)
for i=0,11 do inten[i,*,*]=readfits(filename[i*2])

vel1=fltarr(4096,4096)
vel2=vel1
FOR i=0,4095 do for j=0,4095 do begin & f1c=TOTAL(cosi *inten[[10,8,6,4,2,0],i,j]) & f1s=TOTAL(sini *inten[[10,8,6,4,2,0],i,j]) & f2c=TOTAL(cos2i*inten[[10,8,6,4,2,0],i,j]) & f2s=TOTAL(sin2i*inten[[10,8,6,4,2,0],i,j]) & vel1[i,j] = atan(-f1s,-f1c)*pv1/2.d0/!dpi & vel2[i,j] = atan(-f2s,-f2c)*pv2/2.d0/!dpi & ENDFOR

vel1[2048,*]=0.0
vel1[*,2048]=0.0
vel2[2048,*]=0.0
vel2[*,2048]=0.0
vel1[2047,*]=0.0
vel1[*,2047]=0.0
vel2[2047,*]=0.0
vel2[*,2047]=0.0


distance=fltarr(4096,4096)
for i=0,4095 do for j=0,4095 do distance[i,j]=sqrt( (i-2073.)^2.+(j-1989.)^2. )
a=where(distance le 1870.,complement=b)
distance[a]=1.0
distance[b]=0.0

temp=vel2-vel1
temp2=temp*distance
rop=FINDGEN(4096l*128l)*32.
vel1=vel1*distance
vel2=vel2*distance

SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,/color,bits=24
LOADCT,3
!p.multi=[0,1,2]
john=WHERE(vel1[rop] NE 0.0 AND vel2[rop] NE 0.0)
res=poly_fit(vel1[rop[john]],vel2[rop[john]],1,yfit=y)
plot,vel1[rop],vel2[rop],xst=1,psym=3,xrange=[-2000,2000],tit='!17 From observable sequence',xtit='velocity 1st Fourier (m/s)',ytit='velocity 2nd Fourier (m/s)',charsize=1.5,yrange=[-2000,3000],yst=1
OPLOT,[-2000,2000],[0,0],col=180,thick=2,linestyle=2
OPLOT,[0,0],[-6000,6000],col=180,thick=2,linestyle=2
OPLOT,vel1[rop[john]],y,col=180,thick=2,linestyle=2
OPLOT,[-2000,2000],[-2000,2000],col=180,thick=2,linestyle=2
plot,vela,velb,xrange=[-2000,2000],tit='!17 From look-up tables',xtit='velocity 1st Fourier (m/s)',ytit='velocity 2nd Fourier (m/s)',charsize=1.5,xst=1,yrange=[-2000,3000],yst=1,thick=2
OPLOT,[-2000,2000],[0,0],col=180,thick=2,linestyle=2
OPLOT,[0,0],[-6000,6000],col=180,thick=2,linestyle=2
OPLOT,[-2000,2000],[-2000,2000],col=180,thick=2,linestyle=2
DEVICE,/CLOSE
SET_PLOT,'x'
!p.multi=0

READ,pause

;rangx=0.170
;a=where(lam ge -rangx and lam le rangx)
;
;f1=fltarr(ntest)
;f2=f1
;FOR i=0,ntest-1 DO BEGIN &  shiftlam=round(vtest[i]*dlamdv/dlam) &  f1[i]=atan(-total(sin(lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]),-total(cos(lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]))*2.*rangx/dlamdv/2./!pi & f2[i]=atan(-total(sin(2.*lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]),-total(cos(2.*lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]))*2.*rangx/dlamdv/4./!pi & ENDFOR
;oplot,f1,f2-f1,col=180
;
;
;; create a more asymmetric line to see how it affects f1 and f2
;lineref2=lineref*(1.0-0.1*exp(-(lam-0.075)^2./0.001))
;
;lineref2=1.-0.4*exp(-(lam)^2./0.005)
;lineref2=lineref2*(1.0-0.025*exp(-(lam+0.2)^2./0.001))
;f1b=fltarr(ntest)
;f2b=f1b
;FOR i=0,ntest-1 DO BEGIN &  shiftlam=round(vtest[i]*dlamdv/dlam) &  f1b[i]=atan(-total(sin(lam[a]/2./rangx*2.*!pi)*lineref2[a-shiftlam]),-total(cos(lam[a]/2./rangx*2.*!pi)*lineref2[a-shiftlam]))*2.*rangx/dlamdv/2./!pi & f2b[i]=atan(-total(sin(2.*lam[a]/2./rangx*2.*!pi)*lineref2[a-shiftlam]),-total(cos(2.*lam[a]/2./rangx*2.*!pi)*lineref2[a-shiftlam]))*2.*rangx/dlamdv/4./!pi & ENDFOR
;oplot,f1b,f2b-f1b,col=180,linestyle=2
;

; using kitt-peak line

kp=fltarr(3,5911)
openr,1,'KittPeakSolarAtlas.txt'
readf,1,kp
close,1
lineref = INTERPOL(kp[1,*],kp[0,*]*10.,lam+6173.3433)
;f1c=fltarr(ntest)
;f2c=f1c
;FOR i=0,ntest-1 DO BEGIN &  shiftlam=round(vtest[i]*dlamdv/dlam) &  f1c[i]=atan(-total(sin(lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]),-total(cos(lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]))*2.*rangx/dlamdv/2./!pi & f2c[i]=atan(-total(sin(2.*lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]),-total(cos(2.*lam[a]/2./rangx*2.*!pi)*lineref[a-shiftlam]))*2.*rangx/dlamdv/4./!pi & ENDFOR
;oplot,f1c,f2c-f1c,col=180,linestyle=2

;using the HMI filters instead of the actual integration to calculate
;the Fourier coefficients:
;--------------------------------------------------------------------------------------

; LINE FROM R. K. ULRICH
templineref=FLTARR(7000)
OPENR,1,'/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceFeLine.bin'
READU,1,templineref
CLOSE,1
wavelengthref=FLTARR(7000)
OPENR,1,'/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceWavelength.bin'
READU,1,wavelengthref
CLOSE,1
wavelengthref=[-10.,wavelengthref,10.]
templineref=[1.0,templineref,1.0]
lineref = INTERPOL(templineref,wavelengthref,lam)

f1d=fltarr(ntest)
f2d=f1d
xlam=FINDGEN(nlam)
maxshiftlam=round(vtest[ntest-1]*dlamdv/dlam)
FOR i=0,ntest-1 DO BEGIN &  shiftlam=round(vtest[i]*dlamdv/dlam) &  I0=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,0]) & I1=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,1]) & I2=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,2]) & I3=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,3]) & I4=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,4]) & I5=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,5]) & f1d[i]=atan(-(I0*sini[0]+I1*sini[1]+I2*sini[2]+I3*sini[3]+I4*sini[4]+I5*sini[5]),-total(I0*cosi[0]+I1*cosi[1]+I2*cosi[2]+I3*cosi[3]+I4*cosi[4]+I5*cosi[5]))*pv1/2./!pi & f2d[i]=atan(-(I0*sin2i[0]+I1*sin2i[1]+I2*sin2i[2]+I3*sin2i[3]+I4*sin2i[4]+I5*sin2i[5]),-total(I0*cos2i[0]+I1*cos2i[1]+I2*cos2i[2]+I3*cos2i[3]+I4*cos2i[4]+I5*cos2i[5]))*pv2/2./!pi  & ENDFOR
f2d=((f2d-f1d+10.5*pv2) MOD pv2-pv2/2.0+f1d)
plot,f1d,f2d-f1d,xst=1

; using kitt-peak line

kp=fltarr(3,5911)
openr,1,'KittPeakSolarAtlas.txt'
readf,1,kp
close,1
lineref = INTERPOL(kp[1,*],kp[0,*]*10.,lam+6173.3433)
f1e=fltarr(ntest)
f2e=f1e
xlam=FINDGEN(nlam)
maxshiftlam=round(vtest[ntest-1]*dlamdv/dlam)
FOR i=0,ntest-1 DO BEGIN &  shiftlam=round(vtest[i]*dlamdv/dlam) &  I0=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,0]) & I1=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,1]) & I2=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,2]) & I3=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,3]) & I4=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,4]) & I5=TOTAL(lineref[maxshiftlam-shiftlam:nlam-maxshiftlam-1-shiftlam]*filters[maxshiftlam:nlam-maxshiftlam-1,5]) & f1e[i]=atan(-(I0*sini[0]+I1*sini[1]+I2*sini[2]+I3*sini[3]+I4*sini[4]+I5*sini[5]),-total(I0*cosi[0]+I1*cosi[1]+I2*cosi[2]+I3*cosi[3]+I4*cosi[4]+I5*cosi[5]))*pv1/2./!pi & f2e[i]=atan(-(I0*sin2i[0]+I1*sin2i[1]+I2*sin2i[2]+I3*sin2i[3]+I4*sin2i[4]+I5*sin2i[5]),-total(I0*cos2i[0]+I1*cos2i[1]+I2*cos2i[2]+I3*cos2i[3]+I4*cos2i[4]+I5*cos2i[5]))*pv2/2./!pi  & ENDFOR
f2e=((f2e-f1e+10.5*pv2) MOD pv2-pv2/2.0+f1e)
oplot,f1e,f2e-f1e,col=180,linestyle=2

;measurement of bisector of solar lines (line of R. K. Ulrich at disk
;center, line from Kitt Peak Atlas, and the line of R.K. Ulrich after
;I add some asymmetry to have a bisector closer to what Dravins 1981
;shows in his figure 6b
;-------------------------------------------------------
; LINE FROM R. K. ULRICH
; GRID WE WANT IN WAVELENGTH
dvtest        = 24.d0
dlamdv        = 2.059205672212074294d-5
nlam          = 16001   ; number of wavelength points
dlam          = dvtest*dlamdv ; sampling rate in Angstroms
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam


templineref=FLTARR(7000)
OPENR,1,'/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceFeLine.bin'
READU,1,templineref
CLOSE,1
wavelengthref=FLTARR(7000)
OPENR,1,'/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceWavelength.bin'
READU,1,wavelengthref
CLOSE,1
wavelengthref=[-10.,wavelengthref,10.]
templineref=[1.0,templineref,1.0]
lineref = INTERPOL(templineref,wavelengthref,lam)

inten=findgen(1000)/999. 
bis=fltarr(1000)
aa=where(lam le 0.0 and lam ge -0.45)
bb=where(lam ge 0.0 and lam le 0.45)

for i=0,999 do begin & a=where(abs(lineref[aa]-inten[i]) eq min(abs(lineref[aa]-inten[i]))) & b=where(abs(lineref[bb]-inten[i]) eq min(abs(lineref[bb]-inten[i]))) & bis[i]=(lam[aa[a]]+lam[bb[b]])/2. & endfor
plot,bis/dlamdv,inten,xrange=[-0.007,0.007]/dlamdv,xst=1,yrange=[0.35,0.98],yst=1

; KITT PEAK
kp=fltarr(3,5911)
openr,1,'KittPeakSolarAtlas.txt'
readf,1,kp
close,1
lineref = INTERPOL(kp[1,*],kp[0,*]*10.,lam+6173.3433)
bis2=fltarr(1000)
for i=0,999 do begin  & a=where(abs(lineref[aa]-inten[i]) eq min(abs(lineref[aa]-inten[i]))) & b=where(abs(lineref[bb]-inten[i]) eq min(abs(lineref[bb]-inten[i]))) & bis2[i]=(lam[aa[a]]+lam[bb[b]])/2. & endfor
oplot,bis2/dlamdv,inten,col=180
oplot,[0,0],[0,1]
;CUSTOM LINE
lineref = INTERPOL(templineref,wavelengthref,lam+0.003)
lineref = lineref*(1.0-0.08*exp(-(lam+0.045)^2./0.001))
bis3=fltarr(1000)
for i=0,999 do begin  & a=where(abs(lineref[aa]-inten[i]) eq min(abs(lineref[aa]-inten[i]))) & b=where(abs(lineref[bb]-inten[i]) eq min(abs(lineref[bb]-inten[i]))) & bis3[i]=(lam[aa[a]]+lam[bb[b]])/2. & endfor
oplot,bis3/dlamdv,inten,col=180


READ,pause












a=where(abs(vela) eq min(abs(vela)))
print,vela[a],vtest[a]
print,vela[310],vtest[310]
READ,pause


END
