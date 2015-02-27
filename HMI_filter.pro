; draw the filter positions for the HMI filter profiles


PRO HMI_filter

table= FLTARR(3,20)
table= [ [-30,0,0], $
         [-27,-6,-12], $
         [-24,-12,-24], $
         [-21,-18,24], $
         [-18,-24,12], $
         [-15,-30,0], $
         [-12,24,-12], $
         [-9,18,-24], $
         [-6,12,24], $
         [-3,6,12], $
         [0,0,0], $
         [3,-6,-12], $
         [6,-12,-24], $
         [9,-18,24], $
         [12,-24,12], $
         [15,-30,0], $
         [18,24,-12], $
         [21,18,-24], $
         [24,12,24], $
         [27,6,12] ]


HCMNB=82
HCMWB=58
HCME1=36


for i=0,19 do begin
    table[0,i]=table[0,i]+HCME1
    table[1,i]=table[1,i]+HCMWB
    table[2,i]=table[2,i]+HCMNB
endfor


lam0    = 6173.3433d0
phase   = FLTARR(7)
;phase      =[0.,0.,0.,-6.13,-0.174,-5.34,1.14]*!dpi/180.d0
phase   =[-143.364,4.70287,-140.510,-6.60914,-0.154691,-8.61043,-5.15019]*!dpi/180.d0 ; at disk center /home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/non_tunable_phases_710660_June09_cal_256_2.bin
correction=[0.,0.,0.,0.2,0.2,0.4,-1.1]

contrast   = FLTARR(7)+1.0
contrast   =[0.976811,0.999908,0.976513,0.986981,0.965252,0.989889,0.999402] ; at disk center home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/non_tunable_contrasts_710660_June09_cal_256_2.bin and /home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/tunable_contrasts_710660_June09_cal_256.bin


wmich      = [0.172d0-0.0010576d0,0.344d0-0.00207683d0,0.693d0+0.000483467d0]          ; FSRs of the tunable elements
;wmich      = [0.169,0.337d0,0.695d0]    
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0] ; FSRs non-tunable elements
dtune=wmich[0]*12.d0/360.d0*6.d0
PRINT,"DTUNE=",dtune
; GRID WE WANT IN WAVELENGTH
nlam          = 90000.
dlam          = .5d0/1.75d3
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

dlamdv        = 2.059205672212074294d-5
dvdlam        = 1.d0/dlamdv
vel           = lam*dvdlam   ; DOPPLER shift in cm/s
dvel          = dlam*dvdlam

; BLOCKER FILTER
RESTORE,'frontwindow3.bin' ; front window
blocker       = INTERPOL(transmission/100.d0,wavelength*10.d0-lam0,lam)
q             = READFITS('blocker11.fits')
blocker       = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+2.6-lam0,lam) ; blocker used for look-up tables

;a=where(blocker eq max(blocker))
;blocker=shift(blocker,45000-a[0])

lyot          = blocker
FOR i = 0,3 DO BEGIN
    lyot      = lyot   *(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]+correction[i+3]))/2.d0
ENDFOR


cmich     = 2.d0*!dpi/wmich
signe     = [1,1,-1]
;fuck=[4.,2.,1.]

 filters   = DBLARR(nlam,20)
 maxi=LONARR(20)
 maxi2=FLTARR(20)
 integ=FLTARR(20)
            FOR itune = 0,19 DO BEGIN
                filters[*,itune] = lyot
                ;FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*lam+(table[2-i,itune]-60.d0)*signe[i]*6.d0*!dpi/180.d0+phase[i]-9.*fuck[i]*!dpi/180.d0))/2.d0
                FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*lam+(table[2-i,itune]-60.d0)*signe[i]*6.d0*!dpi/180.d0+phase[i]))/2.d0
                ;filters[*,itune] = shift(filters[*,itune],-8.5*4./360.*0.170/dlam)
                integ[itune]=total(filters[*,itune])*dlam
                temp=where(filters[*,itune] EQ MAX(filters[*,itune]))
                maxi[itune]=temp[0]
                maxi2[itune]=MAX(filters[*,itune])
            ENDFOR

; 6 HMI filters (for BOB STEIN)
filter=fltarr(90000,6)
filter[*,0]=filters[*,15]
filter[*,1]=filters[*,13]
filter[*,2]=filters[*,11]
filter[*,3]=filters[*,9]
filter[*,4]=filters[*,7]
filter[*,5]=filters[*,5]
wavelength=FLOAT(lam)

READ,pause

writefits,'filters.fits',filter
 writefits,'wavelength.fits',wavelength

SET_PLOT,'PS'
device,file='yo.ps',xoffset=0,yoffset=0.5,xsize=21,ysize=15,/color,bits=24
LOADCT,4
plot,lam,filters[*,0],xst=1,xrange=[-0.4,0.4],yrange=[0,0.65],yst=1,xtit="Wavelength (A)"
i=0
xyouts,lam[maxi[i]]-0.01,maxi2[i]+0.01,STRTRIM(STRING(i),1),charsize=1.5
FOR i=1,19 DO BEGIN
    oplot,lam,filters[*,i]
    xyouts,lam[maxi[i]]-0.01,maxi2[i]+0.01,STRTRIM(STRING(i),1),charsize=1.5
ENDFOR
for i=0,19 do oplot,[i-10,i-10]*dtune,[0,1]

oplot,lam,filters[*,5],col=40,thick=3
xyouts,lam[maxi[5]]-0.01,maxi2[5]+0.05,'I5',charsize=2,color=40

oplot,lam,filters[*,7],col=80,thick=3
xyouts,lam[maxi[7]]-0.01,maxi2[7]+0.05,'I4',charsize=2,color=80

oplot,lam,filters[*,9],col=120,thick=3
xyouts,lam[maxi[9]]-0.01,maxi2[9]+0.05,'I3',charsize=2,color=120

oplot,lam,filters[*,11],col=160,thick=3
xyouts,lam[maxi[11]]-0.01,maxi2[11]+0.05,'I2',charsize=2,color=160

oplot,lam,filters[*,13],col=200,thick=3
xyouts,lam[maxi[13]]-0.01,maxi2[13]+0.05,'I1',charsize=2,color=200

oplot,lam,filters[*,15],col=240,thick=3
xyouts,lam[maxi[15]]-0.01,maxi2[15]+0.05,'I0',charsize=2,color=240

device,/close

; CALCULATE THE INTENSITY RETURNED FOR THE FURTHER-MOST FILTERS AS A
; FUNCTION OF DOPPLER VELOCITY
;------------------------------------------------------------------------------------------------

; LINE PROFILE FROM ROGER ULRICH'S WEBSITE
OPENR,1,'Ulrich_Fe_0.txt'
data = FLTARR(2,98)
READF,1,data
lamp = REFORM(data[0,*])
nlamp= N_ELEMENTS(lamp)
line = REFORM(data[1,*])
line[97]=1.0
close,1
ntest         = 101 ; Set test velocities: we will use ntest different velocities, resulting
                                ; in ntest Doppler shifted line profiles
dvtest        = 500.; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1)/2)
lines         = DBLARR(nlam,ntest)
dlinesdv      = DBLARR(nlam,ntest)
q1            = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
q2            = [MIN(lamp)-10,lamp,MAX(lamp)+10] ; Following two lines take care of continuum outside of where profile is given
;FOR i = 0,ntest-1 DO lines(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)

line0 = INTERPOL(q1,q2,lam)

;inten=FLTARR(ntest)
;FOR i=0,ntest-1 DO inten[i]=TOTAL(lines[*,i]*filters[*,0])

device,file='yo2.ps',xoffset=0,yoffset=0.5,xsize=21,ysize=15,/color,bits=24
LOADCT,4
roger=MAX(filters)
plot ,lam,filters[*,5]/roger ,xst=1,xrange=[-0.9,0.9],yrange=[0,1.0],yst=1,xtit="Wavelength (A)",thick=3,charsize=1.5
oplot,lam,filters[*,7]/roger,col=80,thick=3
oplot,lam,filters[*,9]/roger ,col=120,thick=3
oplot,lam,filters[*,11]/roger,col=160,thick=3
oplot,lam,filters[*,13]/roger,col=200,thick=3
oplot,lam,filters[*,15]/roger,col=240,thick=3
oplot,lam,line0,thick=3,linestyle=3
device,/close



SET_plot,'x'



read,pause

END
