; program to write the listLaserXXX files for the HMI_laser2.pro code

PRO writelistlaser,index,n

filename='listLaser070225_'+STRTRIM(STRING(index),1)

OPENW,1,filename
FOR i=0,n-1 DO PRINTF,1,'/scr20/couvidat/HMI/TEMP/070225/0000'+STRTRIM(STRING(index+i),1)+'.fits'
CLOSE,1

END
