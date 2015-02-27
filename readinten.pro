; PROGRAM TO READ THE LASER INTENSITY MONITORING
; FILE

PRO readinten

;nlines = 167738l
;nlines=106263
;nlines = 72605l
nlines=106592
line   = 'lll'

data=STRARR(3,nlines)

;OPENR,1,'~schou/public_html/hmi/suntest/meter/20070221_LASER_NEWPORT.dat'
;OPENR,1,'~schou/public_html/hmi/suntest/meter/20070222_LASER_NEWPORT.dat'
;OPENR,1,'~schou/public_html/hmi/suntest/meter/20070223_LASER_NEWPORT.dat'
;OPENR,1,'~/20070223_LASER_NEWPORT.dat'
OPENR,1,'~schou/public_html/hmi/suntest/meter/20070225_LASER_NEWPORT.dat'

FOR i=0l,nlines-1 DO BEGIN
    READF,1,line
    data[*,i] = STRSPLIT(line,/extract)
ENDFOR
CLOSE,1

time  = STRARR(3,nlines)
intensity = STRARR(nlines)

FOR i=0l,nlines-1 DO BEGIN
    time[*,i] = STRSPLIT(data[1,i],":",/extract)
    intensity[i] = STRSPLIT(data[2,i],":",/extract)
ENDFOR

hour = LONG(time[0,*])
minute=LONG(time[1,*])
second=LONG(time[2,*])

a = WHERE(intensity EQ 'Connect')
IF(a[0] NE -1) THEN intensity[a]='-1000.0'
intensity = FLOAT(intensity)
a = WHERE(intensity EQ -1000.0)
IF(a[0] NE -1) THEN intensity[a] = !VALUES.F_NAN

READ,pause

END
