; from Richard Wachter

function solar_rotation, size,bin, radius, cent_x, cent_y, p0, b0, t, shift=shift

  fac=bin
  n=size/bin

  solar_radius=6.9894e8
  shift=dblarr(n,n,3)
  vrot=dblarr(n,n)


  pang_matrix=double([[1.0, 0.0, 0.0], [0.0, cos(p0), sin(p0)], [0.0, -sin(p0), cos(p0)]])  ; p-angle rotation matrix (rotate around x)
  bang_matrix=double([[cos(b0), 0.0, -sin(b0)], [0.0, 1.0, 0.0], [sin(b0), 0.0, cos(b0)]])  ; b-angle rotation matrix (rotate around y)
  pang_invert=double([[1.0, 0.0, 0.0], [0.0, cos(p0), -sin(p0)], [0.0, sin(p0), cos(p0)]])  ; inverse p-angle rotation matrix (rotate around x)
  bang_invert=double([[cos(b0), 0.0, sin(b0)], [0.0, 1.0, 0.0], [-sin(b0), 0.0, cos(b0)]])  ; inverse b-angle rotation matrix (rotate around y)

  for i=0, n-1 do begin  
   for j=0, n-1 do begin 

    if (j*fac+fac/2) lt size then begin
     
      xy=transpose([((double(i)+0.5)*fac)-cent_x, ((double(j)+0.5)*fac-cent_y)]/radius) ; center of pixel coordinate
      inr=sqrt(xy[0]^2+xy[1]^2) 
 

     if inr le 1.0 then begin 
      xyz=transpose([sqrt(1.d0-xy[0]^2-xy[1]^2), xy[0], xy[1]]) 
   
  

      xyp=pang_matrix ## xyz
      ; b-angle rotation
      xyzb=bang_matrix ## xyp   

      sinlat=xyzb[2]
      Omega = 452.0 - 49.0*sinlat^2. - 84.0*sinlat^4. - 31.7 ; Differential rotation
      rotang=2.*!pi*Omega*1e-9*t ; rotation angle
 
      rot_matrix=double([[cos(rotang), -sin(rotang), 0.0], [sin(rotang), cos(rotang), 0.0], [0.0, 0.0, 1.0]])  ; rotate around z
 
      xyzr=rot_matrix ## xyzb

      xyzrbi=bang_invert ## xyzr ; b-angle
      xyzrpi=pang_invert ## xyzrbi

      vr=xyzrpi[0]
      xyr=transpose(xyzrpi[[1,2]]) 


      shift[i,j,0:1]=(xyr-xy)*radius
      shift[i,j,2]=(vr-xyz[0])*solar_radius/t
    endif
  endif


  endfor
 endfor

 return,shift
 end
