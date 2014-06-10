
function opt_fl, d, rx, rt

dim=size(d, /dim)

g_x=(shift(d,-1,0,0)-shift(d,1,0,0))*0.5
g_y=(shift(d,0,-1,0)-shift(d,0,1,0))*0.5
g_t=(shift(d,0,0,-1)-shift(d,0,0,1))*0.5


ind=array_indices(d,lindgen(dim[0]*dim[1]*dim[2]))
xind=reform(ind[0,*],dim[0],dim[1],dim[2])-dim[0]/2
yind=reform(ind[1,*],dim[0],dim[1],dim[2])-dim[1]/2
zind=reform(ind[2,*],dim[0],dim[1],dim[2])-dim[2]/2

w=exp(-(xind^2.+yind^2.)/rx^2.-zind^2./rt^2.)
w=shift(w,-dim[0]/2,-dim[1]/2,-dim[2]/2)
fw=conj(fft(w))

a11=real_part(fft(fft(g_x^2)*fw,1))
a22=real_part(fft(fft(g_y^2)*fw,1))
a12=real_part(fft(fft(g_x*g_y)*fw,1))


det_a=a11*a22-a12^2


b1=real_part(fft(fft(g_t*g_x)*fw,1))
b2=real_part(fft(fft(g_t*g_y)*fw,1))

p1=(a12*b2-a22*b1)/(det_a>0.1)
p2=(a12*b1-a11*b2)/(det_a>0.1)

return, {vx:p1,vy:p2}

end

;dir='D:\work\fd_M_96m_01d.006018\'
;dir='D:\work\mdi_20090624_out\'
;;dir='F:\fits\'
;fnames=file_search(dir+'*.fit', count=n)
;end

if ~keyword_set(data) then begin
restore, 'D:\work\hi-c\hi-c_1.sav'
;
data=rebin(data,1024,1024,27)
data=extrac(data, 256, 384, 0, 256, 256, 27)
;data=extrac(data, 200, 500, 0, 256, 256, 27)
endif
;goto, q

;hmi=rfits('D:\work\hi-c\hmi.M_45s.20120711_183000_TAI.2.magnetogram.fits', header)
;crpix1=sxpar(header,'crpix1')
;crpix2=sxpar(header,'crpix2')
;crota2=sxpar(header, 'crota2')
;temp1=rot(hmi,crota2,1/1.2,crpix1,crpix2,/int, cubic=-0.5)
;temp1=extrac(rot(temp1,0,1024/176.,400*4+22+90,300*4+150+91.5+12,/int, cubic=-0.5),1536,1536,1024,1024)
;
;end


;res=opt_fl(float(data),5,5)

;vect=fltarr(3,256,256,27)
;
;vect[0,*,*,*]=res.vx*8
;vect[1,*,*,*]=res.vy*8
;vect[2,*,*,*]=replicate(1.,256,256,27)

;vect=fltarr(2,256,256)
;
;vect[0,*,*]=res.vx[*,*,10]*8
;vect[1,*,*]=res.vy[*,*,10]*8

ro=5

;temp=gausmooth(float(data[*,*,10]),ro)
temp=float(data[*,*,10])

ax=derivXY(temp)
ay=derivXY(temp, /y)

bxx=gausmooth(ax^2,ro)
byy=gausmooth(ay^2,ro)
bxy=gausmooth(ax*ay,ro)

q1=((bxx+byy)+sqrt((bxx-byy)^2.+bxy^2.))/2
q2=((bxx+byy)-sqrt((bxx-byy)^2+bxy^2))/2


rx=((bxx-byy)-sqrt((bxx-byy)^2.+bxy^2.))/2
ry=bxy


r=sqrt(rx^2+ry^2)

rx/=r
ry/=r

vect=fltarr(2,256,256)

vect[0,*,*]=rx*8
vect[1,*,*]=ry*8

nx=64l
ny=64l

xind=replicate(1,ny)##indgen(nx)*4
yind=indgen(ny)##replicate(1,nx)*4


window, 0, xs=512, ys=512
tvscl, congrid(data[*,*,10], 512, 512, /int)
velovect, rx[xind,yind]*8, ry[xind,yind]*8, indgen(64)*8, indgen(64)*8, /over, /dev, color='000000'xl


end
q:



;PARTICLE_TRACE, vect, transpose([[reform(xind,nx*ny)],[reform(yind,nx*ny)],[bytarr(nx*ny)]]), Verts, Conn
;verts[2,*]=1
PARTICLE_TRACE, vect, transpose([[reform(xind,nx*ny)],[reform(yind,nx*ny)]]), Verts, Conn
verts=[verts,transpose(replicate(1.,n_elements(verts[0,*])))]



end
i=0
nconn=n_elements(conn)


while i lt nconn do begin
connt=extrac(conn,i,conn[i])
vertt=verts[*,connt]

i+=conn[i]+1
endwhile



end
