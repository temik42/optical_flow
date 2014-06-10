;function opt_fl, d1, d2, d3, ro
;
;dim=size(d1, /dim)
;
;;q1=find_shift(d1,d2)
;;q2=find_shift(d2,d3)
;;
;;d1=float_shift(d1,-q1.dx,-q1.dy)
;;d3=float_shift(d3,q2.dx,q2.dy)
;
;g_x=derivXY(d2)
;g_y=derivXY(d2, /y)
;g_t=(d3-d1)*0.5
;
;xind=replicate(1,dim[1])##(indgen(dim[0])-dim[0]/2)
;yind=(indgen(dim[1])-dim[1]/2)##replicate(1,dim[0])
;
;w=exp(-(xind^2+yind^2)/ro^2.)
;w=shift(w,-dim[0]/2,-dim[1]/2)
;fw=conj(fft(w))
;
;a11=real_part(fft(fft(g_x^2)*fw,1))
;a22=real_part(fft(fft(g_y^2)*fw,1))
;a12=real_part(fft(fft(g_x*g_y)*fw,1))
;
;det_a=a11*a22-a12^2
;
;b1=real_part(fft(fft(g_t*g_x)*fw,1))
;b2=real_part(fft(fft(g_t*g_y)*fw,1))
;
;p1=(a12*b2-a22*b1)/(det_a>1)
;p2=(a12*b1-a11*b2)/(det_a>1)
;
;return, {vx:p1,vy:p2}
;
;end



function opt_fl, d1, d2, d3, x, y, ro

dim=size(d1, /dim)
n=n_elements(x)

g_x=derivXY(d2)
g_y=derivXY(d2, /y)
g_t=(d3-d1)*0.5

a11=g_x^2
a22=g_y^2
a12=g_x*g_y

b1=g_t*g_x
b2=g_t*g_y

xind=replicate(1,dim[1])##(indgen(dim[0])-dim[0]/2)
yind=(indgen(dim[1])-dim[1]/2)##replicate(1,dim[0])

w=exp(-(xind^2+yind^2)/ro^2.)
w/=total(w)

p1=fltarr(n)
p2=fltarr(n)

for i=0, n-1 do begin

w1=shift(w,x[i]-dim[0]/2,y[i]-dim[1]/2)
;fw=conj(fft(w))

;a11=real_part(fft(fft(g_x^2)*fw,1))
;a22=real_part(fft(fft(g_y^2)*fw,1))
;a12=real_part(fft(fft(g_x*g_y)*fw,1))

aw11=total(a11*w1)
aw22=total(a22*w1)
aw12=total(a12*w1)

det_aw=aw11*aw22-aw12^2

;b1=real_part(fft(fft(g_t*g_x)*fw,1))
;b2=real_part(fft(fft(g_t*g_y)*fw,1))

bw1=total(b1*w1)
bw2=total(b2*w1)

p1[i]=(aw12*bw2-aw22*bw1)/(det_aw>1)
p2[i]=(aw12*bw1-aw11*bw2)/(det_aw>1)

endfor

return, {vx:p1,vy:p2}

end

;dir='D:\work\fd_M_96m_01d.006018\'
dir='D:\work\mdi_20090624_out\'
fnames=file_search(dir+'*.fits', count=n)
;end
window, 0, xs=512, ys=512

vect=fltarr(3,256,256,n-2)

for i=1, n-2 do begin

if i eq 1 then begin
d1=rfits(fnames[0]) > (- 1e3)
;d1=extrac(d1, 390, 680, 128, 128)

d2=rfits(fnames[1]) > (- 1e3)
;d2=extrac(d2, 390, 680, 128, 128)

;q1=find_shift(d1,d2)
;d1=float_shift(d1,-q1.dx,-q1.dy)
endif

d3=rfits(fnames[i+1]) > (- 1e3)

;d3=extrac(d3, 390-6.3*i, 680, 128, 128)

;q1=find_shift(d1,d2)
;q2=find_shift(d2,d3)

;d1=float_shift(d1,-q1.dx,-q1.dy)
;d3=float_shift(d3,q2.dx,q2.dy)
;

res=opt_fl(d1,d2,d3,5)
vect[0,*,*,i-1]=res.vx*8
vect[1,*,*,i-1]=res.vy*8
vect[2,*,*,i-1]=replicate(1.,256,256)
tvscl, d2

d2=d3
d1=d2
endfor

nx=32l
ny=32l

xind=replicate(1,ny)##indgen(nx)*8
yind=indgen(ny)##replicate(1,nx)*8

PARTICLE_TRACE, vect, transpose([[reform(xind,nx*ny)],[reform(yind,nx*ny)],[bytarr(nx*ny)]]), Verts, Conn
end
tvscl, congrid(d2,512,512,/interp,cubic=-0.5), channel=3
tvscl, congrid(d3,512,512,/interp,cubic=-0.5), channel=2
tvscl, bytarr(512,512), channel=1
plot, [0],[0], position=[0,0,1,1], /normal, xrange=[0,32], yrange=[0,32], xstyle=1, ystyle=1, /noe
velovect, rebin(res.vx,32,32), rebin(res.vy,32,32),/over
;tvscl, data2

;window, 1, xs=512, ys=512
;tvscl, congrid(d3-d1,512,512,/interp,cubic=-0.5)


end