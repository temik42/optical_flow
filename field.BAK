

xs=1200
ys=800
dx=0
dy=0
scale=.3

base=widget_base(XSIZE=xs, YSIZE=ys,/TLB_KILL_REQUEST_EVENTS)
;itwin=obj_new('IDLitWindow', location=[100,0], dimensions=[600,600])
draw=widget_draw(xsize=xs, ysize=ys, xoffset=0, yoffset=50, base, graphics_level=2, retain=2, uvalue='draw',/EXPOSE_EVENTS, /BUTTON_EVENTS)

slider=CW_FSLIDER(base, xsize=100, ysize=0, uvalue='slide', /drag, value=scale, min=0.01, max=10.0)


WIDGET_CONTROL, BASE, /REALIZE
WIN=OBJ_NEW()
WIDGET_CONTROL, DRAW, GET_VALUE=VALUE
WIN=VALUE

WIDGET_CONTROL, Draw, GET_VALUE=win

model=OBJ_NEW('IDLgrModel')
;model= OBJ_NEW('IDLitVisualization')
;manip= OBJ_NEW('IDLitManipulator')
;itwin->AddWindowEventObserver, manip
;manip->OnMouseDown, itwin, p,q, ButtonMask, Modifiers, NumClicks


VIEW=OBJ_NEW('IDLgrView', COLOR=[0,0,0], VIEW=1/scale*[0,0,xs,ys],zclip=[2000,-1000], eye=2000+0.1)







rtr=xs/2
tx=xs/2
ty=ys/2


Track1 = OBJ_NEW('TrackBall', [tx,ty], rtr)
Track2 = OBJ_NEW('TrackBall', [tx,ty], rtr)

;goto, q

n=256l

nz=64


;bn=fltarr(n,n)
;
;x1=replicate(1.,n)##indgen(n)
;y1=indgen(n)##replicate(1.,n)

x=400
y=400;+32
dx=256
dy=256

;i=4

hmi=rfits('D:\work\hi-c\hmi.M_45s.20120711_183000_TAI.2.magnetogram.fits', header)
crpix1=sxpar(header,'crpix1')
crpix2=sxpar(header,'crpix2')
crota2=sxpar(header, 'crota2')
temp=rot(hmi,crota2,1/1.2,crpix1,crpix2,/int, cubic=-0.5)
temp=extrac(rot(temp,0,1024/176.,400*4+22+90,300*4+150+91.5+12,/int, cubic=-0.5),1536,1536,1024,1024)

;bn=extrac(temp> (-1e3),x-dx/2,y-dy/2,dx,dy)
bn=rebin(temp,256,256)

;bn=1e3*exp(-(x1-100)^2./3^2.-(y1-100)^2./3^2.)-1e3*exp(-(x1-150)^2./3^2.-(y1-150)^2./3^2.);+1*exp(-(x1-20)^2./3^2.-(y1-40)^2./3^2.)

t1=systime(1)

phi=fltarr(n,n,nz)
;bn=fltarr(n,n)
;temp=fltarr(20,20,20)
;bn[8,9]=10
;bn[11,9]=-10
;bn=reform(q[2,*,*,0])

bnpos=bn gt 30
regpos=label_region(bnpos, /all_nei)
nregpos=max(regpos)

bnneg=bn lt -30
regneg=label_region(bnneg, /all_nei)
nregneg=max(regneg)

t=where(abs(bn) le 30)
bn[t]=0

;bn/=100

xind=reform(replicate(1.,2*n*nz)##(findgen(2*n)-n),2*n,2*n,nz)
yind=transpose(xind, [1,0,2])
;zind=transpose(xind, [1,2,0])+1/sqrt(2*!pi)

zind=reform(findgen(nz)##replicate(1.,4*n^2),2*n,2*n,nz)+1/sqrt(2*!pi)

e=1/sqrt(xind^2+yind^2+zind^2)/(2*!pi)


for i=1, nregpos do begin
q=regpos eq i
tw=where(q,ntw)
if ntw gt 6 then begin
arrtw=array_indices(regpos,tw)

xtw=reform(arrtw[0,*])
ytw=reform(arrtw[1,*])

xctw=(max(xtw)+min(xtw))*0.5
yctw=(max(ytw)+min(ytw))*0.5

tbn=(bn*q)[min(xtw):max(xtw),min(ytw):max(ytw)]
if n_elements(dim) eq 2 then begin
dim=size(tbn, /dim)

tbn=reform(tbn, dim[0],dim[1],1)
qqq=convol(e, tbn, /center)
qqq=extrac(qqq, n-xctw,n-yctw,0,n,n,nz)
phi+=qqq
endif
endif
endfor


for i=1, nregneg do begin
q=regneg eq i
tw=where(q,ntw)

if ntw gt 6 then begin
arrtw=array_indices(regneg,tw)

xtw=reform(arrtw[0,*])
ytw=reform(arrtw[1,*])

xctw=(max(xtw)+min(xtw))*0.5
yctw=(max(ytw)+min(ytw))*0.5

tbn=(bn*q)[min(xtw):max(xtw),min(ytw):max(ytw)]
dim=size(tbn, /dim)
if n_elements(dim) eq 2 then begin
tbn=reform(tbn, dim[0],dim[1],1)
qqq=convol(e, tbn, /center)
qqq=extrac(qqq, n-xctw,n-yctw,0,n,n,nz)
phi+=qqq
endif
endif
endfor




t2=systime(1)

print, t2-t1


t3=systime(1)

bx=(shift(phi,-1,0,0)-shift(phi,1,0,0))/2.
by=(shift(phi,0,-1,0)-shift(phi,0,1,0))/2.
bz=(shift(phi,0,0,-1)-shift(phi,0,0,1))/2.

bx[0,*,*]=(-3.*phi[0,*,*] + 4.*phi[1,*,*] - phi[2,*,*])/2.
by[*,0,*]=(-3.*phi[*,0,*] + 4.*phi[*,1,*] - phi[*,2,*])/2.
bz[*,*,0]=(-3.*phi[*,*,0] + 4.*phi[*,*,1] - phi[*,*,2])/2.

bx[n-1,*,*]=(3.*phi[n-1,*,*] - 4.*phi[n-2,*,*] + phi[n-3,*,*])/2.
by[*,n-1,*]=(3.*phi[*,n-1,*] - 4.*phi[*,n-2,*] + phi[*,n-3,*])/2.
bz[*,*,nz-1]=(3.*phi[*,*,nz-1] - 4.*phi[*,*,nz-2] + phi[*,*,nz-3])/2.

q1=fltarr(3,n,n,nz)
q1[0,*,*,*]=bx
q1[1,*,*,*]=by
q1[2,*,*,*]=bz

q:
m=5

xind=replicate(1,n/m)##(m*indgen(n/m))
yind=transpose(xind)
xind=reform(xind,(n/m)^2.)
yind=reform(yind,(n/m)^2.)

;tw=where(bnneg, ntw)
;arrtw=array_indices(bnneg,tw)
;
;xtw=reform(arrtw[0,*])
;ytw=reform(arrtw[1,*])


;vector_field, q, verts, conn
PARTICLE_TRACE, q1, transpose([[xind],[yind],[intarr((n/m)^2)+5]]), Verts, Conn;, integration=1
;q:
fline=OBJ_NEW('IDLgrPolyline', verts, polylines=conn, COLOR=[0,255,0])
t4=systime(1)

print, t4-t3
;PARTICLE_TRACE, q1, transpose([[xind],[yind],[intarr(400)]]), Verts, Conn
;fline1=OBJ_NEW('IDLgrPolyline', verts, polylines=conn, COLOR=[255,0,0])

surf = OBJ_NEW('IDLgrSurface', [[0,0],[0,0]],[[0,n],[0,n]],[[0,0],[n,n]], color=[255,255,255],style=2)

;temp=rfits('D:\TES_20090624_out\Fe132_20090624.090535.fits')
;temp=unsharp_mask(temp, amount=3, radius=3)^0.5
image=OBJ_NEW('IDLgrImage', bytscl(bn))
surf->SetProperty, texture_map=image

model->add,surf
model->add, fline

VIEW->ADD, model

;t3d, matrix=matrix,translate=[-10,-10,0]
t3d, matrix=matrix,scale=[10,10,10]


;t3d, matrix, matrix=matrix,rotate=[-90,-90,0]


model->setproperty, transform=matrix

state={btndown:0b, scale:scale, dx:dx, dy:dy,slider:slider, track1:track1, track2:track2,model:model, Draw: Draw, Win: Win, View: View, rtr:rtr}
WIDGET_CONTROL, Base, SET_UVALUE=State, /NO_COPY
win->draw, view

XMANAGER, 'model', Base, /NO_BLOCK

end

