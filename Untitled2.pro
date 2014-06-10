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



VIEW=OBJ_NEW('IDLgrView', COLOR=[0,0,0], VIEW=1/scale*[0,0,xs,ys],zclip=[2000,-1000], eye=2000+0.1)


rtr=xs/2
tx=xs/2
ty=ys/2


Track1 = OBJ_NEW('TrackBall', [tx,ty], rtr)
Track2 = OBJ_NEW('TrackBall', [tx,ty], rtr)


fline=OBJ_NEW('IDLgrPolyline', verts, polylines=conn, COLOR=[0,255,0])

myimage = OBJ_NEW('IDLgrImage')
myimage->SetProperty, data=bytscl(data[*,*,10])


X=fltarr(2,2)
Y=fltarr(2,2)
Z=fltarr(2,2)

x=[[0,256],[0,256]]
y=[[0,0],[256,256]]
z=[[0,0],[0,0]]

;X=[-1,1,1,-1, -1]*6.5
;Y=[0,0,0,0,0]
;Z=[-1,-1,1,1,-1]*0.3


panel = OBJ_NEW('IDLgrSurface', Z,X,Y,style=2,TEXTURE_COORD = [[0,0], [1,0], [0,1], [1,1]], texture_interp=1, color=[255,255,255],TEXTURE_MAP=myimage)

model->add, fline
model->add, panel

VIEW->ADD, model

;t3d, matrix=matrix,translate=[-10,-10,0]
t3d, matrix=matrix,scale=[10,10,50]


;t3d, matrix, matrix=matrix,rotate=[-90,-90,0]


model->setproperty, transform=matrix

state={btndown:0b, scale:scale, dx:dx, dy:dy,slider:slider, track1:track1, track2:track2,model:model, Draw: Draw, Win: Win, View: View, rtr:rtr}
WIDGET_CONTROL, Base, SET_UVALUE=State, /NO_COPY
win->draw, view

XMANAGER, 'model', Base, /NO_BLOCK

end

