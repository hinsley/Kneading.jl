from PIL import Image, ImageDraw, ImageFont
im=Image.new('RGB',(1400,890),'white'); d=ImageDraw.Draw(im)
font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
f=lambda n: ImageFont.truetype(font,n)
ink='#193549'; red='#c03351'; gray='#8796a5'
d.text((40,22),'One corner has matrix B; the other three have matrix A',font=f(28),fill=ink)
d.text((40,65),'Each panel applies the current rule to one choice of the shared endpoint p.',font=f(21),fill=ink)
def panel(x,y,size,title,tri=None,segs=(),fillloop=False,caption=''):
    pts={'bl':(x,y+size),'br':(x+size,y+size),'tr':(x+size,y),'tl':(x,y)}
    pts.update(L=(x,y+size/2),D=(x+size/2,y+size),C=(x+size/2,y+size/2))
    d.text((x-22,y-46),title,font=f(21),fill=ink)
    if tri:
        d.polygon([pts[t] for t in tri],fill='#e9eff5')
        d.line([pts[t] for t in tri]+[pts[tri[0]]],fill='#b6c3cf',width=2)
    d.line([pts[t] for t in ['bl','br','tr','tl','bl']],fill=gray,width=2)
    if fillloop: d.polygon([pts[t] for t in ['L','D','C']],fill='#fce5eb')
    for a,b in segs: d.line([pts[a],pts[b]],fill=red,width=5)
    for k in ['bl','br','tr','tl']:
        px,py=pts[k];d.ellipse((px-5,py-5,px+5,py+5),fill=ink)
        d.text((px+(-27 if k.endswith('l') else 10),py+(-27 if k.startswith('t') else 7)),'B' if k=='bl' else 'A',font=f(22),fill=ink)
    for k in sorted(set(t for s in segs for t in s)):
        px,py=pts[k];d.ellipse((px-4,py-4,px+4,py+4),fill='white',outline=red,width=2)
        dx,dy={'L':(-25,-14),'D':(-7,9),'C':(10,-25)}[k];d.text((px+dx,py+dy),k,font=f(20),fill=red)
    d.text((x-22,y+size+51),caption,font=f(19),fill=ink)
panel(65,180,210,'1. p = lower left',['bl','br','tl'],[('L','D')],caption='Retained: L, D → draw LD')
panel(415,180,210,'2. p = lower right',['br','bl','tr'],[('D','C')],caption='Retained: D, C → draw DC')
panel(765,180,210,'3. p = upper left',['tl','bl','tr'],[('L','C')],caption='Retained: L, C → draw LC')
panel(1115,180,210,'4. p = upper right',['tr','tl','br'],[],caption='All A → draw nothing')
panel(270,565,230,'Combined output',segs=[('L','D'),('D','C'),('L','C')],fillloop=True,caption='Three segments form an extra triangle.')
panel(900,565,230,'Diagonal-free fix for this case',segs=[('L','D')],caption='Two changed side edges → join their midpoints.')
im.save('/home/carter/Documents/paper-review/contour-triangle-explanation.png')
