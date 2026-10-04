"""R03 structural RC glider. Units mm, densities kg/m3, masses g.
Design assumptions are explicit. No measured flight/strength certification.
Run using the bundled Python with workspace .cad_deps available.
"""
from pathlib import Path
import sys, math, json, csv
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.cad_deps'))
import cadquery as cq
import numpy as np
from scipy.interpolate import PchipInterpolator
OUT=ROOT/'output/Glider_R03'
for d in ['parts_step','cutting','analysis','solidworks','preview']:(OUT/d).mkdir(parents=True,exist_ok=True)
PARTS=[]; SHAPES={}; CUTS=[]
rho={'balsa':160,'dense_balsa':220,'plywood':600,'nylon':1150,'steel':7850,'film':1000,'equipment':1000,'adhesive':1100}
colors={'balsa':(0.88,.72,.43),'dense_balsa':(.65,.40,.18),'plywood':(.75,.5,.27),'nylon':(.35,.38,.45),'steel':(.55,.58,.62),'film':(.30,.65,.85),'equipment':(.18,.23,.30),'adhesive':(.55,.55,.5)}
def box(x0,x1,y0,y1,z0,z1):
    return cq.Solid.makeBox(x1-x0,y1-y0,z1-z0,cq.Vector(x0,y0,z0))
def poly(pts):return cq.Wire.makePolygon([cq.Vector(*p) for p in pts],close=True)
def loft(wires):return cq.Solid.makeLoft(wires,ruled=True)
def extrude(w,vec):return cq.Solid.extrudeLinear(w,[],cq.Vector(*vec))
def rod(a,b,r):
    a,b=cq.Vector(*a),cq.Vector(*b);v=b-a
    return cq.Solid.makeCylinder(r,v.Length,a,v.normalized())
def tube(a,b,ro,ri):return rod(a,b,ro).cut(rod(a,b,ri))
def clean(s):
    try:return s.clean()
    except:return s
def add(name,s,material,group,mass=None,note='',cut=None):
    s=clean(s)
    if not s.isValid() or s.Volume()<1e-5:raise ValueError(f'Invalid/empty: {name}')
    v=s.Volume();c=s.Center();bb=s.BoundingBox()
    m=v*rho[material]*1e-6 if mass is None else mass
    entry=dict(name=name,material=material,group=group,volume_mm3=v,mass_g=m,mass_basis='density_assumption' if mass is None else 'specified_design_allowance',cg_mm=[c.x,c.y,c.z],bounds_mm=[bb.xmin,bb.ymin,bb.zmin,bb.xmax,bb.ymax,bb.zmax],note=note)
    PARTS.append(entry);SHAPES[name]=s
    cq.exporters.export(s,str(OUT/'parts_step'/f'{name}.step'))
    if cut:CUTS.append(dict(name=name,**cut))
    print(name,round(m,2),flush=True)
    return s
def cut_xy(s,name,thickness,material,z):
    section=cq.Workplane('XY').newObject([s]).section(z)
    path=OUT/'cutting'/f'{name}.dxf'
    cq.exporters.export(section,str(path))
    CUTS.append(dict(name=name,material=material,thickness=thickness,source=str(path.name)))

dat=np.loadtxt(ROOT/'source/sd7037.dat',skiprows=1);k=np.argmin(dat[:,0])
upper=PchipInterpolator(np.r_[0,dat[:k+1,0][::-1]],np.r_[0,dat[:k+1,1][::-1]])
lower=PchipInterpolator(np.r_[0,dat[k+1:,0]],np.r_[0,dat[k+1:,1]])
ang=math.radians(2);dh=math.tan(math.radians(4))
def chord(y):return 170-40*abs(y)/850
def pt(u,v,y):return (220+u*math.cos(ang)+v*math.sin(ang),y,27+abs(y)*dh-u*math.sin(ang)+v*math.cos(ang))
def foilwire(y,x0=0,x1=None,topoff=0,botoff=0):
    c=chord(y);x1=c if x1 is None else x1
    xs=np.linspace(x0,x1,35)
    up=[cq.Vector(*pt(float(x),float(c*upper(x/c)-topoff),y)) for x in xs]
    lo=[cq.Vector(*pt(float(x),float(c*lower(x/c)+botoff),y)) for x in xs]
    edges=[cq.Edge.makeSpline(up),cq.Edge.makeSpline(list(reversed(lo)))]
    if (up[-1]-lo[-1]).Length>1e-7:edges.insert(1,cq.Edge.makeLine(up[-1],lo[-1]))
    if (up[0]-lo[0]).Length>1e-7:edges.append(cq.Edge.makeLine(lo[0],up[0]))
    return cq.Wire.assembleEdges(edges)
def bandwire(y,x0,x1,mode,t):
    c=chord(y);xs=np.linspace(x0,x1,18)
    if mode=='top': f=lambda x:c*upper(x/c); g=lambda x:f(x)-t
    else: g=lambda x:c*lower(x/c); f=lambda x:g(x)+t
    up=[cq.Vector(*pt(float(x),float(f(x)),y)) for x in xs]
    lo=[cq.Vector(*pt(float(x),float(g(x)),y)) for x in xs]
    return cq.Wire.assembleEdges([cq.Edge.makeSpline(up),cq.Edge.makeLine(up[-1],lo[-1]),cq.Edge.makeSpline(list(reversed(lo))),cq.Edge.makeLine(lo[0],up[0])])
def webwire(y,x0,x1,skin=1,cap=3):
    c=chord(y);xs=np.linspace(x0,x1,4)
    pts=[pt(float(x),float(c*upper(x/c)-skin-cap),y) for x in xs]
    pts += [pt(float(x),float(c*lower(x/c)+skin+cap),y) for x in xs[::-1]]
    return poly(pts)

# Full-span one-piece removable wing; a central plywood V-joiner carries bending.
print('BUILD WING',flush=True)
wing_env={};wing_members={};rib_meta=[]
bolts=[(245,0),(360,0)]
bolt_cutters=[rod((x,y,-5),(x,y,60),1.7) for x,y in bolts]
for side,label in [(1,'R'),(-1,'L')]:
    env=loft([foilwire(side*y) for y in [0,850]])
    wing_env[side]=env;members=[]
    for which in ['top','bottom']:
        # D-box skin, 1 mm vertical thickness, sand to the supplied section templates.
        s=loft([bandwire(side*y,4,.3*chord(y)+10,which,1) for y in [0,850]])
        for q in bolt_cutters:s=s.cut(q)
        members.append(add(f'W_{label}_Dskin_{which}',s,'balsa','wing',note='1 mm balsa formed skin; chordwise width varies; grain spanwise'))
        # Caps follow the inside skin and are 20 x 3 mm nominal.
        wires=[]
        for y in [0,850]:
            w=bandwire(side*y,.3*chord(y)-10,.3*chord(y)+10,which,3)
            dz=-1 if which=='top' else 1
            w=w.translate((math.sin(ang)*dz,0,math.cos(ang)*dz));wires.append(w)
        members.append(add(f'W_{label}_Cap_{which}',loft(wires),'dense_balsa','wing',note='20 mm wide x 3 mm nominal; dense straight-grain balsa; spanwise grain'))
    for face,sgn in [('front',-1),('rear',1)]:
        members.append(add(f'W_{label}_Web_{face}',loft([webwire(side*y,.3*chord(y)+sgn*9.25-.75,.3*chord(y)+sgn*9.25+.75) for y in [0,850]]),'balsa','wing',note='1.5 mm shear web; grain through depth; fit between spar caps'))
    le=loft([foilwire(side*y,0,4) for y in [0,850]])
    members.append(add(f'W_{label}_Leading_Edge',le,'dense_balsa','wing',note='Carve/sand a 4 mm chordwise leading-edge blank to templates'))
    te=loft([foilwire(side*y,.90*chord(y),chord(y)) for y in [0,850]])
    members.append(add(f'W_{label}_Trailing_Edge',te,'dense_balsa','wing',note='Sand tapered trailing strip; retain 0.3 mm practical edge and fair covering'))
    # Root full sheeting up to y=50, aft of D-box, top and bottom.
    for which in ['top','bottom']:
        s=loft([bandwire(side*y,.3*chord(y)+10,.9*chord(y),which,1) for y in [0,50]])
        for q in bolt_cutters:s=s.cut(q)
        members.append(add(f'W_{label}_Root_Sheet_{which}',s,'balsa','wing',note='1 mm root sheeting; supports wing attachment area'))
    wing_members[side]=members

# Two 4 mm birch-ply V braces outside cap edges; no carbon wing spar.
joiners=[]
for face,sgn in [('front',-1),('rear',1)]:
    halves=[]
    for side in [1,-1]:
        wires=[]
        for y in [0,150]:
            c=chord(y);xc=.3*c+sgn*12
            xs=np.linspace(xc-2,xc+2,4)
            pts=[pt(float(x),float(c*upper(x/c)-1),side*y) for x in xs]
            pts += [pt(float(x),float(c*lower(x/c)+1),side*y) for x in xs[::-1]]
            wires.append(poly(pts))
        halves.append(loft(wires))
    joiners.append(add(f'W_Center_Joiner_{face}',halves[0].fuse(halves[1]),'plywood','wing',note='4 mm birch aircraft plywood; paired 300 mm dihedral braces; bond continuously to caps and webs'))

# Bolt compression blocks cut to actual wing outer surface.
wing_full=wing_env[1].fuse(wing_env[-1]);mountblocks=[]
for i,(x,y) in enumerate(bolts):
    s=wing_full.intersect(box(x-9,x+9,-12,12,10,50))
    for member in sum(wing_members.values(),[]):
        s=s.cut(member)
    for q in bolt_cutters:s=s.cut(q)
    mountblocks.append(add(f'W_Bolt_Block_{i+1}',s,'plywood','wing',note='Laminated plywood compression block, drilled 3.4 mm'))

for side,label in [(1,'R'),(-1,'L')]:
    for i,y in enumerate([1.5,50,110,180,250,320,390,460,530,600,670,740,810,849.25]):
        t=3 if i==0 else (2 if i<3 else 1.5);ya=max(0,y-t/2);yb=min(850,y+t/2)
        a,b=sorted([side*ya,side*yb]);c=chord(y);sp=.3*c
        # Constant-thickness laser-cut rib, chord-plane section. Planar booleans
        # avoid slow nearly-coincident curved-solid booleans and are manufacturable.
        xs=np.unique(np.r_[np.linspace(4,.9*c,110),sp-14,sp-10,sp+10,sp+14])
        xs=xs[(xs>=4)&(xs<=.9*c)]
        def skin(x):return 1 if (x<=sp+10 or y<=50) else 0
        pts=[(float(x),float(c*upper(x/c)-skin(x)),0) for x in xs]
        pts += [(float(x),float(c*lower(x/c)+skin(x)),0) for x in xs[::-1]]
        r=extrude(poly(pts),(0,0,t))
        # Cap seats and vertical shear webs. Small 0.08 mm assembly allowance.
        for xx in np.linspace(sp-10,sp+10,12):
            width=20/11+.08
            r=r.cut(box(xx-width/2,xx+width/2,float(c*upper(xx/c)-4)-.08,30,-1,t+1))
            r=r.cut(box(xx-width/2,xx+width/2,-20,float(c*lower(xx/c)+4)+.08,-1,t+1))
        for xx in [sp-9.25,sp+9.25]:r=r.cut(box(xx-.83,xx+.83,-20,30,-1,t+1))
        if y<=150:
            for xx in [sp-12,sp+12]:r=r.cut(box(xx-2.08,xx+2.08,-20,30,-1,t+1))
        if y<12:
            for xx in [25,140]:r=r.cut(box(xx-9.1,xx+9.1,-20,30,-1,t+1))
        # Aft lightening holes only, leave nose/D-box rib webs intact.
        c=chord(y)
        for frac in [.52,.72]:
            cx,cz=c*frac,float(c*(upper(frac)+lower(frac))/2)
            depth=float(c*(upper(frac)-lower(frac)))
            rad=max(1.0,min(3.0,(depth-4)/2))
            r=r.cut(rod((cx,cz,-1),(cx,cz,t+1),rad))
        pl=cq.Plane(origin=(220,side*y+t/2,27+y*dh),xDir=(math.cos(ang),0,-math.sin(ang)),normal=(0,-1,0))
        r=r.transformShape(pl.rG)
        name=f'W_{label}_Rib_{i:02d}';mat='plywood' if i<3 else 'balsa'
        add(name,r,mat,'wing',note=f'Rib station abs(Y)={y} mm; sheet {t} mm; grain chordwise')
        # Take center section in XZ, then flatten the section for cutting.
        plane=cq.Plane(origin=(0,side*y,0),xDir=(1,0,0),normal=(0,-1,0))
        sec=cq.Workplane(plane).newObject([r]).section(0)
        shapes=[w.transformShape(plane.fG) for w in sec.vals()]
        cq.exporters.export(cq.Compound.makeCompound(shapes),str(OUT/'cutting'/f'{name}.dxf'))
        CUTS.append(dict(name=name,material=mat,thickness=t,source=f'{name}.dxf'))
        rib_meta.append(dict(name=name,y=side*y,chord=c,thickness=t))

# Tail: straight physical hinge axes, two independent elevator halves.
print('BUILD TAIL',flush=True)
mac=2*170/3*(1+(130/170)+(130/170)**2)/(1+130/170)
xac=220+.25*mac*math.cos(ang);htac=xac+540;hingex=htac+39.25
def tail_plan(y):c=105-15*abs(y)/200;return htac-c*.25,htac+c*.75
le0,te0=tail_plan(0);let,tet=tail_plan(200)
stabwire=poly([(le0,0,16),(let,200,16),(hingex-.25,200,16),(hingex-.25,-200,16),(let,-200,16)])
stab=extrude(stabwire,(0,0,4))
for side in [-1,1]:
    for inner,outer in [(30,74),(94,138),(158,192)]:
        y0,y1=sorted((side*inner,side*outer))
        # Use the outermost leading edge to retain at least 8 mm of wood.
        le,_=tail_plan(outer)
        stab=stab.cut(box(le+8,hingex-8,y0,y1,15,21))
# Guard against open tip bays and asymmetric lightening cuts.
if len(stab.Solids())!=1:raise ValueError('Horizontal stabilizer is disconnected')
if stab.cut(stab.mirror('XZ')).Volume()>1e-5:raise ValueError('Asymmetric horizontal stabilizer')
# Keep central mounting region solid, fasteners at x=790 and 815.
tailbolts=[(790,6),(815,-6)]
for x,y in tailbolts:stab=stab.cut(rod((x,y,0),(x,y,30),1.1))
add('T_Horizontal_Stabilizer',stab,'balsa','tail',note='4 mm sheet, grain spanwise; open bays covered with film')
cut_xy(stab,'T_Horizontal_Stabilizer',4,'balsa',18)
elevators=[]
for side,label in [(1,'R'),(-1,'L')]:
    y0,y1=side*24,side*200
    _,tea=tail_plan(y0);_,teb=tail_plan(y1)
    w=poly([(hingex+.25,y0,16),(tea,y0,16),(teb,y1,16),(hingex+.25,y1,16)])
    s=extrude(w,(0,0,4))
    # Bottom bevel: tape hinge on z=20, adequate for +/-15 degrees.
    bevel=extrude(poly([(hingex+.24,min(y0,y1)-1,15),(hingex+4,min(y0,y1)-1,15),(hingex+.24,min(y0,y1)-1,19.8)]),(0,abs(y1-y0)+2,0))
    s=s.cut(bevel);name=f'T_Elevator_{label}'
    add(name,s,'balsa','tail',note='4 mm sheet; lower leading edge bevel; tape hinge on top; neutral geometry')
    cut_xy(s,name,4,'balsa',19.9);elevators.append(s)
    for j,y in enumerate([45,110,175]):
        yy=side*y
        add(f'H_Elevator_{label}_Tape_{j}',box(hingex-5,hingex+5,yy-8,yy+8,20,20.06),'film','hardware',note='Flexible 0.06 mm hinge tape; wrap/seal per assembly guide')

vh=190;vx=849.5
fin=extrude(poly([(780,-2,20),(804,-2,210),(vx-.25,-2,210),(vx-.25,-2,20)]),(0,4,0))
for z0 in [40,95,150]:
    le=780+24*(z0-20)/190
    fin=fin.cut(box(le+9,vx-8,-3,3,z0,z0+38))
add('T_Vertical_Fin',fin,'balsa','tail',note='4 mm sheet; straight hinge at X=849.5; grain vertical')
rudder=extrude(poly([(vx+.25,-2,25),(900-6*5/190,-2,25),(894,-2,210),(vx+.25,-2,210)]),(0,4,0))
bevel=extrude(poly([(vx+.24,2.01,24),(vx+4,2.01,24),(vx+.24,-1.8,24)]),(0,0,188))
rudder=rudder.cut(bevel)
add('T_Rudder',rudder,'balsa','tail',note='4 mm sheet; one-face bevel; lower 5 mm clearance; +/-25 deg design travel')
for j,z in enumerate([45,110,180]):add(f'H_Rudder_Tape_{j}',box(vx-5,vx+5,-2.06,-2,z-8,z+8),'film','hardware')
for name,s in [('T_Vertical_Fin',fin),('T_Rudder',rudder)]:
    pl=cq.Plane(origin=(0,0,0),xDir=(1,0,0),normal=(0,-1,0))
    sec=cq.Workplane(pl).newObject([s]).section(0)
    flat=cq.Compound.makeCompound([v.transformShape(pl.fG) for v in sec.vals()])
    cq.exporters.export(cq.Workplane('XY').newObject([flat]),str(OUT/'cutting'/f'{name}.dxf'))
    CUTS.append(dict(name=name,material='balsa',thickness=4,source=f'{name}.dxf'))

# Buildable four-panel tapered balsa fuselage, not a solid aerodynamic envelope.
print('BUILD FUSELAGE',flush=True)
stations=[(40,25,24),(210,25,24),(395,25,20),(500,17,16),(750,10,9),(900,8,8)]
def dims(x):return tuple(float(np.interp(x,[s[0] for s in stations],[s[i] for s in stations])) for i in [1,2])
def fuseloft(inner=False):
    wires=[]
    for x,w,h in stations:
        if inner:w-=2;h-=2
        wires.append(poly([(x,-w,-h),(x,w,-h),(x,w,h),(x,-w,h)]))
    return loft(wires)
outer=fuseloft();inner=fuseloft(True)
sidewalls=[]
for side,label in [(1,'R'),(-1,'L')]:
    wires=[]
    for x,w,h in stations:
        ya,yb=sorted([side*w,side*(w-2)])
        wires.append(poly([(x,ya,-h),(x,yb,-h),(x,yb,h),(x,ya,h)]))
    s=loft(wires).cut(wing_full)
    sidewalls.append(add(f'F_Side_{label}',s,'balsa','fuselage',note='2 mm sheet; scarf segments if sheet stock shorter than 900 mm; grain longitudinal'))
for which in ['Bottom','Top']:
    wires=[]
    for x,w,h in stations:
        z0,z1=(-h,-h+2) if which=='Bottom' else (h-2,h)
        wires.append(poly([(x,-w+2,z0),(x,w-2,z0),(x,w-2,z1),(x,-w+2,z1)]))
    s=loft(wires)
    if which=='Top':
        # Separate hatch and wing saddle from the roof.
        s=s.cut(box(65,205,-30,30,10,40)).cut(box(215,390,-30,30,10,45)).cut(wing_full)
    add(f'F_{which}',s,'balsa','fuselage',note='2 mm balsa; top openings for battery hatch and removable wing')
# Rounded nose block. Back is open to pod; interior removed for ballast pocket.
nose=cq.Workplane('XY').box(42,50,48,centered=(False,True,True)).edges('|X').fillet(12).faces('<X').edges().fillet(15).val()
nose=nose.cut(box(20,43,-15,15,-14,14))
add('F_Rounded_Nose',nose,'balsa','fuselage',note='Carved balsa nose; 15 mm nose fillet; enclosed ballast cavity')

# Frame and tray clearances include Bowden routing.
formers=[]
for i,x in enumerate([42,65,205,215,300,390,500,650,750,825,898]):
    w,h=dims(x);t=2 if x<400 else 1.5
    f=box(x-t/2,x+t/2,-w+2,w-2,-h+2,h-2)
    margin=4
    if w>8 and h>8:f=f.cut(box(x-3,x+3,-w+2+margin,w-2-margin,-h+2+margin,h-2-margin))
    f=f.cut(wing_full)
    formers.append(add(f'F_Former_{i:02d}',f,'plywood' if x in [205,215,300,390,825] else 'balsa','fuselage',note=f'X={x}; t={t}; former and equipment support'))
    pl=cq.Plane(origin=(x,0,0),xDir=(0,1,0),normal=(1,0,0));sec=cq.Workplane(pl).newObject([f]).section(0)
    name=f'F_Former_{i:02d}';cq.exporters.export(cq.Compound.makeCompound([v.transformShape(pl.fG) for v in sec.vals()]),str(OUT/'cutting'/f'{name}.dxf'))
    CUTS.append(dict(name=name,material='plywood' if x in [205,215,300,390,825] else 'balsa',thickness=t,source=f'{name}.dxf'))
tray=box(68,200,-23,23,-16,-14)
for x in [85,115,145,175]:tray=tray.cut(box(x,x+3,-19,19,-17,-13))
add('F_Battery_Tray',tray,'plywood','fuselage',note='2 mm ply; side edges bonded to inner walls at Y=+/-23; four strap slots; adjustable battery station')
cut_xy(tray,'F_Battery_Tray',2,'plywood',-15)
hatch=box(65,205,-23,23,22,24).cut(rod((198,0,21),(198,0,25),1.1))
add('F_Battery_Hatch',hatch,'balsa','fuselage',note='Removable 2 mm lid, front tongue and rear nylon screw')
add('F_Hatch_Tongue',box(58,78,-8,8,20,22),'plywood','fuselage')
for x in [62,198]:
    ledge=box(x-3,x+3,-23,23,18,20 if x==62 else 22)
    if x==198:ledge=ledge.cut(rod((x,0,17),(x,0,23),.8))
    add(f'F_Hatch_Ledge_{x}',ledge,'plywood' if x==198 else 'balsa','fuselage',note='Rear ledge has M2 pilot hole; front ledge supports hatch tongue')
add('H_Hatch_Screw',rod((198,0,18),(198,0,24),1).fuse(rod((198,0,24),(198,0,25.5),2)),'nylon','hardware',note='M2 removable hatch screw with bearing head; matching tapped plywood pilot')

# Removable wing: paired nylon M3 screws, seat/doublers and captive nuts.
for i,(x,y) in enumerate(bolts):
    seat=box(x-12,x+12,-23,23,14,29).cut(wing_full).cut(rod((x,0,5),(x,0,60),1.7))
    # Only the lower saddle is structural; subtracting the airfoil from a tall
    # box otherwise leaves a disconnected chunk above the aft wing surface.
    seat=min(seat.Solids(),key=lambda solid:solid.Center().z)
    add(f'F_Wing_Seat_{i}',seat,'plywood','fuselage',note='Laminated ply saddle; upper face matches wing; M3 clearance hole')
    # Recess the head into a horizontal bearing surface instead of floating
    # above the airfoil; pocket is cut in every affected wing member below.
    u=(x-220)/math.cos(ang)
    headz=pt(u,float(170*upper(u/170)),0)[2]-.6
    pocket=rod((x,0,headz),(x,0,60),3.1)
    for entry in PARTS:
        if entry['group']!='wing':continue
        n=entry['name'];original=SHAPES[n];bb=original.BoundingBox()
        if bb.xmax<x-3.1 or bb.xmin>x+3.1 or bb.ymin>3.1 or bb.ymax< -3.1 or bb.zmax<headz:continue
        new=original.cut(pocket)
        if new.Volume()<original.Volume()-1e-5:
            SHAPES[n]=clean(new);entry['volume_mm3']=new.Volume();entry['mass_g']=new.Volume()*rho[entry['material']]*1e-6
            ctr=new.Center();entry['cg_mm']=[ctr.x,ctr.y,ctr.z]
            cq.exporters.export(SHAPES[n],str(OUT/'parts_step'/f'{n}.step'))
    add(f'H_Wing_Screw_{i}',rod((x,0,10),(x,0,headz),1.5).fuse(rod((x,0,headz),(x,0,headz+2),3)),'nylon','hardware',note='M3 nylon screw; recessed head seated on machined flat; trim length after fitting')
    nut=box(x-3,x+3,-3,3,11,14).cut(rod((x,0,10),(x,0,15),1.3))
    add(f'H_Wing_Captive_Nut_{i}',nut,'steel','hardware',note='M3 captive nut, internal only')

# Tail saddle, two screw holes and glued fin foot; no floating tail components.
tailseat=box(778,825,-10,10,7,16).cut(outer)
for x,y in tailbolts:tailseat=tailseat.cut(rod((x,y,0),(x,y,22),1.1))
add('F_Tail_Saddle',tailseat,'plywood','fuselage',note='Laminated ply saddle, shaped to pod; horizontal tail at z=16')
for i,(x,y) in enumerate(tailbolts):add(f'H_Tail_Screw_{i}',rod((x,y,4),(x,y,20),1).fuse(rod((x,y,20),(x,y,21.5),2)),'nylon','hardware',note='M2 screw; captive threaded insert in saddle')

# Generic equipment envelopes; purchased dimensions must fit or model must be updated.
print('BUILD EQUIPMENT',flush=True)
add('E_Battery_Envelope',box(80,125,-15,15,-13,-1),'equipment','equipment',mass=32,note='Design envelope only: 45x30x12 mm, 32 g; receiver-compatible regulated supply required')
add('E_Receiver_Envelope',box(155,185,-10,10,-13,-3),'equipment','equipment',mass=6,note='Design envelope 30x20x10, 6 g; elevator Y harness and rudder channels')
servo_centers=[(315,-11),(315,11),(355,0)]
tray=box(300,372,-22,22,2,4)
for x,y in servo_centers:tray=tray.cut(box(x-12,x+12,y-6.5,y+6.5,1,5))
add('F_Servo_Tray',tray,'plywood','fuselage',note='2 mm plywood tray; three 24x13 mm servo cutouts')
cut_xy(tray,'F_Servo_Tray',2,'plywood',3)
for side,label in [(1,'R'),(-1,'L')]:
    a,b=sorted([side*20,side*23]);add(f'F_Servo_Tray_Rail_{label}',box(300,372,a,b,0,2),'dense_balsa','fuselage',note='Support rail bonded to inner fuselage side; tray bonds onto z=2 face')
for i,(x,y) in enumerate(servo_centers):
    add(f'E_Servo_{i+1}',box(x-11.5,x+11.5,y-6,y+6,-16,4),'equipment','equipment',mass=9,note='Generic 23x12x20 mm 9 g servo envelope; specify >=0.8 kg.cm at supply voltage')
    arm=rod((x+5,y,4),(x+5,y,6),3).fuse(box(x+4,x+6,y-8,y+8,4,6))
    add(f'H_Servo_Arm_{i+1}',arm,'nylon','hardware',mass=.5,note='8 mm servo arm radius; neutral position')
    # Mounting tabs and screw positions on rails beyond body cutout.
    for xx in [x-14,x+14]:
        tx0,tx1=(xx-2,x-11.5) if xx<x else (x+11.5,xx+2)
        add(f'H_Servo_{i+1}_Tab_{int(xx)}',box(tx0,tx1,y-6,y+6,4,5.2),'nylon','hardware',note='Generic integral mounting lug touches servo body; drill mounting screw pilot on assembly')
        add(f'H_Servo_{i+1}_Screw_{int(xx)}',rod((xx,y,1),(xx,y,6),.8),'steel','hardware')

# Two elevator control horns below the halves, one rudder horn on left.
horns=[]
for side,label in [(1,'R'),(-1,'L')]:
    y=side*36;x=hingex+7
    s=box(x-4,x+5,y-1,y+1,5,17)
    s=s.cut(rod((x,y-2,7),(x,y+2,7),.6))
    add(f'H_Elevator_Horn_{label}',s,'plywood','hardware',note='2 mm ply horn; bond and pin to elevator; control hole ~13 mm below hinge')
    horns.append((x,y,7))
rudhorn=box(vx+3,vx+12,-14,-1.8,42,44).cut(rod((vx+7,-12,41),(vx+7,-12,45),.6))
add('H_Rudder_Horn',rudhorn,'plywood','hardware',note='2 mm ply horn, 12 mm control lever')
horns.append((vx+7,-12,43))
horns=[horns[1],horns[0],horns[2]]
for i,((x,y),end) in enumerate(zip(servo_centers,horns)):
    start=(x+5,y+(-8 if i==0 else 8),5)
    # A supported flexible steel inner wire, route bows through tube; final clevis straight reach.
    direction=np.array(end)-np.array(start)
    a=np.array(start)+direction*.04;b=np.array(start)+direction*.94
    wire_end=np.array(end)-direction/np.linalg.norm(direction)*4.5
    wire=rod(start,wire_end,.4)
    add(f'H_Pushrod_{i+1}',wire,'steel','hardware',note='0.8 mm steel inner wire; supported Bowden routing; adjustable clevises at ends')
    add(f'H_Bowden_Tube_{i+1}',tube(a,b,.9,.5),'nylon','hardware',note='1.8 mm OD / 1.0 mm ID guide tube; bond to formers; cut exit apertures to route')
    ex,ey,ez=end
    clevis=box(ex-5,ex+2,ey-2.2,ey+2.2,ez-2.2,ez+2.2)
    if i<2:
        clevis=clevis.cut(box(ex-4.2,ex+3,ey-1.1,ey+1.1,ez-3,ez+3))
        clevis=clevis.cut(rod((ex,ey-3,ez),(ex,ey+3,ez),.6))
        pin=rod((ex,ey-2.6,ez),(ex,ey+2.6,ez),.55)
    else:
        clevis=clevis.cut(box(ex-4.2,ex+3,ey-3,ey+3,ez-1.1,ez+1.1))
        clevis=clevis.cut(rod((ex,ey,ez-3),(ex,ey,ez+3),.6))
        pin=rod((ex,ey,ez-2.6),(ex,ey,ez+2.6),.55)
    add(f'H_Control_Clevis_{i+1}',clevis,'nylon','hardware',note='Clevis envelope with 2.2 mm horn slot; 1.1 mm pin; retain pin and threaded adjustment')
    add(f'H_Control_Pin_{i+1}',pin,'steel','hardware')
    # Model the necessary exit bore in side/top panels; done below for collision-free routing.

# Remove routing clearances from intersected fuselage panels and formers.
for entry in list(PARTS):
    if entry['group']!='fuselage' or entry['name'].startswith('F_Battery'):continue
    name=entry['name'];s=SHAPES[name]
    if name in ['F_Top','F_Former_09']:
        for x,y in tailbolts:s=s.cut(rod((x,y,0),(x,y,24),1.1))
    for i,((x,y),end) in enumerate(zip(servo_centers,horns)):
        s=s.cut(rod((x+5,y+(-8 if i==0 else 8),5),end,1.1))
    if s.Volume()<SHAPES[name].Volume()-1e-5:
        s=clean(s);assert s.isValid();SHAPES[name]=s
        entry['volume_mm3']=s.Volume();entry['mass_g']=s.Volume()*rho[entry['material']]*1e-6
        c=s.Center();entry['cg_mm']=[c.x,c.y,c.z]
        cq.exporters.export(s,str(OUT/'parts_step'/f'{name}.step'))

# Real mating pockets replace geometric overlaps at bonded wood joints.
joint_cuts={'F_Rounded_Nose':['F_Side_R','F_Side_L','F_Top','F_Bottom','F_Former_00'],
            'F_Former_01':['F_Hatch_Tongue','F_Hatch_Ledge_62'],
            'F_Former_04':['F_Servo_Tray','F_Servo_Tray_Rail_R','F_Servo_Tray_Rail_L'],
            'T_Elevator_R':['H_Elevator_Horn_R'],
            'T_Elevator_L':['H_Elevator_Horn_L'],
            'T_Rudder':['H_Rudder_Horn'],
            'H_Rudder_Tape_0':['H_Rudder_Horn']}
for entry in PARTS:
    name=entry['name'];cutters=list(joint_cuts.get(name,[]))
    if name.startswith('F_Former_'):cutters+=['F_Side_R','F_Side_L','F_Top','F_Bottom']
    if not cutters:continue
    s=SHAPES[name]
    for other in cutters:s=s.cut(SHAPES[other])
    s=clean(s)
    if not s.isValid() or not s.Solids():raise ValueError('Joint pocket failed '+name)
    SHAPES[name]=s;entry['volume_mm3']=s.Volume();entry['mass_g']=s.Volume()*rho[entry['material']]*1e-6
    ctr=s.Center();entry['cg_mm']=[ctr.x,ctr.y,ctr.z]
    cq.exporters.export(s,str(OUT/'parts_step'/f'{name}.step'))
# Update profiles affected by the final mating pockets.
for name in ['T_Elevator_R','T_Elevator_L']:
    sec=cq.Workplane('XY').newObject([SHAPES[name]]).section(19.9)
    cq.exporters.export(sec,str(OUT/'cutting'/f'{name}.dxf'))
pl=cq.Plane(origin=(0,0,0),xDir=(1,0,0),normal=(0,-1,0))
sec=cq.Workplane(pl).newObject([SHAPES['T_Rudder']]).section(0)
cq.exporters.export(cq.Compound.makeCompound([v.transformShape(pl.fG) for v in sec.vals()]),str(OUT/'cutting/T_Rudder.dxf'))
# Refresh former profiles after routing bores; cutting drawings track final CAD.
for i,x in enumerate([42,65,205,215,300,390,500,650,750,825,898]):
    name=f'F_Former_{i:02d}';pl=cq.Plane(origin=(x,0,0),xDir=(0,1,0),normal=(1,0,0))
    sec=cq.Workplane(pl).newObject([SHAPES[name]]).section(0)
    flat=cq.Compound.makeCompound([v.transformShape(pl.fG) for v in sec.vals()])
    cq.exporters.export(cq.Workplane('XY').newObject([flat]),str(OUT/'cutting'/f'{name}.dxf'))

# Film is a mass allowance and external presentation reference, not a solid wing core.
# Area: both wing faces + both tail faces with 10% overlap. Fin and tail open bays sealed.
film_area=2*(.255+.039+.01995)*1.10
film_mass=film_area*25 # g/m2 target film, not a purchased product specification
# Do not model micron shells: separate mass-only bill-of-material entries below.
extras=[dict(name='Covering_film',mass_g=film_mass,cg_mm=[340,0,35],material='film',group='allowance',note='25 g/m2, 10% overlap; distributed-CG approximation'),
        dict(name='Glue_finish_wiring',mass_g=22,cg_mm=[330,0,5],material='adhesive',group='allowance',note='Weigh glue, finish, wire harness and retainers; allowance not CAD volume')]

# Balance via an enclosed ballast pocket. Design mass constrained by 14-20 g/dm2.
mass0=sum(p['mass_g'] for p in PARTS+extras);moment0=sum(p['mass_g']*p['cg_mm'][0] for p in PARTS+extras)
targetcg=220+.30*mac*math.cos(ang)
# Minimum ballast for CG at a fixed forward pocket x=30. Additional distributed margin cannot conceal overweight.
ballast=max(0,(moment0-targetcg*mass0)/(targetcg-30))
if ballast>0:
    # Brass is enclosed ballast, not external structure.
    v=ballast/(8500e-6);h=v/(20*20)
    s=box(20,40,-10,10,-h/2,h/2)
    add('H_Enclosed_Ballast',s,'steel','hardware',mass=ballast,note='Equivalent brass ballast mass; fully enclosed and bonded in nose; do not infer alloy from display material')
mass=sum(p['mass_g'] for p in PARTS+extras)
cg=[sum(p['mass_g']*p['cg_mm'][j] for p in PARTS+extras)/mass for j in range(3)]
for entry in PARTS:
    bb=SHAPES[entry['name']].BoundingBox()
    entry['bounds_mm']=[bb.xmin,bb.ymin,bb.zmin,bb.xmax,bb.ymax,bb.zmax]
assembly=cq.Assembly(name='Glider_R03_Structure')
for p in PARTS:assembly.add(SHAPES[p['name']],name=p['name'],color=cq.Color(*colors[p['material']]))
assembly.export(str(OUT/'Glider_R03_Structure.step'))
(OUT/'analysis/parts.json').write_text(json.dumps(PARTS,indent=2),encoding='utf-8')
(OUT/'analysis/mass_balance.json').write_text(json.dumps(dict(mass_g=mass,cg_mm=cg,target_cg_x_mm=targetcg,ballast_g=ballast,allowances=extras,wing_loading_g_dm2=mass/25.5,material_density_assumptions=rho),indent=2),encoding='utf-8')
(OUT/'analysis/geometry.json').write_text(json.dumps(dict(span=1700,root_chord=170,tip_chord=130,area_mm2=255000,MAC=mac,wing_ac_x=xac,tail_ac_x=htac,tail_arm=540,tail_area_mm2=39000,tail_hinge_x=hingex,rudder_hinge_x=vx,ribs=rib_meta),indent=2),encoding='utf-8')
(OUT/'cutting/cut_manifest.json').write_text(json.dumps(CUTS,indent=2),encoding='utf-8')
with (OUT/'Bill_of_Materials.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['Part','Material','Group','Mass_g','CG_X_mm','CG_Y_mm','CG_Z_mm','Basis','Notes'])
    for p in PARTS+extras:w.writerow([p['name'],p['material'],p['group'],round(p['mass_g'],3),*p['cg_mm'],p.get('mass_basis','allowance'),p.get('note','')])
print('SUMMARY',json.dumps(dict(parts=len(PARTS),mass_g=mass,cg=cg,ballast_g=ballast)),flush=True)
