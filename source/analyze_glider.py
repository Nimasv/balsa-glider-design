"""Measured-section nonlinear lifting-line + tail trim and beam screening.
This is an engineering estimate, not CFD, a calibrated flight simulator or certification.
UIUC data produced under the UIUC Low-Speed Airfoil Test program; included license applies.
"""
from pathlib import Path
import sys,json,csv,math
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'.cad_deps'))
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = ROOT / 'results'
OUT.mkdir(parents=True, exist_ok=True)
def readpol(filename):
    lines=(ROOT/'source/aero_data'/filename).read_text().splitlines();out=[];i=0
    while i<len(lines):
        if lines[i].strip().startswith('Average Reynolds'):
            re=float(lines[i+1]);n=int(lines[i+3]);rows=[]
            for j in range(n):rows.append([float(x) for x in lines[i+5+j].split()][:3])
            a=np.array(rows);d=np.where(np.diff(a[:,0])<0)[0]
            if len(d):a=a[:d[0]+1]
            out.append((re,a));i+=5+n
        else:i+=1
    return sorted(out,key=lambda q:q[0])
LIFT=readpol('SD7037D.LFT');DRAG=readpol('SD7037D.DRG')
def interp(data,re,alpha,col):
    vals=[np.interp(alpha,a[:,0],a[:,col]) for _,a in data]
    return float(np.interp(re,[r for r,_ in data],vals))
def arrayinterp(data,re,alpha,col):return np.array([interp(data,r,a,col) for r,a in zip(np.atleast_1d(re),np.atleast_1d(alpha))])
massdata=json.loads((OUT/'mass_balance.json').read_text());geo=json.loads((OUT/'geometry.json').read_text())
M=massdata['mass_g']/1000;rho=1.225;mu=1.789e-5;g=9.80665
b=1.7;S=.255;AR=b*b/S;MAC=geo['MAC']/1000;Sh=.039;Bh=.4;ARh=Bh*Bh/Sh
xcg=massdata['cg_mm'][0]/1000;hac=(xcg-.220)/(MAC*np.cos(np.deg2rad(2)));htail=geo['tail_ac_x']/1000
ltcg=htail-xcg;eta=.9;it=0.;iw=np.deg2rad(2);tau=.55
ah=2*np.pi*ARh/(2+np.sqrt(4+ARh**2));N=14;n=np.arange(1,2*N,2)
theta=np.linspace(np.pi/(2*N),np.pi/2,N);y=.5*b*np.cos(theta)
c=.17-.04*y/(b/2);B=np.sin(theta[:,None]*n[None,:]);D=B*n[None,:]/np.sin(theta[:,None])
yd=np.linspace(0,b/2,301);td=np.arccos(2*yd/b);cd=.17-.04*yd/(b/2)
Bd=np.sin(td[:,None]*n[None,:]);sint=np.maximum(np.sin(td),1e-6)
Dd=Bd*n[None,:]/sint[:,None]
def wing(A,alpha,V):
    CL=np.pi*AR*A[0];cl=4*b*(Bd@A)/cd
    eff=np.rad2deg(alpha+iw-Dd@A);Re=rho*V*cd/mu
    cm=arrayinterp(LIFT,Re,eff,2)
    Cmw=2*np.trapezoid(cm*cd**2,yd)/(S*MAC)
    Cdi=np.pi*AR*sum(n*A*A)
    # Spanwise section drag, no extrapolation outside tested alpha/Re domain.
    cdp=arrayinterp(DRAG,Re,eff,2)
    Cdp=2*np.trapezoid(cdp*cd,yd)/S
    return CL,Cmw,Cdi,Cdp,cl,eff,Re
def solve(V,initial=None):
    Req=M*g/(.5*rho*V*V*S);Re=rho*V*c/mu
    def residual(x):
        A=x[:N];alpha=x[N];CLh=x[N+1]
        cl=4*b*(B@A)/c;eff=np.rad2deg(alpha+iw-D@A)
        cldata=arrayinterp(LIFT,Re,eff,1)
        CL,Cm,*_=wing(A,alpha,V)
        return np.r_[cl-cldata,CL+eta*Sh/S*CLh-Req,Cm+CL*(hac-.25)-eta*Sh/S*CLh*ltcg/MAC]
    cl0=arrayinterp(LIFT,Re,np.zeros(N),1)
    slope=(arrayinterp(LIFT,Re,np.ones(N)*5,1)-cl0)/np.deg2rad(5)
    mat=np.zeros((N+2,N+2));rhs=np.zeros(N+2)
    mat[:N,:N]=4*b*B/c[:,None]+slope[:,None]*D
    mat[:N,N]=-slope;rhs[:N]=cl0+slope*iw
    mat[N,0]=np.pi*AR;mat[N,N+1]=eta*Sh/S;rhs[N]=Req
    mat[N+1,0]=(hac-.25)*np.pi*AR;mat[N+1,N+1]=-eta*Sh/S*ltcg/MAC;rhs[N+1]=.08
    x0=np.linalg.solve(mat,rhs)
    sol=least_squares(residual,x0,method='lm',max_nfev=400,xtol=2e-10,ftol=2e-10,gtol=2e-10)
    residual_error=float(max(abs(residual(sol.x))))
    A=sol.x[:N];alpha=sol.x[N];CLh=sol.x[N+1]
    CL,Cm,Cdi,Cdp,cl,eff,Re=wing(A,alpha,V)
    epsilon=2*CL/(np.pi*AR)
    delta=(CLh/ah-(alpha+it-epsilon))/tau
    tailCD=eta*Sh/S*(.012+CLh**2/(np.pi*.8*ARh))
    finCD=.01995/S*.012
    # Body Cf=0.004, FF=1.3, wet area estimated from CAD panels.
    bodyCD=.004*1.3*.15/S
    extrasCD=.003 # seams, rough covering and external horns/guide tubes
    CD=Cdi+Cdp+tailCD+finCD+bodyCD+extrasCD
    CLtot=CL+eta*Sh/S*CLh;LD=CLtot/CD;gamma=math.atan(1/LD)
    q=.5*rho*V*V
    row=dict(V_m_s=V,CL_total=CLtot,CL_wing=CL,Cm_wing=Cm,CL_tail=CLh,CD_total=CD,
             CD_profile_wing=Cdp,CD_induced_wing=Cdi,alpha_body_deg=float(np.rad2deg(alpha)),
             elevator_trim_deg=float(np.rad2deg(delta)),L_D=LD,sink_m_s=V*math.sin(gamma),
             Re_min=float(Re.min()),Re_max=float(Re.max()),max_section_cl=float(max(cl)),
             section_data_Re_clamped=bool(min(Re)<min(r for r,_ in LIFT)),
             residual=residual_error,lift_N=q*S*CLtot,tail_lift_N=q*eta*Sh*CLh)
    # Fixed-elevator local static derivative with measured section response.
    def solvewing(al):
        def res(A):return 4*b*(B@A)/c-arrayinterp(LIFT,Re[:N] if False else rho*V*c/mu,np.rad2deg(al+iw-D@A),1)
        return least_squares(res,A,max_nfev=60,gtol=1e-10).x
    da=np.deg2rad(.15);ap=solvewing(alpha+da);am=solvewing(alpha-da)
    wp=wing(ap,alpha+da,V);wm=wing(am,alpha-da,V)
    clprime=(wp[0]-wm[0])/(2*da);cmprime=(wp[1]-wm[1])/(2*da)
    de=2*clprime/(np.pi*AR);hprime=ah*(1-de)
    totalprime=clprime+eta*Sh/S*hprime
    cmalpha=cmprime+clprime*(hac-.25)-eta*Sh/S*hprime*ltcg/MAC+.08
    row.update(CL_alpha_per_rad=totalprime,Cm_alpha_per_rad=cmalpha,static_margin=-cmalpha/totalprime,
               neutral_point_x_mm=(xcg-cmalpha/totalprime*MAC)*1000,downwash_gradient=de)
    return row,sol.x,(yd,cl,cd,A)

rows=[];last=None;solutions={}
for V in np.arange(6.5,10.51,.25):
    row,last,dist=solve(float(V),last);rows.append(row);solutions[round(float(V),2)]=dist
    print('AERO',V,round(row['L_D'],2),round(row['elevator_trim_deg'],2),round(row['static_margin'],3),'res',row['residual'],'Reclamp',row['section_data_Re_clamped'],flush=True)
valid=[r for r in rows if r['residual']<.002 and abs(r['elevator_trim_deg'])<=15 and not r['section_data_Re_clamped']]
if not valid:raise RuntimeError('No valid trim point inside modeled elevator travel/data Re range')
best=max(valid,key=lambda r:r['L_D']);cruise=min(valid,key=lambda r:abs(r['V_m_s']-7.5))
with (OUT/'performance.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

# Structural screening of two caps, D-box skins and shear webs, elastic Euler beam.
foil=np.loadtxt(ROOT/'source/sd7037.dat',skiprows=1);split=np.argmin(foil[:,0])
def u(x):return np.interp(x,foil[:split+1,0][::-1],foil[:split+1,1][::-1])
def l(x):return np.interp(x,foil[split:,0],foil[split:,1])
ys=np.linspace(0,850,201);ch=.17-.04*ys/850
# Elliptic distribution provides a transparent conservative design case near this taper.
ell=np.sqrt(np.maximum(0,1-(ys/850)**2))
wing_load_factor=1+abs(min(0,cruise['tail_lift_N']))/(M*g)
load=ell/(2*np.trapezoid(ell,ys))*M*g*wing_load_factor
Vshear=np.array([np.trapezoid(load[i:],ys[i:]) for i in range(len(ys))])
mom=np.array([np.trapezoid(load[i:]*(ys[i:]-ys[i]),ys[i:]) for i in range(len(ys))])
EI=[];extreme=[];torsionGJ=[];boxarea=[]
for yy,ccm in zip(ys,ch):
    cc=ccm*1000;sp=.3*cc
    # Integrate E-weighted rectangular fibres over cap/skin cross section.
    elements=[]
    for mode in ['top','bottom']:
        for x in np.linspace(sp-10,sp+10,41):
            surf=cc*(u(x/cc) if mode=='top' else l(x/cc))
            z=surf+(-2.5 if mode=='top' else 2.5)
            elements.append((3400,20/41*3,z,20/41*3**3/12,'cap'))
        for x in np.linspace(4,sp+10,81):
            surf=cc*(u(x/cc) if mode=='top' else l(x/cc));z=surf+(-.5 if mode=='top' else .5)
            dx=(sp+6)/81;elements.append((2200,dx,z,dx/12,'skin'))
    for x in [sp-9.25,sp+9.25]:
        top=cc*u(x/cc)-4;bottom=cc*l(x/cc)+4;h=max(.2,top-bottom)
        elements.append((300,1.5*h,(top+bottom)/2,1.5*h**3/12,'web'))
    zn=sum(E*a*z for E,a,z,I,kind in elements)/sum(E*a for E,a,z,I,kind in elements)
    stiff=sum(E*(I+a*(z-zn)**2) for E,a,z,I,kind in elements);EI.append(stiff)
    extreme.append(max(3400*(abs(z-zn)+1.5)/stiff for E,a,z,I,kind in elements if kind=='cap'))
    xx=np.linspace(2,sp+10,81);top=cc*np.array([u(x/cc) for x in xx])-.5;bot=cc*np.array([l(x/cc) for x in xx])+.5
    area=np.trapezoid(top-bot,xx);perim=np.sum(np.sqrt(np.diff(xx)**2+np.diff(top)**2))+np.sum(np.sqrt(np.diff(xx)**2+np.diff(bot)**2))+(top[-1]-bot[-1])
    torsionGJ.append(4*area**2/(perim/150));boxarea.append(area)
EI=np.array(EI);extreme=np.array(extreme);GJ=np.array(torsionGJ)
struct={}
for nload in [1,4,6]:
    curvature=mom*nload/EI;slope=cumulative_trapezoid(curvature,ys,initial=0);defl=cumulative_trapezoid(slope,ys,initial=0)
    stress=mom*nload*extreme
    struct[str(nload)]=dict(root_moment_Nm=mom[0]*nload/1000,root_shear_N=Vshear[0]*nload,
                          max_cap_stress_MPa=max(stress),tip_deflection_mm=defl[-1])
joiner_stress=struct['4']['root_moment_Nm']*1000/(2*4*11**2/6)
# Aerodynamic twisting moment at 10.5 m/s, Cm=-0.10 and 0.05 chord lift arm.
q=.5*rho*10.5**2;torque_per_mm=q*((ch**2)*(.10+.05*1.0)) # N.mm/mm
T=np.array([np.trapezoid(torque_per_mm[i:],ys[i:]) for i in range(len(ys))])
twist=np.trapezoid(T/GJ,ys)
summary=dict(model='Nonlinear lifting line with UIUC SD7037(D) section data + analytical tail; elastic beam screening',
             mass=massdata,geometry=geo,cruise=cruise,best_glide_in_scanned_range=best,
             data_reynolds=[r for r,_ in LIFT],structural_cases=struct,
             wing_load_factor_including_tail=wing_load_factor,joiner_stress_4g_MPa=joiner_stress,joiner_stress_6g_MPa=1.5*joiner_stress,estimated_twist_at_10_5m_s_deg=float(np.rad2deg(twist)),
             design_allowables_assumed=dict(cap_service_MPa=7.5,cap_ultimate_MPa=12,plywood_service_MPa=20,plywood_ultimate_MPa=30,cap_E_MPa=3400,skin_E_MPa=2200),
             criteria=dict(mass_in_range=320<=massdata['mass_g']<=480,wing_loading_in_range=14<=massdata['wing_loading_g_dm2']<=20,
                           cg_in_range=.30-1e-4<=hac<=.38,positive_static_margin=cruise['static_margin']>0,
                           cap_service=struct['4']['max_cap_stress_MPa']<=7.5,cap_ultimate=struct['6']['max_cap_stress_MPa']<=12,
                           joiner_service=joiner_stress<=20,joiner_ultimate=joiner_stress*1.5<=30),
             limitations=['No flight calibration; no CFD or full dynamic mode analysis.',
                          'Balsa/ply allowables and equipment masses require sample verification.',
                          'Tail effectiveness, downwash, body moment/drag and roughness remain assumptions.',
                          'L=W approximation neglects small cos(gamma) correction.',
                          'Section data clamped outside Reynolds bounds; such points excluded from selected best glide.',
                          'No guarantee of travel distance without launch height/speed and course definition.'])
(OUT/'engineering_results.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item()),encoding='utf-8')
fig,axs=plt.subplots(2,2,figsize=(11,7));vs=[r['V_m_s'] for r in rows]
for ax,key,label in zip(axs.flat,['L_D','sink_m_s','elevator_trim_deg','static_margin'],['Lift / drag','Sink rate [m/s]','Elevator trim [deg, down +]','Static margin [MAC fraction]']):
    ax.plot(vs,[r[key] for r in rows],'-o',ms=3);ax.set_xlabel('Airspeed [m/s]');ax.set_ylabel(label);ax.grid(alpha=.3)
fig.suptitle('R03 analytical predictions — measured section data, assumed tail/body terms');fig.tight_layout();fig.savefig(OUT/'performance.png',dpi=160);plt.close(fig)
print(json.dumps(summary['criteria'],indent=2,default=lambda x:x.item()));print('CRUISE',cruise);print('STRUCTURE',struct)
