"""Minimal Poisson bracket verification using only stdlib math."""
import math

MU = 1.0

def solve_kepler(M, e):
    E = M + e * math.sin(M)
    for _ in range(100):
        dE = (M - E + e * math.sin(E)) / (1 - e * math.cos(E))
        E += dE
        if abs(dE) < 1e-15: break
    return E

def true_anom(E, e):
    return 2 * math.atan2(math.sqrt(1+e)*math.sin(E/2), math.sqrt(1-e)*math.cos(E/2))

def S1(L,G,H,l,g,h):
    a=L**2/MU; eta=G/L; e=math.sqrt(max(1-eta**2,0)); theta=H/G
    E=solve_kepler(l,e); f=true_anom(E,e); r=a*(1-e*math.cos(E))
    Gam=MU**2/(a**3*eta**3); B0=-0.5+1.5*theta**2; s2=1-theta**2; phi=f-l
    return -Gam*(B0*(phi+e*math.sin(f))+0.75*s2*math.sin(2*f+2*g)+0.75*e*s2*math.sin(f+2*g)+0.25*e*s2*math.sin(3*f+2*g))

def p1(L,G,H,l,g,h):
    a=L**2/MU; e=math.sqrt(max(1-(G/L)**2,0)); E=solve_kepler(l,e); r=a*(1-e*math.cos(E))
    return L*(2*a/r-1)

def p2(L,G,H,l,g,h): return G
def p3(L,G,H,l,g,h): return H

def pb_fd(fp,fs,L,G,H,l,g,h,eps=1e-7):
    v=[L,G,H,l,g,h]; res=0.0
    for j in range(3):
        vp=list(v); vm=list(v); vp[3+j]+=eps; vm[3+j]-=eps
        dp_dlj=(fp(*vp)-fp(*vm))/(2*eps)
        vp=list(v); vm=list(v); vp[j]+=eps; vm[j]-=eps
        ds_dLj=(fs(*vp)-fs(*vm))/(2*eps)
        vp=list(v); vm=list(v); vp[j]+=eps; vm[j]-=eps
        dp_dLj=(fp(*vp)-fp(*vm))/(2*eps)
        vp=list(v); vm=list(v); vp[3+j]+=eps; vm[3+j]-=eps
        ds_dlj=(fs(*vp)-fs(*vm))/(2*eps)
        res+=dp_dlj*ds_dLj-dp_dLj*ds_dlj
    return res

def dS1_dg_an(L,G,H,l,g,h):
    a=L**2/MU; eta=G/L; e=math.sqrt(max(1-eta**2,0)); theta=H/G
    E=solve_kepler(l,e); f=true_anom(E,e)
    Gam=MU**2/(a**3*eta**3); B1p=1.5*(1-theta**2)
    return -Gam*B1p*(math.cos(2*f+2*g)+e*math.cos(f+2*g)+(e/3)*math.cos(3*f+2*g))

# Analytical {p1,S1} assembly
def pb_p1_S1_an(L,G,H,l,g,h):
    a=L**2/MU; eta=G/L; e=math.sqrt(max(1-eta**2,0)); theta=H/G
    E=solve_kepler(l,e); f=true_anom(E,e); r=a*(1-e*math.cos(E))
    Gam=MU**2/(a**3*eta**3); B0=-0.5+1.5*theta**2; B1p=1.5*(1-theta**2); s2=1-theta**2; phi=f-l
    # dp1/dl
    dr_dl=a*e*math.sin(E)/(1-e*math.cos(E))
    dp1dl=-2*L*a/r**2*dr_dl
    # dp1/dL
    de_dL=eta**2/(e*L) if e>1e-14 else 0
    dE_dL=de_dL*math.sin(E)/(1-e*math.cos(E))
    da_dL=2*a/L
    dr_dL=da_dL*(1-e*math.cos(E))+a*(-de_dL*math.cos(E)+e*math.sin(E)*dE_dL)
    dp1dL=(2*a/r-1)+L*(2*da_dL/r-2*a*dr_dL/r**2)
    # dp1/dG
    de_dG=-eta/(e*L) if e>1e-14 else 0
    dE_dG=de_dG*math.sin(E)/(1-e*math.cos(E))
    dr_dG=a*(-de_dG*math.cos(E)+e*math.sin(E)*dE_dG)
    dp1dG=-2*L*a/r**2*dr_dG
    # dS1/dl
    dS1dl=MU**2*B0/(a**3*eta**3)-MU**2/r**3*(B0+B1p*math.cos(2*f+2*g))
    # dS1/dg
    dS1dg=-Gam*B1p*(math.cos(2*f+2*g)+e*math.cos(f+2*g)+(e/3)*math.cos(3*f+2*g))
    # dS1/dL (full chain rule)
    dGam_dL=-3*Gam/L
    df_dE=a*eta/r
    if abs(math.sin(E))<1e-14: df_de=0.0
    else: df_de=math.sin(E)/(eta*(1-e*math.cos(E)))
    df_dL=df_dE*dE_dL+df_de*de_dL
    F=B0*(phi+e*math.sin(f))+0.75*s2*math.sin(2*f+2*g)+0.75*e*s2*math.sin(f+2*g)+0.25*e*s2*math.sin(3*f+2*g)
    dF_dL=(B0*(df_dL+de_dL*math.sin(f)+e*math.cos(f)*df_dL)
           +0.75*s2*2*math.cos(2*f+2*g)*df_dL
           +0.75*s2*(de_dL*math.sin(f+2*g)+e*math.cos(f+2*g)*df_dL)
           +0.25*s2*(de_dL*math.sin(3*f+2*g)+3*e*math.cos(3*f+2*g)*df_dL))
    dS1dL=-dGam_dL*F-Gam*dF_dL
    return dp1dl*dS1dL-dp1dL*dS1dl-dp1dG*dS1dg

print("="*90)
print("{p1, S1}: Analytical vs Finite Differences")
print("="*90)
print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'Analytical':>18} {'FD':>18} {'RelErr':>12} {'St':>4}")
print("-"*80)
mx=0; pc=0; fc=0
for ev in [0.01,0.1,0.3]:
    for Id in [30,60,85]:
        for gd in [0,45,90]:
            for lv in [0.5,1.5,3.0]:
                etv=math.sqrt(1-ev**2); Lv=math.sqrt(MU); Gv=Lv*etv; Hv=Gv*math.cos(math.radians(Id))
                gv=math.radians(gd)
                an=pb_p1_S1_an(Lv,Gv,Hv,lv,gv,0)
                fd=pb_fd(p1,S1,Lv,Gv,Hv,lv,gv,0)
                ref=max(abs(an),abs(fd),1e-20)
                er=abs(an-fd)/ref; mx=max(mx,er)
                st="PASS" if er<1e-5 else "FAIL"
                if st=="PASS": pc+=1
                else: fc+=1
                print(f"{ev:6.2f} {Id:5d} {gd:5d} {lv:5.1f} | {an:+18.10e} {fd:+18.10e} {er:12.2e} {st:>4}")
print(f"\n{{p1,S1}}: {pc} PASS, {fc} FAIL, max err={mx:.2e}")

print("\n"+"="*90)
print("{p2, S1} = -dS1/dg: Analytical vs FD")
print("="*90)
print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'Analytical':>18} {'FD':>18} {'RelErr':>12} {'St':>4}")
print("-"*80)
mx2=0; pc2=0; fc2=0
for ev in [0.01,0.1,0.3]:
    for Id in [30,60,85]:
        for gd in [0,45,90]:
            for lv in [0.5,1.5,3.0]:
                etv=math.sqrt(1-ev**2); Lv=math.sqrt(MU); Gv=Lv*etv; Hv=Gv*math.cos(math.radians(Id))
                gv=math.radians(gd)
                an=-dS1_dg_an(Lv,Gv,Hv,lv,gv,0)
                fd=pb_fd(p2,S1,Lv,Gv,Hv,lv,gv,0)
                ref=max(abs(an),abs(fd),1e-20)
                er=abs(an-fd)/ref; mx2=max(mx2,er)
                st="PASS" if er<1e-5 else "FAIL"
                if st=="PASS": pc2+=1
                else: fc2+=1
                print(f"{ev:6.2f} {Id:5d} {gd:5d} {lv:5.1f} | {an:+18.10e} {fd:+18.10e} {er:12.2e} {st:>4}")
print(f"\n{{p2,S1}}: {pc2} PASS, {fc2} FAIL, max err={mx2:.2e}")

print("\n"+"="*90)
print("{p3, S1} = 0: FD verification")
print("="*90)
mx3=0
for ev in [0.01,0.1,0.3]:
    for Id in [30,60,85]:
        for gd in [0,45,90]:
            for lv in [0.5,1.5,3.0]:
                etv=math.sqrt(1-ev**2); Lv=math.sqrt(MU); Gv=Lv*etv; Hv=Gv*math.cos(math.radians(Id))
                gv=math.radians(gd)
                fd=pb_fd(p3,S1,Lv,Gv,Hv,lv,gv,0)
                mx3=max(mx3,abs(fd))
print(f"Max |{{p3,S1}}| = {mx3:.2e} [{'PASS' if mx3<1e-10 else 'FAIL'}]")
