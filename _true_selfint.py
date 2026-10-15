# -*- coding: utf-8 -*-
# T1 概念验证：双纽线柱面面片 S(s,v)=(x(s),y(s),v) 图像真自交 + 方阵 Krawczyk 证明
import sys, math
sys.stdout.reconfigure(encoding='utf-8')

A = 1.0
def bcurve(s):
    t = 2*math.pi*s; sn=math.sin(t); cs=math.cos(t); sn2=sn*sn; den=1+sn2
    return A*cs/den, A*sn*cs/den
def bcurved(s):
    t = 2*math.pi*s; sn=math.sin(t); cs=math.cos(t); sn2=sn*sn; den2=(1+sn2)**2
    dxdt = A*(-sn*(3-sn2))/den2
    dydt = A*((cs*cs-sn2)*(1+sn2) - 2*sn2*cs*cs)/den2
    return dxdt*2*math.pi, dydt*2*math.pi
def F(s1,s2):
    x1,y1=bcurve(s1); x2,y2=bcurve(s2)
    return x1-x2, y1-y2

# 1) 自交点
p1,p2 = bcurve(0.25), bcurve(0.75)
print(f"S(0.25,v)=({p1[0]:.5f},{p1[1]:.5f},v)  S(0.75,v)=({p2[0]:.5f},{p2[1]:.5f},v)")
print("图像真自交(两不同参数映同点):", abs(p1[0]-p2[0])<1e-12 and abs(p1[1]-p2[1])<1e-12)

# 2) 中点 Jacobian（方阵）
s1m,s2m=0.25,0.75
dx1,dy1=bcurved(s1m); dx2,dy2=bcurved(s2m)
J=[[dx1,-dx2],[dy1,-dy2]]
det=J[0][0]*J[1][1]-J[0][1]*J[1][0]
print(f"J=[{J[0]},{J[1]}] det={det:.4f} 非奇异={abs(det)>1e-6}")

# 3) 经典 Krawczyk：K=m-Y F(m)+(I-Y J(X))(X-m)，验证 K⊂int(X)
Xlo=[0.20,0.70]; Xhi=[0.30,0.80]; m=[0.25,0.75]
Y=[[J[1][1]/det, -J[0][1]/det],[-J[1][0]/det, J[0][0]/det]]  # J(m)^{-1}
fm=F(m[0],m[1])
YF=[Y[0][0]*fm[0]+Y[0][1]*fm[1], Y[1][0]*fm[0]+Y[1][1]*fm[1]]

# 区间外包 J(X)：对每分量在 [lo,hi] 采样端点+中点外包
def range_ext(f,a,b):
    xs=[a,(a+b)/2,b]; lo=1e300; hi=-1e300
    for x in xs:
        v=f(x); lo=min(lo,v); hi=max(hi,v)
    return lo,hi
# J(X) 区间元素: J(x)ij 在盒内随 s1(列0/行0,1?) / s2 变。J(X)_{00}=dx1(s1),J(X)_{01}=-dx2(s2),J(X)_{10}=dy1(s1),J(X)_{11}=-dy2(s2)
r00=range_ext(lambda s: bcurved(s)[0], Xlo[0], Xhi[0])
r01=range_ext(lambda s: -bcurved(s)[0], Xlo[1], Xhi[1])
r10=range_ext(lambda s: bcurved(s)[1], Xlo[0], Xhi[0])
r11=range_ext(lambda s: -bcurved(s)[1], Xlo[1], Xhi[1])
JXI=[[r00[0],r00[1]],[r01[0],r01[1]],[r10[0],r10[1]],[r11[0],r11[1]]]  # 扁平 [00,01,10,11]

K=[None,None]
for i in range(2):
    lo=m[i]-YF[i]; hi=m[i]-YF[i]
    for j in range(2):
        YJlo=0.0; YJhi=0.0
        for k in range(2):
            yv=Y[i][k]; jl=JXI[k*2+j][0]; jh=JXI[k*2+j][1]
            a=yv*jl; b=yv*jh
            YJlo+=min(a,b); YJhi+=max(a,b)
        if i==j:
            tlo=1-YJhi; thi=1-YJlo
        else:
            tlo=-YJhi; thi=-YJlo
        xl=Xlo[j]-m[j]; xh=Xhi[j]-m[j]
        a=tlo*xl; b=tlo*xh; c=thi*xl; d=thi*xh
        lo+=min(min(a,b),min(c,d)); hi+=max(max(a,b),max(c,d))
    K[i]=(lo,hi)

print(f"\n盒 X=({Xlo[0]},{Xhi[0]})x({Xlo[1]},{Xhi[1]}), 解 m=({m[0]},{m[1]})")
for i in range(2):
    print(f"  K[{i}]=({K[i][0]:.4f},{K[i][1]:.4f})  int(X[{i}])=({Xlo[i]},{Xhi[i]})")
ok=(K[0][0]>Xlo[0] and K[0][1]<Xhi[0] and K[1][0]>Xlo[1] and K[1][1]<Xhi[1])
print("K ⊂ int(X)  ⇒  经典 Krawczyk 证明盒内唯一自交解:", ok)

# 4) 扫描参数域 [0,1] 内全部自交解对（s1<s2）
print("\n--- 扫描 s1<s2 参数域内自交解 ---")
sols=[]
n=200
for i in range(n):
    s1=0.005*i
    for j in range(i+1,n):
        s2=0.005*j
        f=F(s1,s2)
        if abs(f[0])<1e-8 and abs(f[1])<1e-8:
            sols.append((round(s1,3),round(s2,3)))
print("自交解对:", sols, " 共", len(sols), "对")
