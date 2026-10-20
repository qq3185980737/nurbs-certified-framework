// bench_step_case2.cpp — 真实 STEP 零件数据上的认证曲线-曲面求交（牛顿定位 + 3x3 Krawczyk 认证）
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>

struct IV { double lo, hi; IV(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} double wid()const{return hi-lo;} };
static IV ival(double l,double h){ return IV(l,h); }
static IV iadd(const IV&a,const IV&b){ return IV(std::nextafter(a.lo+b.lo,-1e300),std::nextafter(a.hi+b.hi,1e300)); }
static IV isub(const IV&a,const IV&b){ return IV(std::nextafter(a.lo-b.hi,-1e300),std::nextafter(a.hi-b.lo,1e300)); }
static IV imul(const IV&a,const IV&b){ double p[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi}; double mn=p[0],mx=p[0]; for(int i=1;i<4;++i){mn=std::min(mn,p[i]);mx=std::max(mx,p[i]);} return IV(std::nextafter(mn,-1e300),std::nextafter(mx,1e300)); }
static bool idiv(const IV&a,const IV&b,IV&out){ if(b.lo<=0&&b.hi>=0) return false; double p[4]={a.lo/b.lo,a.lo/b.hi,a.hi/b.lo,a.hi/b.hi}; double mn=p[0],mx=p[0]; for(int i=1;i<4;++i){mn=std::min(mn,p[i]);mx=std::max(mx,p[i]);} out=IV(std::nextafter(mn,-1e300),std::nextafter(mx,1e300)); return true; }
static bool has0(const IV&a){ return a.lo<=0.0&&a.hi>=0.0; }

using H4 = std::array<double,4>;

struct Curve {
    int deg=0; std::vector<H4> P;
    H4 pt(double t) const {
        if(P.size()<2) return P.empty()?H4{0,0,0,1}:P[0];
        std::vector<H4> c=P;
        for(size_t r=1;r<c.size();++r) for(size_t i=0;i+r<c.size();++i) for(int k=0;k<4;++k)
            c[i][k]=(1-t)*c[i][k]+t*c[i+1][k];
        return c[0];
    }
    std::array<IV,4> iv(const IV& t) const {
        std::array<IV,4> z; for(int k=0;k<4;++k) z[k]=ival(0,0);
        if(P.size()<2) return z;
        std::vector<std::array<IV,4>> c(P.size());
        for(size_t i=0;i<P.size();++i) for(int k=0;k<4;++k) c[i][k]=ival(P[i][k],P[i][k]);
        IV one(1.0,1.0);
        for(size_t r=1;r<c.size();++r) for(size_t i=0;i+r<c.size();++i) for(int k=0;k<4;++k)
            c[i][k]=iadd(imul(c[i][k],isub(one,t)),imul(c[i+1][k],t));
        return c[0];
    }
    std::array<IV,4> div(const IV& t) const {
        std::array<IV,4> z; for(int k=0;k<4;++k) z[k]=ival(0,0);
        int p=deg; if(p<1||(int)P.size()<2) return z;
        std::vector<std::array<IV,4>> c(p);
        for(int i=0;i<p;++i) for(int k=0;k<4;++k)
            c[i][k]=imul(ival((double)p,(double)p), isub(ival(P[i+1][k],P[i+1][k]), ival(P[i][k],P[i][k])));
        IV one(1.0,1.0);
        for(size_t r=1;r<c.size();++r) for(size_t i=0;i+r<c.size();++i) for(int k=0;k<4;++k)
            c[i][k]=iadd(imul(c[i][k],isub(one,t)),imul(c[i+1][k],t));
        return c[0];
    }
};

struct Surf {
    int du=0,dv=0,nu=0,nv=0; std::vector<H4> P;
    H4 pt(double s,double t) const {
        if(nu<1||nv<1) return H4{0,0,0,1};
        std::vector<H4> col(nu);
        for(int i=0;i<nu;++i){
            std::vector<H4> row(nv);
            for(int j=0;j<nv;++j) row[j]=P[(size_t)i*nv+j];
            for(int r=1;r<nv;++r) for(int j=0;j+r<nv;++j) for(int k=0;k<4;++k)
                row[j][k]=(1-t)*row[j][k]+t*row[j+1][k];
            col[i]=row[0];
        }
        for(int r=1;r<nu;++r) for(int i=0;i+r<nu;++i) for(int k=0;k<4;++k)
            col[i][k]=(1-s)*col[i][k]+s*col[i+1][k];
        return col[0];
    }
    std::array<IV,4> iv(const IV& s,const IV& t) const {
        std::array<IV,4> z; for(int k=0;k<4;++k) z[k]=ival(0,0);
        if(nu<1||nv<1) return z;
        IV one(1.0,1.0);
        std::vector<std::array<IV,4>> col(nu);
        for(int i=0;i<nu;++i){
            std::vector<std::array<IV,4>> row(nv);
            for(int j=0;j<nv;++j) for(int k=0;k<4;++k) row[j][k]=ival(P[(size_t)i*nv+j][k],P[(size_t)i*nv+j][k]);
            for(int r=1;r<nv;++r) for(int j=0;j+r<nv;++j) for(int k=0;k<4;++k)
                row[j][k]=iadd(imul(row[j][k],isub(one,t)),imul(row[j+1][k],t));
            col[i]=row[0];
        }
        for(int r=1;r<nu;++r) for(int i=0;i+r<nu;++i) for(int k=0;k<4;++k)
            col[i][k]=iadd(imul(col[i][k],isub(one,s)),imul(col[i+1][k],s));
        return col[0];
    }
    std::array<IV,4> ds(const IV& s,const IV& t) const {
        std::array<IV,4> z; for(int k=0;k<4;++k) z[k]=ival(0,0);
        if(nu<1||nv<1||du<1) return z;
        Surf d; d.du=du-1; d.dv=dv; d.nu=nu-1; d.nv=nv; d.P.resize((size_t)(nu-1)*nv);
        for(int i=0;i<nu-1;++i) for(int j=0;j<nv;++j) for(int k=0;k<4;++k)
            d.P[(size_t)i*nv+j][k]=(double)du*(P[(size_t)(i+1)*nv+j][k]-P[(size_t)i*nv+j][k]);
        return d.iv(s,t);
    }
    std::array<IV,4> dt(const IV& s,const IV& t) const {
        std::array<IV,4> z; for(int k=0;k<4;++k) z[k]=ival(0,0);
        if(nu<1||nv<1||dv<1) return z;
        Surf d; d.du=du; d.dv=dv-1; d.nu=nu; d.nv=nv-1; d.P.resize((size_t)nu*(nv-1));
        for(int i=0;i<nu;++i) for(int j=0;j<nv-1;++j) for(int k=0;k<4;++k)
            d.P[(size_t)i*(nv-1)+j][k]=(double)dv*(P[(size_t)i*nv+j+1][k]-P[(size_t)i*nv+j][k]);
        return d.iv(s,t);
    }
};

static bool rat3(const std::array<IV,4>& h, IV out[3]){
    IV w=h[3]; if(has0(w)) return false;
    for(int i=0;i<3;++i) if(!idiv(h[i],w,out[i])) return false;
    return true;
}
static bool ratd(const std::array<IV,4>& h,const std::array<IV,4>& dh,const IV S[3],IV out[3]){
    for(int i=0;i<3;++i){ IV num=isub(dh[i],imul(S[i],dh[3])); if(!idiv(num,h[3],out[i])) return false; }
    return true;
}

static int g_cert=0, g_refused=0, g_newton_fail=0;

static bool certify(const Curve& C,const Surf& S,double u0,double s0,double t0,double outPt[3]){
    double x[3]={u0,s0,t0};
    for(int it=0;it<60;++it){
        H4 a=C.pt(x[0]), b=S.pt(x[1],x[2]);
        if(std::fabs(a[3])<1e-14||std::fabs(b[3])<1e-14) return false;
        double F[3], nrm=0;
        for(int k=0;k<3;++k){ F[k]=a[k]/a[3]-b[k]/b[3]; nrm+=F[k]*F[k]; }
        nrm=std::sqrt(nrm);
        if(nrm<1e-14) break;
        double J[3][3];
        const double h=1e-8;
        for(int c=0;c<3;++c){
            double xp[3]={x[0],x[1],x[2]}, xm[3]={x[0],x[1],x[2]};
            xp[c]=std::min(1.0,x[c]+h); xm[c]=std::max(0.0,x[c]-h);
            H4 ap=C.pt(xp[0]), bp=S.pt(xp[1],xp[2]);
            H4 am=C.pt(xm[0]), bm=S.pt(xm[1],xm[2]);
            if(std::fabs(ap[3])<1e-14||std::fabs(bp[3])<1e-14||std::fabs(am[3])<1e-14||std::fabs(bm[3])<1e-14) return false;
            double Fp[3],Fm[3];
            for(int k=0;k<3;++k){ Fp[k]=ap[k]/ap[3]-bp[k]/bp[3]; Fm[k]=am[k]/am[3]-bm[k]/bm[3]; }
            double den=xp[c]-xm[c];
            for(int r=0;r<3;++r) J[r][c]=(Fp[r]-Fm[r])/den;
        }
        double A[3][6];
        for(int r=0;r<3;++r){ for(int c=0;c<3;++c){A[r][c]=J[r][c];A[r][3+c]=(r==c)?1.0:0.0;} }
        bool sing=false;
        for(int col=0;col<3;++col){
            int piv=col; double bst=std::fabs(A[col][col]);
            for(int r=col+1;r<3;++r) if(std::fabs(A[r][col])>bst){bst=std::fabs(A[r][col]);piv=r;}
            if(bst<1e-13){sing=true;break;}
            if(piv!=col) for(int j=0;j<6;++j) std::swap(A[piv][j],A[col][j]);
            double d=A[col][col]; for(int j=0;j<6;++j) A[col][j]/=d;
            for(int r=0;r<3;++r) if(r!=col){ double f=A[r][col]; for(int j=0;j<6;++j) A[r][j]-=f*A[col][j]; }
        }
        if(sing) return false;
        for(int c=0;c<3;++c){ double dx=0; for(int r=0;r<3;++r) dx+=A[c][3+r]*(-F[r]);
            x[c]=std::min(1.0,std::max(0.0,x[c]+dx)); }
    }
    // 残差诊断
    {
        H4 a=C.pt(x[0]), b=S.pt(x[1],x[2]);
        double r=0; for(int k=0;k<3;++k){ double d=a[k]/a[3]-b[k]/b[3]; r+=d*d; }
        static int dbg=0;
        if(++dbg<=8){ std::printf("    newton res=%.2e x=(%.4f,%.4f,%.4f)\n",std::sqrt(r),x[0],x[1],x[2]); std::fflush(stdout); }
    }
    // 认证盒
    const double w=1e-5;
    double lo[3]={std::max(0.0,x[0]-w),std::max(0.0,x[1]-w),std::max(0.0,x[2]-w)};
    double hi[3]={std::min(1.0,x[0]+w),std::min(1.0,x[1]+w),std::min(1.0,x[2]+w)};
    if(hi[0]-lo[0]<1e-5||hi[1]-lo[1]<1e-5||hi[2]-lo[2]<1e-5) return false;
    IV u(lo[0],hi[0]), s(lo[1],hi[1]), t(lo[2],hi[2]);
    auto hc=C.iv(u); auto hs=S.iv(s,t);
    IV Cc[3],Sc[3];
    if(!rat3(hc,Cc)||!rat3(hs,Sc)) return false;
    for(int k=0;k<3;++k){ IV F=isub(Cc[k],Sc[k]); if(!has0(F)) return false; }
    IV dC[3],dSs[3],dSt[3];
    if(!ratd(hc,C.div(u),Cc,dC)) return false;
    if(!ratd(hs,S.ds(s,t),Sc,dSs)) return false;
    if(!ratd(hs,S.dt(s,t),Sc,dSt)) return false;
    // 中点
    double mu=u.mid(), ms=s.mid(), mt=t.mid();
    H4 a=C.pt(mu), b=S.pt(ms,mt);
    if(std::fabs(a[3])<1e-14||std::fabs(b[3])<1e-14) return false;
    double Fm[3]={a[0]/a[3]-b[0]/b[3],a[1]/a[3]-b[1]/b[3],a[2]/a[3]-b[2]/b[3]};
    // 浮点 J，数值差分
    double Jm[3][3];
    const double h=1e-8;
    for(int c=0;c<3;++c){
        double xp[3]={mu,ms,mt}, xm[3]={mu,ms,mt};
        xp[c]=std::min(1.0,xp[c]+h); xm[c]=std::max(0.0,xm[c]-h);
        H4 ap=C.pt(xp[0]), bp=S.pt(xp[1],xp[2]);
        H4 am=C.pt(xm[0]), bm=S.pt(xm[1],xm[2]);
        double Fp[3],Fmm[3];
        for(int k=0;k<3;++k){ Fp[k]=ap[k]/ap[3]-bp[k]/bp[3]; Fmm[k]=am[k]/am[3]-bm[k]/bm[3]; }
        double den=xp[c]-xm[c];
        for(int r=0;r<3;++r) Jm[r][c]=(Fp[r]-Fmm[r])/den;
    }
    double M[3][6];
    for(int r=0;r<3;++r){ for(int c=0;c<3;++c){M[r][c]=Jm[r][c];M[r][3+c]=(r==c)?1.0:0.0;} }
    for(int col=0;col<3;++col){
        int piv=col; double bst=std::fabs(M[col][col]);
        for(int r=col+1;r<3;++r) if(std::fabs(M[r][col])>bst){bst=std::fabs(M[r][col]);piv=r;}
        if(bst<1e-13) return false;
        if(piv!=col) for(int j=0;j<6;++j) std::swap(M[piv][j],M[col][j]);
        double d=M[col][col]; for(int j=0;j<6;++j) M[col][j]/=d;
        for(int r=0;r<3;++r) if(r!=col){ double f=M[r][col]; for(int j=0;j<6;++j) M[r][j]-=f*M[col][j]; }
    }
    double Y[3][3];
    for(int r=0;r<3;++r) for(int c=0;c<3;++c) Y[r][c]=M[r][3+c];
    IV JX[3][3];
    JX[0][0]=dC[0]; JX[0][1]=IV(-dSs[0].hi,-dSs[0].lo); JX[0][2]=IV(-dSt[0].hi,-dSt[0].lo);
    JX[1][0]=dC[1]; JX[1][1]=IV(-dSs[1].hi,-dSs[1].lo); JX[1][2]=IV(-dSt[1].hi,-dSt[1].lo);
    JX[2][0]=dC[2]; JX[2][1]=IV(-dSs[2].hi,-dSs[2].lo); JX[2][2]=IV(-dSt[2].hi,-dSt[2].lo);
    IV v[3]={isub(u,IV(mu,mu)),isub(s,IV(ms,ms)),isub(t,IV(mt,mt))};
    IV Xc[3]={u,s,t};
    double m[3]={mu,ms,mt};
    for(int i=0;i<3;++i){
        double yf=0; for(int j=0;j<3;++j) yf+=Y[i][j]*Fm[j];
        IV K=isub(IV(m[i],m[i]),IV(yf,yf));
        for(int j=0;j<3;++j){
            IV yj(0,0);
            for(int l=0;l<3;++l) yj=iadd(yj,imul(IV(Y[i][l],Y[i][l]),JX[l][j]));
            IV mij=isub(IV((i==j)?1.0:0.0,(i==j)?1.0:0.0),yj);
            K=iadd(K,imul(mij,v[j]));
        }
        if(!(Xc[i].lo<K.lo && K.hi<Xc[i].hi)) { static int dbg2=0; if(++dbg2<=8){ std::printf("    K fail i=%d X=[%.2e,%.2e] K=[%.2e,%.2e]\n",i,Xc[i].lo,Xc[i].hi,K.lo,K.hi); std::fflush(stdout);} return false; }
    }
    {
        H4 a=C.pt(x[0]);
        outPt[0]=a[0]/a[3]; outPt[1]=a[1]/a[3]; outPt[2]=a[2]/a[3];
    }
    ++g_cert;
    return true;
}

int main(int argc,char** argv){
    const char* path = argc>1?argv[1]:"case_data.txt";
    std::ifstream f(path);
    if(!f){ std::printf("cannot open %s\n",path); return 1; }
    std::vector<std::vector<Curve>> curves;
    std::vector<std::vector<Surf>> surfs;
    struct Pair{int e,s; std::vector<std::array<double,3>> pts;};
    std::vector<Pair> pairs;
    std::string tok;
    while(f>>tok){
        if(tok=="NSURF"){int n;f>>n;surfs.resize(n);}
        else if(tok=="NCURVE"){int n;f>>n;curves.resize(n);}
        else if(tok=="SURF"){
            int fid,pi,du,dv,nu,nv; f>>fid>>pi>>du>>dv>>nu>>nv;
            Surf s; s.du=du; s.dv=dv; s.nu=nu; s.nv=nv; s.P.resize((size_t)nu*nv);
            for(int i=0;i<nu;++i) for(int j=0;j<nv;++j){ double x,y,z,w; f>>x>>y>>z>>w; s.P[(size_t)i*nv+j]=H4{x*w,y*w,z*w,w}; }
            surfs[fid].push_back(s);
            std::string e; f>>e;
        } else if(tok=="CURVE"){
            int eid,si,deg,n; f>>eid>>si>>deg>>n;
            Curve c; c.deg=deg; c.P.resize(n);
            for(int i=0;i<n;++i){ double x,y,z,w; f>>x>>y>>z>>w; c.P[i]=H4{x*w,y*w,z*w,w}; }
            curves[eid].push_back(c);
            std::string e; f>>e;
        } else if(tok=="NPAIR"){int n;f>>n;pairs.resize(n);}
        else if(tok=="PAIR"){
            static size_t pi=0; int e,s,m; f>>e>>s>>m;
            pairs[pi].e=e; pairs[pi].s=s;
            for(int i=0;i<m;++i){ double x,y,z; f>>x>>y>>z; pairs[pi].pts.push_back({x,y,z}); }
            std::string en; f>>en; ++pi;
        }
    }
    int nexp=0; for(auto&p:pairs) nexp+=(int)p.pts.size();
    std::printf("parsed: curves=%zu surfs=%zu pairs=%zu expected=%d\n",curves.size(),surfs.size(),pairs.size(),nexp);
    auto t0=std::chrono::high_resolution_clock::now();
    int matched=0,missing=0,extra=0,uniqTotal=0;
    for(auto&p:pairs){
        std::vector<std::array<double,3>> roots;
        for(auto&c:curves[p.e]) for(auto&s:surfs[p.s]){
            for(int iu=0;iu<=6;++iu) for(int is=0;is<=6;++is) for(int it=0;it<=6;++it){
                double pt3[3];
                if(certify(c,s,iu/6.0,is/6.0,it/6.0,pt3)){
                    bool dup=false;
                    for(auto&q:roots){ double d=std::sqrt((pt3[0]-q[0])*(pt3[0]-q[0])+(pt3[1]-q[1])*(pt3[1]-q[1])+(pt3[2]-q[2])*(pt3[2]-q[2])); if(d<1e-6){dup=true;break;} }
                    if(!dup) roots.push_back({pt3[0],pt3[1],pt3[2]});
                }
            }
        }
        std::vector<bool> usedE(p.pts.size(),false), usedG(roots.size(),false);
        for(size_t i=0;i<p.pts.size();++i) for(size_t j=0;j<roots.size();++j){
            if(usedG[j]) continue;
            double d=std::sqrt((p.pts[i][0]-roots[j][0])*(p.pts[i][0]-roots[j][0])+(p.pts[i][1]-roots[j][1])*(p.pts[i][1]-roots[j][1])+(p.pts[i][2]-roots[j][2])*(p.pts[i][2]-roots[j][2]));
            if(d<1e-4){ usedE[i]=true; usedG[j]=true; ++matched; break; }
        }
        int miss=0, ext=0;
        for(size_t i=0;i<p.pts.size();++i) if(!usedE[i]) ++miss;
        for(size_t j=0;j<roots.size();++j) if(!usedG[j]) ++ext;
        missing+=miss; extra+=ext; uniqTotal+=(int)roots.size();
        std::printf("  pair(edge %d, face %d): expected %zu, certified %zu, miss %d, extra %d\n",p.e,p.s,p.pts.size(),roots.size(),miss,ext);
    }
    auto t1=std::chrono::high_resolution_clock::now();
    double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    std::printf("\nsummary: pairs=%zu expected=%d certified_unique=%d matched=%d missing=%d extra=%d starts_certified=%d time=%.1f ms\n",
        pairs.size(),nexp,uniqTotal,matched,missing,extra,g_cert,ms);
    return 0;
}
