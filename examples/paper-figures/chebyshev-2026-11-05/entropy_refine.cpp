#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
using P=std::vector<long long>;
using B=std::vector<long double>;
using boost::multiprecision::cpp_int;
using MP=boost::multiprecision::cpp_dec_float_100;
static std::vector<B> weights;
extern "C" void refine_init(int n){
 weights.clear(); weights.resize(n+1); for(int k=0;k<=n;k++) weights[k].resize(k+1);
 for(int k=0;k<=n;k++){weights[k][0]=1;for(int i=1;i<=k;i++)weights[k][i]=weights[k][i-1]*(k-i+1)/(n-i+1);}
}
extern "C" void refine_coeff(double u,double v,int N,long long*out){
 std::vector<int> a(N+1),b(N+1),c(N+1),d(N+1);a[0]=1;c[0]=-1;d[0]=1;
 for(int r=0;r<2;r++){
  double x=r?v:u;int s=r?1:-1;
  for(int n=1;n<=N;n++){
   int lap=x<-.5?0:x>.5?2:(x==-.5?(s<0?0:1):(x==.5?(s<0?1:2):1));
   if(lap==1)(r?c:a)[n]=2*s;if(lap==2)(r?d:b)[n]=2*s;
   s*=lap==1?-1:1;x=((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2;
  }
 }
 long long prev=0;
 for(int n=0;n<=N;n++){long long z=0;for(int i=0;i<=n;i++)z+=a[i]*d[n-i]-b[i]*c[n-i];out[n]=prev+=z;}
}
long double val(const P&p,long double x,long double&bound){
 long double y=p.back(),a=std::abs(y);for(int i=int(p.size())-2;i>=0;i--){y=y*x+p[i];a=a*x+std::abs((long double)p[i]);}bound=a*(p.size()+2)*2e-19L;return y;
}
int precise_sign(const P&p,long double x){MP a=x,y=p.back();for(int i=int(p.size())-2;i>=0;i--)y=y*a+p[i];return (y>0)-(y<0);}
// Exact sign of the Bernstein coefficient; used when cancellation obscures it.
long double exact_bernstein(const P&p,int k){
 int n=p.size()-1;cpp_int bin=1;
 for(int j=1;j<=k;j++)bin=bin*(n-j+1)/j;
 cpp_int denom=bin,sum=cpp_int(p[0])*bin;
 for(int i=1;i<=k;i++){bin=bin*(k-i+1)/(n-i+1);sum+=cpp_int(p[i])*bin;}
 if(sum==0)return 0;
 return (MP(sum)/MP(denom)).convert_to<long double>();
}
long double search(const P&p,const B&b,long double lo,long double hi,int depth){
 int n=p.size()-1,prev=0,v=0;for(auto z:b){int s=(z>0)-(z<0);if(s){v+=prev&&prev!=s;prev=s;}}
 if(!v)return -1;
 if(v==1){
  long double err;auto fl=val(p,lo,err);int sl=std::abs(fl)>err?((fl>0)-(fl<0)):precise_sign(p,lo);
  auto fh=val(p,hi,err);int sh=std::abs(fh)>err?((fh>0)-(fh<0)):precise_sign(p,hi);
  if(!sl&&lo>0)return lo;if(!sh&&hi<1)return hi;if(sl==sh)return NAN;
  for(int j=0;j<60&&hi-lo>1e-11L*std::max(lo,.1L);j++){
   auto mid=(lo+hi)/2,f=val(p,mid,err);int s=std::abs(f)>err?((f>0)-(f<0)):precise_sign(p,mid);
   if(!s)return mid;if(s==sl)lo=mid;else hi=mid;
  }return (lo+hi)/2;
 }
 if(depth>45)return NAN;
 B work=b,left(n+1),right(n+1);left[0]=b[0];right[n]=b[n];
 for(int j=1;j<=n;j++){for(int k=0;k<=n-j;k++)work[k]=(work[k]+work[k+1])/2;left[j]=work[0];right[n-j]=work[n-j];}
 auto mid=(lo+hi)/2,r=search(p,left,lo,mid,depth+1);if(r!=-1)return r;
 if(left[n]==0)return mid;
 return search(p,right,mid,hi,depth+1);
}
extern "C" double refine_entropy(const long long*ptr,int n){
 P p(ptr,ptr+n+1);while(p.size()>1&&!p.back())p.pop_back();
 while(p.size()>1){long long sum=0;for(auto a:p)sum+=a;if(sum)break;for(size_t i=1;i<p.size();i++)p[i]+=p[i-1];p.pop_back();}
 n=p.size()-1;if(!n)return 0;
 B b(n+1);
 bool cached=weights.size()==p.size();
 for(int k=0;k<=n;k++){
  long double z=p[0],a=std::abs(z),w=1;
  for(int i=1;i<=k;i++){if(cached)w=weights[k][i];else w=w*(k-i+1)/(n-i+1);z+=p[i]*w;a+=std::abs(p[i]*w);}
  if(std::abs(z)<=a*(n+2)*2e-19L)z=exact_bernstein(p,k);
  b[k]=z;
 }
 auto r=search(p,b,0,1,0);return r<0?0:-std::log((double)r);
}

extern "C" void refine_signature(double u,double v,int N,unsigned char*out){
 for(int r=0;r<2;r++){double x=r?v:u;int s=r?1:-1;
  for(int n=0;n<N;n++){int lap=x<-.5?0:x>.5?2:(x==-.5?(s<0?0:1):(x==.5?(s<0?1:2):1));
   out[r*N+n]=lap+(s>0?3:0);s*=lap==1?-1:1;x=((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2;
  }
 }
}
