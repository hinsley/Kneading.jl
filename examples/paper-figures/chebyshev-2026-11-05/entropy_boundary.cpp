// Degree-64 finite kneading determinant, matching the site's documented method.
// This region has u>=v and invariant interval [-1,1], so no escape continuation is needed.
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <omp.h>
using Poly=std::array<int,65>;
using Wide=__int128_t;
Wide choose[65][65];
void init(){for(int n=0;n<=64;n++){choose[n][0]=choose[n][n]=1;for(int k=1;k<n;k++)choose[n][k]=choose[n-1][k-1]+choose[n-1][k];}}
Poly determinant(double u,double v){
 int rows[2][3][65]={};rows[0][0][0]=-1;rows[0][1][0]=1;rows[1][1][0]=-1;rows[1][2][0]=1;
 for(int r=0;r<2;r++){
  double x=r?v:u;int s=r?1:-1;
  for(int n=1;n<=64;n++){
   int lap=x<-.5?0:x>.5?2:(x==-.5?(s<0?0:1):(x==.5?(s<0?1:2):1));
   rows[r][lap][n]=2*s;s*=lap==1?-1:1;
   x=((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2;
  }
 }
 Poly p={};for(int n=0;n<=64;n++){
  for(int i=0;i<=n;i++)p[n]+=rows[0][1][i]*rows[1][2][n-i]-rows[0][2][i]*rows[1][1][n-i];
  if(n)p[n]+=p[n-1]; // formal division by 1-t, truncated through degree64
 }return p;
}
long double eval(const Poly&p,int n,long double t,long double &err){
 long double x=p[n],a=std::abs(x);for(int i=n-1;i>=0;i--){x=x*t+p[i];a=a*t+std::abs(p[i]);}err=a*2e-17L;return x;
}
// Search left to right using Bernstein sign variations and de Casteljau splits.
// Ambiguous/multiple-root cases are returned as NaN for independent Arb isolation.
long double search(const Poly&p,int n,const std::array<long double,65>&b,long double lo,long double hi,int depth){
 int prev=0,variations=0;for(int i=0;i<=n;i++){int s=(b[i]>0)-(b[i]<0);if(s){variations+=prev&&s!=prev;prev=s;}}
 if(!variations)return -1;
 if(variations==1){
  long double e,flo=eval(p,n,lo,e),fhi=eval(p,n,hi,e);
  if(flo==0&&lo>0)return lo;if(fhi==0&&hi<1)return hi;
  if((flo>0)==(fhi>0))return NAN;
  for(int i=0;i<64&&hi-lo>1e-12L*std::max(lo,0.1L);i++){
   long double mid=(lo+hi)/2,f=eval(p,n,mid,e);
   if(std::abs(f)<=e)return NAN;
   if((f>0)==(flo>0)){lo=mid;flo=f;}else hi=mid;
  }return (lo+hi)/2;
 }
 if(depth>45)return NAN;
 auto work=b;std::array<long double,65>left{},right{};left[0]=b[0];right[n]=b[n];
 for(int j=1;j<=n;j++){for(int k=0;k<=n-j;k++)work[k]=(work[k]+work[k+1])/2;left[j]=work[0];right[n-j]=work[n-j];}
 long double mid=(lo+hi)/2;
 if(left[n]==0&&mid>0){ // exact subdivision root; search strictly left first
  auto result=search(p,n,left,lo,mid,depth+1);return result==-1?mid:result;
 }
 auto result=search(p,n,left,lo,mid,depth+1);if(result!=-1)return result;
 return search(p,n,right,mid,hi,depth+1);
}
double entropy(Poly p,int n){
 while(n&&p[n]==0)n--;
 // Remove exact roots at t=1, which are outside the open search interval.
 while(n){long long sum=0;for(int i=0;i<=n;i++)sum+=p[i];if(sum)break;
  for(int i=1;i<=n;i++)p[i]+=p[i-1];n--;
 }
 if(!n)return 0;
 std::array<long double,65>b{};
 for(int k=0;k<=n;k++){Wide sum=0;for(int i=0;i<=k;i++)sum+=Wide(p[i])*choose[n-i][k-i];b[k]=(long double)sum/(long double)choose[n][k];}
 auto r=search(p,n,b,0,1,0);return r<0?0:-std::log((double)r);
}
int main(int argc,char**argv){
 init();if(argc==4){double u=std::stod(argv[2]),v=std::stod(argv[3]);auto p=determinant(u,v);for(int a:p)std::cout<<a<<' ';std::cout<<'\n'<<std::setprecision(17)<<entropy(p,64)<<' '<<entropy(p,32)<<'\n';return 0;}
 int N=argc>1?std::stoi(argv[1]):3000;std::vector<double>us(N),vs(N),h64(N*N),h32(N*N);
 std::ifstream axes("entropy_axes_"+std::to_string(N)+".bin",std::ios::binary);axes.read((char*)us.data(),N*8);axes.read((char*)vs.data(),N*8);if(!axes)return 2;
 std::atomic<int>done=0,bad=0;double start=omp_get_wtime();
 #pragma omp parallel for schedule(dynamic,1)
 for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){long k=long(i)*N+j;if(us[i]==vs[j]){h64[k]=h32[k]=0;continue;}
   auto p=determinant(us[i],vs[j]);h64[k]=entropy(p,64);h32[k]=entropy(p,32);if(!std::isfinite(h64[k])||!std::isfinite(h32[k]))bad++;
  }
  int n=++done;if(n%100==0||n==N){
   #pragma omp critical
   std::cout<<n<<'/'<<N<<" columns, "<<omp_get_wtime()-start<<" s; ambiguous "<<bad<<std::endl;
  }
 }
 std::ofstream out("entropy_raw_"+std::to_string(N)+".bin",std::ios::binary);out.write((char*)h64.data(),h64.size()*8);out.write((char*)h32.data(),h32.size()*8);
}
