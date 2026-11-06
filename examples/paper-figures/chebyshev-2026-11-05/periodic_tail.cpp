#include <vector>
#include <algorithm>
#include <cstdint>
// Infer a tail only after 32,768 iterates; require >=8 repeats in the final 8192 symbols.
template<class T> std::vector<unsigned char> sequence(double ud,double vd,int r,int N){
 T u=ud,v=vd,x=r?v:u;int s=r?1:-1;std::vector<unsigned char>a(N+1);a[0]=0;
 for(int n=1;n<=N;n++){int lap=x<T(-.5)?0:x>T(.5)?2:(x==T(-.5)?(s<0?0:1):(x==T(.5)?(s<0?1:2):1));a[n]=lap+(s>0?3:0);s*=lap==1?-1:1;x=((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2;}return a;
}
std::pair<int,int> period(const std::vector<unsigned char>&s){
 int N=s.size()-1,L=8192,start=N-L+1;std::vector<int>pi(L);
 for(int i=1;i<L;i++){int j=pi[i-1];while(j&&s[start+i]!=s[start+j])j=pi[j-1];if(s[start+i]==s[start+j])j++;pi[i]=j;}
 int p=L-pi.back();if(p>1024)return {-1,-1};
 int begin=N-p+1;while(begin>1&&s[begin-1]==s[begin-1+p])begin--;
 return {begin,p};
}
extern "C" int periodic_numerator(double u,double v,int N,int capacity,long long*out,int*meta){
 auto s0=sequence<double>(u,v,0,N),s1=sequence<double>(u,v,1,N);auto [m0,p0]=period(s0);auto [m1,p1]=period(s1);
 if(p0<0||p1<0)return -1;
 // Independent long-double itinerary must give the same rational kneading rows.
 auto t0=sequence<long double>(u,v,0,N),t1=sequence<long double>(u,v,1,N);auto [mt0,pt0]=period(t0);auto [mt1,pt1]=period(t1);
 if(pt0!=p0||pt1!=p1||mt0!=m0||mt1!=m1)return -2;
 for(int i=1;i<m0+p0;i++)if(s0[i]!=t0[i])return -2;
 for(int i=1;i<m1+p1;i++)if(s1[i]!=t1[i])return -2;
 int n0=m0+p0,n1=m1+p1;if(n0+n1-1>capacity)return -3;
 std::vector<long long>a(n0),b(n0),c(n1),d(n1);a[0]=1;c[0]=-1;d[0]=1;
 for(int i=1;i<n0;i++){int l=s0[i]%3,sg=s0[i]>=3?2:-2;if(l==1)a[i]=sg;if(l==2)b[i]=sg;}
 for(int i=1;i<n1;i++){int l=s1[i]%3,sg=s1[i]>=3?2:-2;if(l==1)c[i]=sg;if(l==2)d[i]=sg;}
 for(int i=n0-1;i>=p0;i--){a[i]-=a[i-p0];b[i]-=b[i-p0];}
 for(int i=n1-1;i>=p1;i--){c[i]-=c[i-p1];d[i]-=d[i-p1];}
 std::fill(out,out+n0+n1-1,0);for(int i=0;i<n0;i++)for(int j=0;j<n1;j++)out[i+j]+=a[i]*d[j]-b[i]*c[j];
 meta[0]=m0;meta[1]=p0;meta[2]=m1;meta[3]=p1;
 int n=n0+n1-1;while(n>1&&!out[n-1])n--;return n;
}
