extern "C" void symmetric_factors(double u,int N,long long*B,long long*Q){
 B[0]=Q[0]=1;double x=u;int s=-1;
 for(int n=1;n<=N;n++){
  int lap=x<-.5?0:x>.5?2:(x==-.5?(s<0?0:1):(x==.5?(s<0?1:2):1));
  B[n]=lap==1?2*s:0;Q[n]=lap==2?2*s:(lap==0?-2*s:0);s*=lap==1?-1:1;x=u*(4*x*x*x-3*x);
 }
}
