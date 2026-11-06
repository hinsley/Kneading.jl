#define main entropy_cli_main
#include "entropy_boundary.cpp"
extern "C" void entropy_init(){init();}
extern "C" void entropy_coefficients(double u,double v,int*out){auto p=determinant(u,v);for(int i=0;i<65;i++)out[i]=p[i];}
extern "C" double entropy_value(int*in,int n){Poly p;for(int i=0;i<65;i++)p[i]=in[i];return entropy(p,n);}
