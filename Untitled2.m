clear all;

n_max=10;
indc=2;

N=wksp.Nband(indc);
N_layer=wksp.Nlayer(indc);
N_size=wksp.size_H(indc); %determine the size of the part of a hamiltonian for n=0

num_bands=N_size+2*N_layer*n_max;

gamma0=3;
gamma1=0.3;

B=30;
eta=0;
hw=0;
v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);
unit=4*wksp.e^2/wksp.h;
[E,V]=DiagH_SC_Mag2(gamma0,gamma1,indc,n_max,B);   

EfS=-0.5:0.001:0.5;
cnt=0;

for Ef=EfS   
    tic
    cnt=cnt+1;
    y(cnt)=nonzeroopdc_mag(indc,E,V,B,v0,eta,num_bands,hw,n_max,Ef);
    toc
end

plot(EfS/gamma0,y/unit)