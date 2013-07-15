clear all;

n_max=30;
B=28;


indc=2;
N=wksp.Nband(indc);
N_layer=wksp.Nlayer(indc);
N_size=wksp.size_H(indc); %determine the size of the part of a hamiltonian for n=0

E=zeros(N_size+2*N_layer*n_max);
V=zeros(N_size+2*N_layer*n_max,N_size+2*N_layer*n_max);
A=zeros(1,10);

gamma1=0.39;
gamma0=3.16;
eta=10/1000;
hwS=0:0.001:gamma1*2.5;
% hwS=0.39:0.001:0.468;

[E,V]=DiagH_SC_Mag2(gamma0,gamma1,indc,n_max,B);    

v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);
cnt=0;
unit=4*wksp.e^2/wksp.hbar;
Ef=0;
for cnt=1:size(hwS,2);
    hw=hwS(cnt);
    tic
    [y(cnt)]=nonzeroopdc_mag_xx(indc,E,V,B,v0,eta,N_size+2*N_layer*n_max,hw,n_max,Ef)/unit;
    toc
end
% graph_draw


plot(hwS/gamma1,real(y))
hold on;
plot(hwS/gamma1,imag(y),'-r')


