clear all;

eta=0/1000;
n_max=30;
BS=0:0.01:10;
indc=1;
N=wksp.Nband(indc);
N_layer=wksp.Nlayer(indc);
N_size=wksp.size_H(indc); %determine the size of the part of a hamiltonian for n=0
E=zeros(N_size+2*N_layer*n_max,size(BS,2));
V=zeros(N_size+2*N_layer*n_max,N_size+2*N_layer*n_max,size(BS,2));
A=zeros(1,10);

gamma1=0.3;
gamma0=3;
cnt=0;

for B=BS
    tic
    cnt=cnt+1;
    [E(:,cnt),V(:,:,cnt)]=DiagH_SC_Mag2(gamma0,gamma1,indc,n_max,B);    
    toc
end

E_real=E(:,1);
V_real=V(:,:,1);

for bndidx=1:N_layer+2*N_layer*n_max
    plot(BS,squeeze(E(bndidx,:))/gamma0);
    hold on;
end
pause
hw=0;
EfS=-0.2:0.001:0.2;
v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);
cnt=0;
unit=4*wksp.e^2/wksp.h;
for Ef=EfS
    tic
    cnt=cnt+1;
    EfS(cnt)
    y(cnt)=nonzeroopdc_mag(indc,E_real,V_real,B,v0,eta,N_size+2*N_layer*n_max,hw,n_max,Ef)/unit;
    toc
end

plot(EfS/gamma0,y)




% v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);
% graph=zeros(size(BS,2),1);
% graph2=zeros(size(BS,2),1);
% for r=1:N/2
%     for n=1:n_max
%         cnt=0;
%         for B=BS
%             cnt=cnt+1;
%             length=sqrt(wksp.hbar/(wksp.e)/B);        
%             graph(cnt)=(sqrt(2*n)*wksp.hevbar*v0/length+2*gamma1*cos(r*pi/(N/2+1)))/gamma0;
%             graph2(cnt)=(-sqrt(2*n)*wksp.hevbar*v0/length+2*gamma1*cos(r*pi/(N/2+1)))/gamma0;
%         end    
%         hold on;
%         plot(BS,graph,'-r');    
%         plot(BS,graph2,'-r');    
%     end
% end
% axis([0 20 -0.04 0.04]);