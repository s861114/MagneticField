clear all;


n_max=10;
BS=0:0.01:20;
indc=2;
N=wksp.Nband(indc);

E=zeros(N,n_max,size(BS,2));
V=zeros(N,N,n_max,size(BS,2));
A=zeros(1,10);

gamma1=0.3;
gamma0=3;
cnt=0;
for B=BS
    cnt=cnt+1;
    [E(:,:,cnt),V(:,:,:,cnt)]=DiagH_SC_Mag2(gamma0,gamma1,indc,n_max,B);    
end

for i=1:n_max
    for bndidx=1:N
        plot(BS,squeeze(E(bndidx,i,:))/gamma0);
        hold on;
    end
end




v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);


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
axis([0 20 -0.04 0.04]);