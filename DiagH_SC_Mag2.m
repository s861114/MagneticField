function [outputE,outputV] = DiagH_SC_Mag2(gamma0,gamma1,indc,n_max,B)
                                              

v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);

N=wksp.Nband(indc);
N_layer=wksp.Nlayer(indc);
N_size=wksp.size_H(indc);

H=zeros(N_size+2*N_layer*n_max,N_size+2*N_layer*n_max);
outputE=zeros(N_size+2*N_layer*n_max,1);
% 
% outputE=zeros(N,n_max);
% outputV=zeros(N,N,n_max);
%     

length=sqrt(wksp.hbar/(wksp.e)/B);


H_0=zero_func(indc);
H(1:N_size,1:N_size)=H_0;

for n=1:n_max
    switch indc
    case 1
                    Htemp=zeros(N,N);
                    Htemp(1,2)=Offdiag(n);Htemp(2,1)=Offdiag(n); %Mono  
        
    case 2
                    Htemp=zeros(N,N);
                    Htemp(1,2)=Offdiag(n);Htemp(2,1)=Offdiag(n);
                    Htemp(3,4)=Offdiag(n+1);Htemp(4,3)=Offdiag(n+1);
                    Htemp(2,3)=gamma1;Htemp(3,2)=gamma1; %Bi AB
    
    case 3
                    Htemp=zeros(N,N);
                    Htemp(1,2)=Offdiag(n);Htemp(2,1)=Offdiag(n);
                    Htemp(3,4)=Offdiag(n);Htemp(4,3)=Offdiag(n);                    
                    Htemp(2,4)=gamma1;Htemp(4,2)=gamma1;
                    Htemp(1,3)=gamma1;Htemp(3,1)=gamma1; %Bi AA

    case 4
                    Htemp=zeros(N,N);
                    Htemp(1,2)=Offdiag(n);Htemp(2,1)=Offdiag(n);
                    Htemp(3,4)=Offdiag(n+1);Htemp(4,3)=Offdiag(n+1);
                    Htemp(5,6)=Offdiag(n+2);Htemp(6,5)=Offdiag(n+2);                        

                    Htemp(2,3)=gamma1;Htemp(3,2)=gamma1;
                    Htemp(4,5)=gamma1;Htemp(5,4)=gamma1; %Tri ABC


    case 5
                    Htemp=zeros(N,N);
                    Htemp(1,2)=Offdiag(n);Htemp(2,1)=Offdiag(n);
                    Htemp(3,4)=Offdiag(n+1);Htemp(4,3)=Offdiag(n+1);
                    Htemp(5,6)=Offdiag(n);Htemp(6,5)=Offdiag(n);                

                    Htemp(2,3)=gamma1;Htemp(3,2)=gamma1;
                    Htemp(3,6)=gamma1;Htemp(6,3)=gamma1;%Tri ABA

     case 6
                    Htemp=zeros(N,N);
                    Htemp(1,2)=Offdiag(n);Htemp(2,1)=Offdiag(n);                    
                    Htemp(3,4)=Offdiag(n);Htemp(4,3)=Offdiag(n);
                    Htemp(5,6)=Offdiag(n);Htemp(6,5)=Offdiag(n);                    
                    Htemp(7,8)=Offdiag(n);Htemp(8,7)=Offdiag(n);                    
                    
                    Htemp(1,3)=gamma1;Htemp(3,1)=gamma1;
                    Htemp(2,4)=gamma1;Htemp(4,2)=gamma1;
                    
                    Htemp(3,5)=gamma1;Htemp(5,3)=gamma1;
                    Htemp(4,6)=gamma1;Htemp(6,4)=gamma1;
                    
                    Htemp(5,7)=gamma1;Htemp(7,5)=gamma1;
                    Htemp(6,8)=gamma1;Htemp(8,6)=gamma1;
                    %Tetra AAAA
    end    
    
    %( n=0 case : This code works only for AA stacked graphene including monolayer.)  
    beginidx=N_size+(n-1)*2*N_layer;
    H(beginidx+1:beginidx+2*N_layer,beginidx+1:beginidx+2*N_layer)=Htemp;   

end
    
    [V,E]=eig(H);    
    outputV(:,:)=V;

    for ii=1:N_size+2*N_layer*n_max
        outputE(ii)=E(ii,ii);
    end

%--------------------------Inline functions-------------------------------
    function[outputH_0]=zero_func(indc)
        switch indc
            case 1                
                
                
                outputH_0=zeros(N_size);
                outputH_0(1,1)=Offdiag(0);                
            case 2
                
                
                outputH_0=zeros(N_size);                
                outputH_0(1,2)=gamma1;outputH_0(2,1)=gamma1;
                outputH_0(2,3)=Offdiag(1);outputH_0(3,2)=Offdiag(1);
                
            case 4
            case 5
        end
    end
        
    function [output]=Offdiag(n)
        output=sqrt(2*n)*wksp.hevbar*v0/length;        
    end
end

    
