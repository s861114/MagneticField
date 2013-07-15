function [ outputE,outputV] = DiagH_SC_Mag(gamma0,gamma1,A,indc,n_max,B)
                                              

v0=wksp.a*gamma0*sqrt(3)/(2*wksp.hevbar);

N=wksp.Nband(indc);


outputE=zeros(N,n_max);
outputV=zeros(N,N,n_max);
    

length=sqrt(wksp.hbar/(wksp.e)/B);


for n=1:n_max
    switch indc

    case 1
                    H=zeros(N,N);
                    H(1,2)=Offdiag(n);H(2,1)=Offdiag(n); %Mono  
    case 2

                    H=zeros(N,N);
                    H(1,2)=Offdiag(n);H(2,1)=Offdiag(n);
                    H(3,4)=Offdiag(n);H(4,3)=Offdiag(n);                    
                    H(2,4)=gamma1;H(4,2)=gamma1;
                    H(1,3)=gamma1;H(3,1)=gamma1; %Bi AA
        

%                       H(1,1)=A(1);
%                       H(2,2)=A(1);
%                       H(3,3)=A(2);
%                       H(4,4)=A(2);
    case 3

                    H=zeros(N,N);
                    H(1,2)=Offdiag(n);H(2,1)=Offdiag(n);
                    H(3,4)=Offdiag(n+1);H(4,3)=Offdiag(n+1);
                    H(2,3)=gamma1;H(3,2)=gamma1; %Bi AB
        

%                       H(1,1)=A(1);
%                       H(2,2)=A(1);
%                       H(3,3)=A(2);
%                       H(4,4)=A(2);

    case 4
                    H=zeros(N,N);
                    H(1,2)=Offdiag(n);H(2,1)=Offdiag(n);
                    H(3,4)=Offdiag(n+1);H(4,3)=Offdiag(n+1);
                    H(5,6)=Offdiag(n+2);H(6,5)=Offdiag(n+2);                        

                    H(2,3)=gamma1;H(3,2)=gamma1;
                    H(4,5)=gamma1;H(5,4)=gamma1; %Tri ABC

                    H(1,1)=A(1);
                    H(2,2)=A(1);
                    H(3,3)=A(2);
                    H(4,4)=A(2);
                    H(5,5)=A(3);
                    H(6,6)=A(3);


    case 5
                    H=zeros(N,N);
                    H(1,2)=Offdiag(n);H(2,1)=Offdiag(n);
                    H(3,4)=Offdiag(n+1);H(4,3)=Offdiag(n+1);
                    H(5,6)=Offdiag(n);H(6,5)=Offdiag(n);                

                    H(2,3)=gamma1;H(3,2)=gamma1;
                    H(3,6)=gamma1;H(6,3)=gamma1;%Tri ABA

                    H(1,1)=A(1);
                    H(2,2)=A(1);
                    H(3,3)=A(2);
                    H(4,4)=A(2);
                    H(5,5)=A(3);
                    H(6,6)=A(3);
     case 6
                    H=zeros(N,N);
                    H(1,2)=Offdiag(n);H(2,1)=Offdiag(n);                    
                    H(3,4)=Offdiag(n);H(4,3)=Offdiag(n);
                    H(5,6)=Offdiag(n);H(6,5)=Offdiag(n);                    
                    H(7,8)=Offdiag(n);H(8,7)=Offdiag(n);                    
                    
                    H(1,3)=gamma1;H(3,1)=gamma1;
                    H(2,4)=gamma1;H(4,2)=gamma1;
                    
                    H(3,5)=gamma1;H(5,3)=gamma1;
                    H(4,6)=gamma1;H(6,4)=gamma1;
                    
                    H(5,7)=gamma1;H(7,5)=gamma1;
                    H(6,8)=gamma1;H(8,6)=gamma1;
                    %Tetra AAAA


    end

    [V,E]=eig(H);
    outputV(:,:,n+1)=V;

    for ii=1:N
        outputE(ii,n+1)=E(ii,ii);
    end
    
end

%--------------------------Inline functions-------------------------------
    function [output]=Offdiag(n)
        output=sqrt(2*n)*wksp.hevbar*v0/length;
        
    end

end

    
