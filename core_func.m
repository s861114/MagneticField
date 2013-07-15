function[output_R, output_I]=core_func(n,m,V,Vx,eta,hw,E)
    output_R=0;
    output_I=0;
        
    cntr=V(:,n)'*Vx*V(:,m);             
    temp1=E(n)-E(m);
    temp2=( (temp1+hw)^2+eta^2 ) ;
    deno=temp1*temp2;
    front_r=eta/( deno );
    front_i=(temp1+hw)/( deno );

    %1i/denomb
    temp2=( (temp1-hw)^2+eta^2 );
    deno=temp1*temp2;
    back_r=-eta/( deno );
    back_i=(temp1-hw)/( deno );

%             denomf=( E(n)-E(m) )*( E(n)-E(m)+ hw+1i*eta);
%             denomb=( E(n)-E(m) )*( E(n)-E(m)- hw-1i*eta);
    A=cntr^2;            
    output_R=output_R+A*front_r-A*back_r;
    output_I=output_I+A*front_i-A*back_i;    
end