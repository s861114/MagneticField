function[output]=nonzeroopdc_mag_xx(indc,E,V,B,v0,eta,num_bands,hw,n_max,Ef)
    
    [HO_state,HO_spin]=HO_func(indc);
    N=wksp.Nband(indc);
    hevbar=wksp.hevbar;    
    [Vx,Vy]=vlcty_op_create(indc);    
%     Vx
%     imag(Vy)
%     pause

    length=sqrt(wksp.hbar/(wksp.e)/B);
    dgnrcy=2*2; %spin+valley degeneracy    
    const=-dgnrcy*wksp.e^2/(2*pi*wksp.hbar*length^2);  %- 부호 다름(xy인 경우)    
    sumR=0;
    sumI=0;
    for n=1:num_bands
        if E(n)<Ef
            for m=1:num_bands            
                [tempR,tempI]=core(n);
                sumR=sumR+tempR;
                sumI=sumI+tempI;
                if E(n)==Ef
                    sumR=sumR-tempR/2;
                    sumI=sumI-tempI/2;
                end
            end
        end
    end
    output=(sumR+1i*sumI)*const;
    
%-------------------------inline functions---------------------------------
function[output_R, output_I]=core(n)
    output_R=0;
    output_I=0;
    for m=1:num_bands
        if m~=n
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
    end
end

function[output,outputSpin]=HO_func(indc)
    output=zeros(num_bands,1);
    outputSpin=zeros(num_bands,1);
    switch indc
        
        case 1
            output(1)=0;
            for i=1:n_max
                output(1+2*i-1)=i-1;
                output(1+2*i)=i;
            end
            
            outputSpin(1)=1;
            for i=1:n_max
                outputSpin(2*i+1)=1;
            end
            
        case 2
            
            output(1:3)=[0 0 1];
            for cntidx=1:n_max
                output(3+4*cntidx-3)=cntidx-1;
                output(3+4*cntidx-2)=cntidx;
                output(3+4*cntidx-1)=cntidx;
                output(3+4*cntidx)=cntidx+1;
            end
            
            outputSpin(1:3)=[1 -2 2];
            cnt=0;            
            while(cnt<n_max)
                cnt=cnt+1;                
                outputSpin(4*cnt:4*cnt+3)=[-1 1 -2 2];                
            end
            
        case 3
        case 4
        case 5            
    end
end

function[outputVx, outputVy]=vlcty_op_create(indc)
    switch indc
        case 1
            outputVx=zeros(num_bands);
            outputVy=zeros(num_bands);
           

            for ii=1:num_bands
                for jj=ii+1:num_bands
                    if HO_spin(ii)~=HO_spin(jj)
                        if HO_state(ii)==HO_state(jj)
                            outputVx(ii,jj)=1;
                            outputVx(jj,ii)=1;
                            
                            outputVy(ii,jj)=-1i;
                            outputVy(jj,ii)=1i;
                        end
                    end
                end
            end
            
            outputVx=outputVx*v0*hevbar;
            outputVy=outputVy*v0*hevbar;
        case 2
            outputVx=zeros(num_bands);
            outputVy=zeros(num_bands);
           

            for ii=1:num_bands
                for jj=ii+1:num_bands
                    if HO_spin(ii)==-HO_spin(jj) 
                        if HO_state(ii)==HO_state(jj)
                            outputVx(ii,jj)=1;
                            outputVx(jj,ii)=1;
                            
                            outputVy(ii,jj)=-1i;
                            outputVy(jj,ii)=1i;
                        end
                    end
                end
            end
            
            outputVx=outputVx*v0*hevbar;
            outputVy=outputVy*v0*hevbar;

        case 3
        case 4
        case 5
    end
end

function[output]=contraction(V1,V_vel,V2)
    output=0;
    for i=1:num_bands
        for j=1:num_bands
            if HO_state(i)==HO_state(j)
                output=output+V1(i)*V_vel(i,j)*V2(j);
            end
        end
    end

end
        
    
end