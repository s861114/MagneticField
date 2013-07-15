function[output]=nonzeroopdc_mag(indc,E,V,B,v0,eta,num_bands,hw,n_max,Ef)
    
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
    
    sum=0;
    
    for n=1:num_bands
        if E(n)<=Ef
%         fn=1./(1.+exp((E(n)-Ef+0.00000000000000000000000001)*inf));        
            sum=sum+core(n);
        end
    end
    output=sum*const;
    
%-------------------------inline functions---------------------------------
function[output_in]=core(n)
    output_in=0;
    for m=1:num_bands
        if m~=n
            cntr1_f=V(:,n)'*Vx*V(:,m);
            cntr2_f=V(:,m)'*Vy*V(:,n);
             
%             cntr1_f=contraction(V(:,n),Vx,V(:,m));
%             cntr2_f=contraction(V(:,m),Vy,V(:,n));  
            cntr1_b=V(:,m)'*Vx*V(:,n);
            cntr2_b=V(:,n)'*Vy*V(:,m);  
            
            
            denomf=( E(n)-E(m) )*( E(n)-E(m)+ hw+1i*eta);
            denomb=( E(n)-E(m) )*( E(n)-E(m)- hw-1i*eta);
            sum_temp=cntr1_f*cntr2_f/denomf-cntr1_b*cntr2_b/denomb;
            output_in=output_in+1i*sum_temp;
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
            
            outputSpin(1)=1;
            for i=1:n_max
                outputSpin(2*i+1)=1;
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