function [ output ] = check( R,value,num_bands )
    output=0;
    for i=1:num_bands
        for j=1:num_bands
            if R(i,j)~=0
                if (abs(R(i,j))-abs(value))<1
                    output=1;
                    return;
                end
            end
        end
    end


end

