function [ No ] = No_Const( Smaster )
%Finds number of contraints given the eigenvalue matrix produced by IPCA
%for all the runs
thresh = 0.8;
[m,n] = size(Smaster);
check = 0;
No = [];
iter = 10;
while (iter>1)
        if(thresh>0.60)
        thresh = thresh-0.01;
        else
            break
        end
        No = [];
       iter = 0;
    
    for i = m:-1:1
    check = 0;
    for j=n:-1:1
        if(abs(Smaster(i,j)-1)< thresh)
        check = check +1;
        end
    end
    if (check == (i))
        iter = iter +1;
        No(iter) = i;
        
    end
    
    end
end
if(iter>1)
    for i = 1:length(No)
       M(i) = mean(abs(Smaster(No(i),(m-No(i)+1):end)-1),2);
    end
  No =   No(abs(M)==min(abs(M)));
end


end

