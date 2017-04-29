function [L,U,P]=LU_pivot(A)


[n,n]=size(A);
L=eye(n); P=L; U=A;
for k=1:n
    %fprintf('col %d\n',k);
    [pivot m]=max(abs(U(k:n,k)));%find the pivot
    m=m+k-1;
    if m~=k
        % interchange rows m and k in U
        temp=U(k,:);
        U(k,:)=U(m,:);
        U(m,:)=temp;
        % interchange rows m and k in P
        temp=P(k,:);
        P(k,:)=P(m,:);
        P(m,:)=temp;
        if k >= 2
            temp=L(k,1:k-1);
            L(k,1:k-1)=L(m,1:k-1);
            L(m,1:k-1)=temp;
        end
    end
    for j=k+1:n
        if(U(j,k)~=0)%ignore calculation of 0 increse speed
        L(j,k)=U(j,k)/U(k,k);
        U(j,:)=U(j,:)-L(j,k)*U(k,:);
        else
            continue
        end
    end
    
end