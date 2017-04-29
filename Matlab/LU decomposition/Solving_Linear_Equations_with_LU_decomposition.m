function x = Solving_Linear_Equations_with_LU_decomposition(A, b)

 [L,U,P]=LU_pivot(A);
 
 newb=P*b;
 %Ld=b slove d
 [n,n] = size(L);
 d=zeros(n,1);
 d(1,1) = newb(1,1)/L(1,1);
 for k1 = 2:n
    j = 1:k1-1;
    d(k1,1) = (newb(k1,1) - L(k1,j)*d(j,1))/L(k1,k1);
 end
%Ux=d slove x
 x = zeros(n,1);
 x(n,1) = d(n,1)/U(n,n);
 for k2 = 1:n-1
    j = n-k2;
    x(j,1) = (d(j,1) - U(j,j+1:n)*x(j+1:n,1))/U(j,j);
 end

% x=A\b; 
end

