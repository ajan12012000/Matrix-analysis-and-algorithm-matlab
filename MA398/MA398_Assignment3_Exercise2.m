
A = diag(1:1:100) + diag(ones(1,99),1) + diag(ones(1,99),-1);
b = ones(100,1);
k2 = norm(A)*norm(inv(A));
[x, r] = SD(A,b);
x
residualnormSD = norm(r, 2)
upperboundSD = (sqrt(1-1/k2))^100
[x,flag,relres] = pcg(A, b);
x
residualnormCG = relres
upperboundCG = 2*inv((sqrt(k2)+1)/(sqrt(k2)-1)+(sqrt(k2)-1)/(sqrt(k2)+1))


function [x, r] = SD(A, b)
x = zeros(100,1);
for k = 1:100
    r = b - A*x;
    alpha=norm(r,2)^2/dot(r,A*r);
    x= x + alpha*r;
end
end
  