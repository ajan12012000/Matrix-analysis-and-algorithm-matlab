jacobianlist = [];
sorlist = [];
n = 256;
A = diag(4.1*ones(1,n)) + diag(2*ones(1,n-1),1) + diag(2*ones(1,n-1),-1);
D = diag(4.1*ones(1,n));
L = diag(2*ones(1,n-1),-1);
U = diag(2*ones(1,n-1),1);
b = ones(n,1);
r = ones(n,1);
x = zeros(n,1);
counter = 0;
for i = 0:3
    while norm(r, inf)> 10^-(2^i)
        [r, x]= jacobi(D, L, U, b, x);
        counter= counter+1;
    end
    jacobianlist(end+1) = counter;
    counter = 0; r = ones(n,1); x = zeros(n,1);
    while norm(r, inf)> 10^-(2^i)
        [r, x]= sor(D, L, U, b, x, 0.6);
        counter= counter+1;
    end
    sorlist(end+1) = counter;
end
jacobianlist
sorlist
function [r, x] = jacobi(D, L, U, b, x)
x = D\(b - (U+L)*x);
r = b - (L+D+U)*x;
end

function [r, x] = sor(D, L, U, b, x, w)
x = (L+(w*D))\(b - (U+(1-w)*D)*x);
r = b - (L+D+U)*x;
end