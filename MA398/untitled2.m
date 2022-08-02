y = [];
for n = 1:100
A = diag(4.1*ones(1,n)) + diag(2*ones(1,n-1),1) + diag(2*ones(1,n-1),-1);
b = ones(n,1);
x = A\b;
y(end+1) = norm(x, inf);
end
x = linspace(1,100);
plot(x,y)


jaco = [];
GS = [];
for n = 128:992:4096
A = diag(4.1*ones(1,n)) + diag(2*ones(1,n-1),1) + diag(2*ones(1,n-1),-1);
D = diag(4.1*ones(1,n));
L = diag(2*ones(1,n-1),-1);
U = diag(2*ones(1,n-1),1);
b = ones(n,1);
r = ones(n,1);
x = zeros(n,1);
counter = 0;
while norm(r, inf)> 10^-9
    [r, x]= jacobi(D, L, U, b, x);
    counter= counter+1;
end
jaco(end+1)=counter
counter = 0; r = ones(n,1); x = zeros(n,1);
while norm(r, inf)>= 10^-9
    [r, x]= GaussSeidel(D, L, U, b, x);
    counter= counter+1;
end
GS(end+1)=counter
end
x1 = linspace(128,4096,5);
plot(x1, jaco)
    
function [r, x] = jacobi(D, L, U, b, x)
x = D\(b - (U+L)*x);
r = b - (L+D+U)*x;
end

function [r, x] = GaussSeidel(D, L, U, b, x)
x = (L+D)\(b - U*x);
r = b - (L+D+U)*x;
end