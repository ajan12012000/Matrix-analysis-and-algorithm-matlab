y = [];
for n = 1:100
A = diag(4.1*ones(1,n)) + diag(2*ones(1,n-1),1) + diag(2*ones(1,n-1),-1);
b = ones(n,1);
x = A\b;
y(end+1) = norm(x, inf);
end
x = linspace(1,100);
plot(x,y)
