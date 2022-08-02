
A = diag(4.1*ones(1,4)) + diag(2*ones(1,4-1),1) + diag(2*ones(1,4-1),-1)
x = sor(A, ones(4,1), 10, 1)

function r = sor(A, b, n, w)
    D = diag(diag(A));
    L = tril(A)-D;
    U = triu(A)-D;
    r = zeros(1,n);
    x = zeros(size(b));
    for k = 1:n
        x = (L+(w*D))\(b - (U+(1-w)*D)*x);
        r(k) = norm(b - A*x);
    end
end

