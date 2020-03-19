function Aq = matvec(q)

% Multiply A * x for the elliptic problem A*q approx q''(x)

Aq = zeros(size(q));

n = length(q);

dx2 = dx*dx;

for i = 1:n
    if (i == 1)
        Aq(i) =  (-2*q(1) + q(2))/dx2;
    elseif (i == n)
        Aq(i) = (q(n-1) - 2*q(n))/dx2;
    else
        Aq(i) = (q(i-1) - 2*q(n) + x(i+1))/dx2;
    end
end







end