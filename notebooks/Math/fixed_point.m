function pout = fixed_point(A,b,method,tol,kmax,fignum)

global N

if (nargin < 3)
    method = 'jacobi';
end

switch method
    case 'Jacobi'
        M = diag(diag(A));
    case 'Gauss-Seidel'
        M = tril(A);
end

g = @(x) x + M\(b - A*x);


figure(fignum);
clf;

title(sprintf('%s',method),'fontsize',18);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'fontsize',16);
hold on;


xbar = A\b;

ax = xbar(1) - 2;
bx = xbar(1) + 2;
ay = xbar(2) - 2;
by = xbar(2) + 2;

N = 500;
xe = linspace(ax,bx,N);
ye = linspace(ay,by,N);
[xm,ym] = meshgrid(xe,ye);
X = [reshape(xm,1,N*N); reshape(ym,1,N*N)];

F = g(X) - X;
Fn = reshape(sum(abs(F)),N,N);
cv = linspace(min(Fn(:)),max(Fn(:)),20);
contour(xm,ym,Fn,50);

hold on;

axis([ax,bx,ay,by]);
daspect([1,1,1]);

plot(xbar(1),xbar(2),'r.','markersize',40);

fprintf('Input an initial guess : \n');
xk = ginput(1);
xk = xk(:);
hold on;
plot(xk(1),xk(2),'r.','markersize',20);
axis([ax,bx,ay,by]);
daspect([1,1,1]);


for k = 1:kmax
    
    xkp1 = g(xk);
    
    p(k) = norm(xkp1 - xk,1);
    fprintf('%4d %16.6e\n',k,p(k));    
        
    plot([xk(1) xkp1(1)],[xk(2) xkp1(2)],'r','linewidth',2);
    hold on;
    plot([xk(1) xkp1(1)],[xk(2) xkp1(2)],'k.','markersize',30);
    axis([ax,bx,ay,by]);
    daspect([1,1,1]);
    
    drawnow;

    % pause;
    
    xk = xkp1;
    if (p(k) < tol)
        break;
    elseif isnan(p(k))
        break;
    end
end

plot(xkp1(1),xkp1(2),'k.','markersize',30);
    
if nargout > 0
    pout = p;
end

end

function plot_fixed_point_contours(g,xm,ym)

global N

X = [reshape(xm,1,N*N); reshape(ym,1,N*N)];

F = g(X) - X;
Fn = reshape(sum(abs(F)),N,N);
contour(xm,xm,Fn,50,'k');

end