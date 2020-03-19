function plot_contours(A,b,s)

xbar = A\b;

% [evec,~] = eig(A);

xlim = [xbar(1)-4 xbar(1)+4];
ylim = [xbar(2)-4 xbar(2)+4];
xc = linspace(xlim(1),xlim(2),500);
yc = linspace(ylim(1),ylim(2),500);
[xm,ym] = meshgrid(xc,yc);

t1 = xm.*(A(1,1)*xm + A(1,2).*ym) + ym.*(A(2,1)*xm + A(2,2)*ym);
t2 = b(1)*xm + b(2)*ym;
f = t1 - 2*t2;

if (nargin == 3)
    cv = [s,s];
    lw = 2;
else
    cv = 20;
    lw = 1;
end
contour(xm,ym,f,cv,'linewidth',lw);
hold on;

plot(xbar(1),xbar(2),'r.','markersize',30);

drawnow

daspect([1,1,1]);

%%{
[evec,~] = eig(A);
t = linspace(-2,2,20);
p1 = xbar + 5*evec(:,1);
p2 = xbar - 5*evec(:,1);
plot([p1(1) p2(1)],[p1(2) p2(2)],'k--');
p1 = xbar + 5*evec(:,2);
p2 = xbar - 5*evec(:,2);
plot([p1(1) p2(1)],[p1(2) p2(2)],'k--');
axis([xlim ylim]);
%}

