function pout = cj(A,b,tol,kmax,fignum)

figure(fignum);
clf;

title('Conjugate Gradient','fontsize',18);
set(gca,'fontsize',16);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);

hold on;


plot_single_contour = true;

F = @(x) x'*A*x - 2*b'*x;

plot_contours(A,b);
fprintf('Choose starting guess\n');
xk = ginput(1);
xk = xk(:);
if (plot_single_contour)
    Fx = F(xk);
    plot_contours(A,b,Fx);
    pause(0.5);
end

rk = b - A*xk;
dk = rk;

clear p;
for k = 1:kmax
    
    Adk = A*dk;  % This is the only matrix vector multiply that we need.
    % Adk = matvec(dk);  % Matrix-free method!
    
    % Update x. 
    a1 = dot(rk,rk)/dot(dk,Adk);
    xkp1 = xk + a1*dk;     % blas operation : gaxpy : genearlized ax + y
    % rkp1 = b - A*xkp1;
    rkp1 = rk - a1*Adk;    % gaxpy
    
    % Get A-conjugate directions
    a2 = dot(rkp1,rkp1)/dot(rk,rk);   % dot product x'*y  AllReduce
    dkp1 = rkp1 + a2*dk;              % dot product
    
    % Compute the norm of the residual.
    p(k) = norm(rkp1,1);
    fprintf('%4d %16.6e\n',k,p(k));
        
    dk = dkp1;
    rk = rkp1;
    plot([xk(1) xkp1(1)],[xk(2) xkp1(2)],'r','linewidth',2);
    hold on;
    plot([xk(1) xkp1(1)],[xk(2) xkp1(2)],'k.','markersize',30);
    if (plot_single_contour)
        Fx = F(xk);
        plot_contours(A,b,Fx);
        pause(0.5);
    end
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
