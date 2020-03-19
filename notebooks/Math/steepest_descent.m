function pout = steepest_descent(A,b,tol,kmax,fignum)

plot_single_contour = true;

F = @(x) x'*A*x - 2*b'*x;

figure(fignum);
clf;

title('Steepest-Descent','fontsize',18);
set(gca,'fontsize',16);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);


plot_contours(A,b);
fprintf('Choose starting guess\n');
xk = ginput(1);
xk = xk(:);
if (plot_single_contour)
    Fx = F(xk);
    plot_contours(A,b,Fx);
    % input('Hit enter to see next iterate :');
    pause(0.5);
end

rk = b - A*xk;
dk = rk;

for k = 1:kmax
    
    Adk = A*dk;
    
    % Update x. 
    a1 = dot(rk,rk)/dot(dk,Adk);
    xkp1 = xk + a1*dk;
    rkp1 = b - A*xkp1;
    dkp1 = rkp1;
    
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
        % input('Hit enter to see next iterate :');
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

if (nargout > 0)
    pout = p;
end

end


