function [pout,eout] = steepest_descent(A,b,tol,kmax,fignum)

plot_single_contour = true;

F = @(x) 0.5*x'*A*x - b'*x;

figure(fignum);
clf;

title('Steepest-Descent','fontsize',18);
set(gca,'fontsize',16);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
hold on;

hout = plot_contours(A,b);
set(hout,'visible','off');
fprintf('Choose starting guess\n');
xk = ginput(1);
xk = xk(:);


if (plot_single_contour)
    Fx = F(xk);
    hout = plot_contours(A,b,Fx);
    set(hout,'color','k');
    % input('Hit enter to see next iterate :');
    pause(0.5);
end

%{
% If the error is an eigenvector, SD will converges in 1 step.
[evec,eval] = eig(A);

xbar = A\b;
xk = xbar + evec(:,1);  
%}

rk = b - A*xk;
pk = rk;

resid = zeros(kmax,1);
err = zeros(kmax,1);

resid(1) = norm(rk,1);
for k = 1:kmax
    
    Apk = A*pk;
    
    % Update x. 
    a1 = dot(rk,rk)/dot(pk,Apk);
    xkp1 = xk + a1*pk;
    rkp1 = rk - a1*Apk;
    pkp1 = rkp1;
    
    % Compute the norm of the residual.
    resid(k) = norm(rkp1,1);
    err(k) = norm(xkp1-xk,1);
    fprintf('%4d %16.6e\n',k,resid(k));    
        
    rk = rkp1;
    pk = pkp1;
    
    plot([xk(1) xkp1(1)],[xk(2) xkp1(2)],'r','linewidth',2);
    hold on;
    plot([xk(1) xkp1(1)],[xk(2) xkp1(2)],'k.','markersize',30);
    if (plot_single_contour)
        Fx = F(xk);
        hout = plot_contours(A,b,Fx);
        set(hout,'color','k');
        % input('Hit enter to see next iterate :');
        pause(0.5);
    end
    xk = xkp1;
    if (resid(k) < tol)
        break;
    elseif isnan(resid(k))
        break;
    end
end
plot(xkp1(1),xkp1(2),'k.','markersize',30);

if (nargout > 0)
    pout = resid;
    eout = err;
end

end


