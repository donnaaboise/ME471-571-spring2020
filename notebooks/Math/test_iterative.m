function test_iterative()

close all

N = 2;
B = -1 + 2*rand(N,N);
b = -1 + 2*rand(N,1);

tol = 1e-6;
kmax = 1000;

figure(2);
clf;

% ------------------------------ Jacobi -----------------------------------
figure(1);
mstr = 'Jacobi';
A = get_matrix(B,'SPD');
disp(A)
r = fixed_point(A,b,mstr,tol,kmax);

figure(2);
p(1) = semilogy(r,'b.-','markersize',20);
hold on;
lstr{1} = mstr;

% ---------------------------- Gauss-Seidel -------------------------------
fprintf('\n');
input('Hit enter to see Gauss-Seidel : ');

mstr = 'Gauss-Seidel';
A = get_matrix(B,'SPD');
r = fixed_point(A,b,mstr,tol,kmax);

figure(2);
p(2) = semilogy(r,'r.-','markersize',20);
lstr{2} = mstr;

% -------------------------- Steepest-Descent -----------------------------
fprintf('\n');
A = get_matrix(B,'SPD');
input('Hit enter to see Steepest-descent : ');

mstr = 'Steepest-descent';
r = steepest_descent(A,b,tol,kmax);
figure(2);
p(3) = semilogy(r,'m.-','markersize',20);
lstr{3} = mstr;

% -------------------------- Conjugate Gradient ---------------------------
fprintf('\n');
input('Hit enter to see Conjugate gradient : ');

mstr = 'Conjugate-gradient';
A = get_matrix(B,'SPD');
r = cj(A,b,tol,kmax);
figure(2);
p(4) = semilogy(r,'c.-','markersize',20);
lstr{4} = mstr;


% ---------------------------- Finish plotting ----------------------------
figure(2)

plot(xlim,[tol,tol],'k--');

yl = ylim;
set(gca,'ylim',[tol/10,yl(2)]);

lh = legend(p,lstr);
set(lh,'fontsize',16,'AutoUpdate','off');

title('Convergence behavior','fontsize',18);
xlabel('Number of iterations');
ylabel('Error');
set(gca,'fontsize',16);
set(gca,'yscale','log');

shg;

end

function A = get_matrix(B,method)


switch method
    case 'DD'
        % Diagonally dominant
        beta = 0.01;   % -1 < beta < 1 : strictly diagonally dominant
        d = diag(B);
        D = diag(d);
        A1 = B - D;          % Subtract diagonal from B
        S = sum(abs(A1'))';  % compute absolution row  sum of off diagonals
        v = beta*d./S;       % Compute scaling factors for off-diagonal entries
        A = A1.*v + D;       
    case 'SPD'
        A = B'*B;
        
    otherwise
        A = B;

end
    
end