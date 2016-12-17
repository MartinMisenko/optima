%% Example for Sequential Quadratic Programming procedure
clc; clear; close all;


fun = @(x, y) [ (x.^4) - 3.*(x.^3) + 1.1.*(x.^2) - 0.25.*x  - 0.7.*(x.^2).*(y.^2) + 1.15.*(y.^4) + 2.*x.*y - 3.*(y.^3) - 14.*(y.^2) + 2.3.*y - 14 ];

a = [-6 : 0.1 : 6];
b = [-6 : 0.1 : 6];
[X, Y] = meshgrid(a,b);
Z = fun(X, Y);

%% EXAMPLE WITH YALMIP

clear x y
sdpvar x y
fx = fun(x, y);
cst = [ (x-2)^2 + (y-2)^2 <= 16 ; 
        y - 0.8*x <= 4;
        y - (x.^2) + 6*x >= 7.3;
        x >= 0; x <= 6;
        y >= -2; x <= 7];
optimize(cst, fx)

x_opt = value(x);
y_opt = value(y);
z_opt = value(fx);

% plot yalmip result
figure
mesh(X,Y,Z);
hold on
plot3(x_opt, y_opt, z_opt, 'kx', 'LineWidth', 3, 'MarkerSize',15)
tit = sprintf('Result with YALMIP\n x* = %.2f, y* = %.2f, f(x*,y*) = %.2f', x_opt, y_opt, z_opt);
title(tit);
xlabel('x')
ylabel('y')
zlabel('f(x,y)')

%---- plot constraints into 3D graph ----
% --- h1 ---
x1 = [-2 : 0.1 :6 ];
y1 = sqrt( 16 - (x1-2).^2) + 2;
y11 = -sqrt( 16 - (x1-2).^2) + 2;
zz1 = fun(x1, y1);
zz11 = fun(x1, y11);
plot3(x1, y1,zz1, 'g', 'LineWidth',3);
plot3(x1, y11, zz11, 'g', 'LineWidth',3);
% --- h2 ---
x2 = [-2 : 0.1 : 2.5 ];
y2 = 0.8.*x2 + 4;
z2 = fun(x2, y2);
plot3(x2,y2,z2, 'g', 'LineWidth',3);
% --- h3 ---
x3 = [0.3 : 0.1 : 5.5 ];
y3 = (x3.^2) - 6.*x3 + 7.3;
z3 = fun(x3, y3);
plot3(x3,y3,z3, 'g', 'LineWidth',3);
% ------------------------------------


%% Example with SQP
clear x y
syms x y real

fx = x^4 - 3*x^3 + 1.1*x^2 - 0.25*x  - 0.7*x^2*y^2 + 1.15*y^4 + 2*x*y - 3*y^3 - 14*y^2 + 2.3*y - 14;
hx = [ (x-2)^2 + (y-2)^2 - 16 ; 
        y - 0.8*x - 4;
        -y + (x.^2) - 6*x + 7.3;
        -x;     x - 6;
        -y + 2; y - 7; ];

eps = 1e-6;
kmax = 5000;

% ========== CHANGE x0 ============
%x0 = feasibleX(hx,[x, y]);
x0 = [1; 4];
% =================================

% ---- call function -----
tic;
convergence = SQP(fx,hx,x0,kmax,eps);
time = toc;
n = size(convergence,2);    % number of iterations
msg = sprintf('Successfully solved with SQP\nElapsed time:  %.3f\nNNumber of iterations:  %d',time, n);
disp(msg)

% --- optimal values ------
x_opt = convergence(:, end);
fx_opt = fun(x_opt(1), x_opt(2));

%% plot results
clear x y

figure
mesh(X,Y,Z);
tit = strcat('',sprintf('Result with SQP from x_0 = [%.1f, %.1f]\n x* = %.2f, y* = %.2f, f(x*,y*) = %.2f', x0(1), x0(2), x_opt(1), x_opt(2), fx_opt));
title(tit);
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
hold on
plot3(x_opt(1), x_opt(2), fx_opt, 'kx', 'LineWidth',3, 'MarkerSize',15);
plot3(x0(1), x0(2), fun(x0(1), x0(2)), 'rx', 'LineWidth',3, 'MarkerSize',15);

%---- plot constraints into 3D graph ----
% --- h1 ---
plot3(x1, y1,zz1, 'g', 'LineWidth',3);
plot3(x1, y11, zz11, 'g', 'LineWidth',3);
% --- h2 ---
plot3(x2,y2,z2, 'g', 'LineWidth',3);
% --- h3 ---
plot3(x3,y3,z3, 'g', 'LineWidth',3);
% ----------


%% plot contour to show convergence

figure
contour(a,b,Z)
hold on
%plot(x_opt(1), x_opt(2),'kx', 'LineWidth', 3)
plot(convergence(1, 1:(end-1)), convergence(2, 1:(end-1)),'kx', 'LineWidth', 3)
plot(x_opt(1), x_opt(2),'rx', 'LineWidth', 3)
plot(x1, y1, 'g');
plot(x1, y11, 'g');
plot(x2,y2, 'g');
plot(x3,y3, 'g');
title('Contour plot - Optimal solution with individual iterations');
leg = legend('f(x,y) contour plot','iterations of f*(x,y)','optimal f(x,y)', 'constraints');
set(leg, 'Location','northwest');
xlabel('x');
ylabel('y');

%% plot only convergence during iterations


t = [1 : 1 : n];    % n is number of iterations
figure
% --- variable value convergence ----
subplot(2,1,1)
hold on
plot(t,convergence(1,:),'r--o', 'LineWidth', 2,'MarkerSize', 12);
plot(t,convergence(2,:),'b--o', 'LineWidth', 2,'MarkerSize', 12 );
title('Convergence of optimization variables x, y during the SQP procedure');
xlabel('Iteration no.');
ylabel('Variable value');
legend('value of x', 'value of y');

% --- corresponding function value ---
fconv = fun(convergence(1,:), convergence(2,:));
subplot(2,1,2)
plot(t,fconv,'g--o', 'LineWidth', 2,'MarkerSize', 12 );
title('Corresponding function value');
xlabel('Iteration no.');
ylabel('f(x,y)');
