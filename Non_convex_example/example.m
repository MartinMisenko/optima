%% Example
% This script shows calculating of local optimum of non-convex function 
% with two non-convex constraints and box constraints using GDM method 
% with interior point barrier method and feasible start function.

% There is no commnon setup of parameters for GDM and Interior point. 

clear all
close all

syms x y real
fx = 3*(1-x)^2*exp(-(x^2) - (y+1)^2) ... 
   - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) ... 
   - 1/3*exp(-(x+1)^2 - y^2);

hx = [-(x+0.3)^2 - (y+0.3)^2 + 0.25;-(x+0.2)^2 - (y-1.6)^2 + 0.36;...
      y-3;-y-3;x-1.5;-x-3];

vp = [x,y];

%% PHASE 1

x0 = feasibleX(hx,vp)

t = 15;             % parameters for GDM and Interior point
kmax = 30;
tbt = 0.001;
alpha = 0.2;
beta = 0.5;
eps = 1e-4;
mi = 0.05;
%% GDM

[xx,xy] = InteriorPointGDM(fx, hx, x0, kmax, t, tbt, alpha, beta,eps,mi,vp) 


%% fmincon

% fun = @(x) (3*(1-x(1))^2*exp(-(x(1)^2) - (x(2)+1)^2) ... 
%    - 10*(x(1)/5 - x(1)^3 - x(2)^5)*exp(-x(1)^2-x(2)^2) ... 
%    - 1/3*exp(-(x(1)+1)^2 - x(2)^2));
% x0 = [1;-2];
% A = [];
% b = [];
% 
% x = fmincon(fun,x0,A,b)

%% plotting
[x,y] = meshgrid(-4:0.1:3,-4:0.1:4);
z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2);

figure
meshc(x,y,z)
xlabel('x')
ylabel('y')
title('Non-convex function')

figure
contour(x,y,z),colorbar,hold on 
xlabel('x')
ylabel('y')
title('Local optimum')
plot(xy(1,:),xy(2,:),'r*'),hold on
plot(x0(1),x0(2),'r*'),hold on

x = 0.3;                % non convex constraint
y = 0.3;
d = 0.6;
ang=0:0.01:2*pi; 
xp=d*cos(ang);
yp=d*sin(ang);
plot(x+xp,y+yp,'r','LineWidth',2),axis equal,grid on,hold on

x = 0.2;                % non convex constraint
y = -1.6;
d = 0.5;
ang=0:0.01:2*pi; 
xp=d*cos(ang);
yp=d*sin(ang);
plot(x+xp,y+yp,'r','LineWidth',2), hold on

k = [-3:0.1:3];         % box constraints
l = ones(length(k))*1.5;
plot(l,k,'r--','LineWidth',2),hold on

k = [-3:0.1:3];
l = ones(length(k))*(-3);
plot(l,k,'r--','LineWidth',2),hold on

k = [-3:0.1:1.5];
l = ones(length(k))*(-3);
plot(k,l,'r--','LineWidth',2),hold on

k = [-3:0.1:1.5];
l = ones(length(k))*3;
plot(k,l,'r--','LineWidth',2),hold on



n = [2:1:length(xy)+1];         % convergence of x,y
figure
plot(1,x0(1),'bx'),hold on
plot(1,x0(2),'rx'),hold on
plot(n,xy(1,:),'bx',n,xy(2,:),'rx')
xlabel('Iterations')
ylabel('x, y')
title('Convergence of x, y')
grid on
legend('x','y')


