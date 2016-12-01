clear all 
clc

syms x1 x2 real
fx = x1+x2;
hx = ([x1-1;-x1-1;x2-2;-x2]);
vp = [x1,x2];

t = 0.5;
kmax = 100;
tbt = 0.1;
alpha = 0.2;
beta = 0.8;
eps = 1e-6;
mi = 0.5;

% function with feasible starting point calculation for convex constraints 
xopt = IP_GDM_convex_feasible_start(fx, hx, kmax, t, tbt, alpha, beta,eps,mi,vp)
