clear all 
clc

syms x1 x2 real
fx = x1+x2;
hx = ([x1-1;-x1-1;x2-2;-x2]);
vp = [x1,x2];

t = 0.5;
x0 = [0;1];
kmax = 100;
tbt = 0.1;
alpha = 0.2;
beta = 0.8;
eps = 1e-6;
mi = 0.5;

xopt = InteriorPointGDM(fx, hx, x0, kmax, t, tbt, alpha, beta,eps,mi,vp)
