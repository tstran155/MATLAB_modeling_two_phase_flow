close all;
clear all;
clc;

global swc sor no nw kro0 krw0 mu_o mu_w;
global c_w c_o b_w_0 b_o_0 po_0 q_o_0 q_w_0;
global A L Nx Nt max_time max_iter;
global phi ki xi ti swi poi dx dt plots;
global boi bwi krwi kroi sw_xt po_xt;
global max_diff dp info lambda_w lambda_o;
global q_c pvs q_r bmin bmax loop_plot;
global qrmin qrmax pmin pmax;

swc=0.2;
sor=0.1;
no=2;
nw=3;
kro0=0.8;
krw0=0.7;
mu_o=1e-3;
mu_w=1e-3;

c_w=1e-10;
c_o=1e-8;

po_0=1e7;
po_sc=1e5;
b_w_0=exp(c_w*(po_sc-po_0));
%allowing for disolved gas
b_o_0=1.3;
q_o_0=0;
q_w_0=12e-6/60;

A=40e-4;
L=0.2;

max_diff=1e-6;
dp=1e-6;
phi=0.2;
pvs=4;
max_time=pvs*(1-swc-sor)*A*L*phi/(q_w_0*b_w_0+q_o_0*b_o_0);
max_iter=1e4;
info=0;
k_0=1e-12;
plots=1;

bmin=0.5;
bmax=1.5;
qrmin=0;
qrmax=20e-6;
pmin=9.995e6;
pmax=1.02e7;

Nx=20;
Nt=100;
if(Nx<3)
    disp('Increase Nx');
    keyboard;
end
if(Nt<2)
    disp('Increase Nt');
    keyboard;
end

dx=L/(Nx-1);
dt=max_time/(Nt-1);

ti=linspace(0,max_time,Nt);
xi=linspace(0,L,Nx);
poi(1:Nx,1)=po_0;
swi(1,1)=swc+(1-swc-sor)*(q_w_0*b_w_0)/(q_w_0*b_w_0+q_o_0*b_o_0);
swi(2:Nx,1)=swc;
ki(1:Nx,1)=k_0;

po_xt=zeros(Nx,Nt);
sw_xt=zeros(Nx,Nt);
q_c=zeros(Nt,2);
q_r=zeros(Nt,2);
loop_plot=1;

tic
two_phase_flow_simulator();
if(abs(1-sor-max(max(sw_xt)))>max_diff || abs(swc-min(min(sw_xt)))>max_diff)
    disp('try using a smaller time step for stability of saturation');
end
toc