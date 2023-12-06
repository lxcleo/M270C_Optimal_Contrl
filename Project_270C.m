clear all;
% State equations
syms x1 x2 p1 p2 u;
Dx1 = x2;
Dx2 = -x2 + u;
% Cost function inside the integral
syms g;
g = 0.5*u^2;
% Hamiltonian
syms p1 p2 H;
H = g + p1*Dx1 + p2*Dx2;
% Costate equations
Dp1 = -diff(H,x1);
Dp2 = -diff(H,x2);
% solve for control u
du = diff(H,u);
sol_u = solve(du,u);
% Substitute u to state equations
Dx2 = subs(Dx2,u,sol_u);
%%
syms x y z psi theta phi alpha beta ...
    x_d y_d z_d psi_d theta_d phi_d alpha_d beta_d ...
    p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 ...
    p12 p13 p14 p15 p16 u1 u2 u3 u4 real;
q = [x y z psi theta phi alpha beta]';
zeta = [x y z]';
neta = [psi theta phi]';
miu = [alpha beta]';
q_d = [x_d y_d z_d psi_d theta_d phi_d  alpha_d beta_d]';
zeta_d = [x_d y_d z_d]';
neta_d = [psi_d theta_d phi_d]';
miu_d = [alpha_d beta_d]';
ml = 10; % Mass of the load
M_quad = 4; % Mass of the quadrotor
l = 0.8; % Length of the quadrotor 
m11 = M_quad + ml;
m22 = m11;
m33 = m11;
m17 = ml * l * cos(alpha) * sin(beta);
m71 = m17;
m18 = -ml * l * sin(alpha) * sin(beta);
m81 = m18;
m27 = ml * l * cos(alpha) * sin(beta);
m72 = m27;
m28 = ml * l * sin(alpha) * cos(beta);
m82 = m28;
m37 = ml * l * sin(alpha);
m73 = m37;
I_psi = 1;
I_theta = 0.5;
I_phi = 0.5;

m44 = I_psi * sin(theta)^2 + cos(theta)^2 * (I_theta * sin(theta)^2 ...
    + I_phi * cos(phi)^2);
m45 = (I_theta - I_phi)* (cos(theta) * sin(phi) * cos(phi));
m54 = m45;
m55 = I_theta * cos(phi)^2 + I_phi * sin(phi)^2;
m77 = ml * l^2;
m88 = ml * l^2 * sin(alpha)^2;
M = [m11 0 0 0 0 0 m17 m18; ...
    0 m22 0 0 0 0 m27 m28; ...
    0 0 m33 0 0 0 m37 0; ...
    0 0 0 m44 m45 -I_psi * sin(alpha) 0 0;...
    0 0 0 m54 m55 0 0 0;...
    0 0 0 -I_psi * sin(theta) 0 I_psi 0 0;...
    m71 m72 m73 0 0 0 m77 0;...
    m81 m82 0 0 0 0 0 m88];
c17 = -ml * l * (cos(alpha) * sin(beta) * beta_d + sin(alpha)...
    * cos(beta) * alpha_d);
c18 = -ml * l *(cos(alpha) * sin(beta) * alpha_d + sin(alpha)...
    * cos(beta) * beta_d);
c27 = ml * l * (cos(alpha) * cos(beta) * beta_d - sin(alpha) ...
     * sin(beta) * alpha_d);
c28 = ml * l * (cos(alpha) * cos(beta) * alpha_d - sin(alpha) ...
     * sin(beta) * beta_d);
c44 = I_psi * theta_d * sin(theta) * cos(theta) - (I_theta + ...
    I_phi) * (theta_d * sin(theta) * cos(theta) * sin(theta)^2);
c45 = I_psi * psi_d * sin(theta) * cos(theta) - (I_theta + I_phi)...
    *(psi_d^2 * sin(theta) * cos(theta) * cos(theta)^2 - phi_d * ...
    cos(theta) * cos(phi)^2);
c46 = -(I_psi * theta_d * cos(theta));
c54 = psi_d * sin(theta) * cos(theta) * (-I_psi + I_theta * sin(theta)^2 ...
    + I_phi * cos(phi)^2);
c55 = 0;
c56 = I_psi * psi_d * cos(theta);
c64 = 0;
c65 = -I_psi * psi_d * cos(theta);
C = [0 0 0 0 0 0 c17 c18;...
    0 0 0 0 0 0 c27 c28; ...
    0 0 0 0 0 0 ml * l * cos(alpha) * alpha_d 0;...
    0 0 0 c44 c46 c46 0 0;...
    0 0 0 c54 c55 c56 0 0;...
    0 0 0 c64 c65 0 0 0;...
    0 0 0 0 0 0 0 -ml * l^2 * sin(alpha) * cos(beta) * beta_d;...
    0 0 0 0 0 0 ml * l^2 * sin(alpha) * cos(beta) * beta_d ml ...
    * l^2 * sin(alpha) * cos(beta) * alpha_d];
g = 9.8;
G = [0 0 (M_quad + ml) * g 0 0 0 ml * l * g * sin(alpha) 0]';
b11 = sin(alpha) * sin(psi) + cos(psi)* cos(phi) * sin(theta);
b21 = sin(psi)* cos(phi) * sin(theta) - cos(psi) * sin(phi);
b = [b11 0 0 0;...
    b21 0 0 0;...
    cos(theta) * cos(phi) 0 0 0;...
    0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 0 0 0;...
    0 0 0 0];
dydt(1:8) = q_d;
u = [u1 u2 u3 u4]';
U = b * u;
dydt(9:16) = -inv(M) * (C * q_d + G) + inv(M) * U;
dydt = simplify(dydt');
P = [p1 p2 p3 p4 p5 p6 p7 p8 ...
    p9 p10 p11 p12 p13 p14 p15 p16]';
%%
syms J real;
Q = eye(8);
J = 1/2 * q' * Q * q + 1/2 * q_d' * Q * q_d + 1/2 * U' * Q * U;
syms H real;
H = J + P' * dydt;
du1 = diff(H,u1);
sol_u1 = solve(du1,u1);
%%
function dydt = pd(t,xx)
dydt = zeros(16,1);
x = xx(1);
y = xx(2);
z = xx(3);
psi = xx(4);
theta = xx(5);
phi = xx(6);
alpha = xx(7);
beta = xx(8);
x_d = xx(9);
y_d = xx(10);
z_d = xx(11);
psi_d = xx(12);
theta_d = xx(13);
phi_d = xx(14);
alpha_d = xx(15);
beta_d = xx(16);
q = [x y z psi theta phi  alpha beta]';
zeta = [x y z]';
neta = [psi theta phi]';
miu = [alpha beta]';
q_d = [x_d y_d z_d psi_d theta_d phi_d  alpha_d beta_d]';
zeta_d = [x_d y_d z_d]';
neta_d = [psi_d theta_d phi_d]';
miu_d = [alpha_d beta_d]';
ml = 10; % Mass of the load
M_quad = 4; % Mass of the quadrotor
l = 0.8; % Length of the quadrotor 
m11 = M_quad + ml;
m22 = m11;
m33 = m11;
m17 = ml * l * cos(alpha) * sin(beta);
m71 = m17;
m18 = -ml * l * sin(alpha) * sin(beta);
m81 = m18;
m27 = ml * l * cos(alpha) * sin(beta);
m72 = m27;
m28 = ml * l * sin(alpha) * cos(beta);
m82 = m28;
m37 = ml * l * sin(alpha);
m73 = m37;
I_psi = 1;
I_theta = 0.5;
I_phi = 0.5;

m44 = I_psi * sin(theta)^2 + cos(theta)^2 * (I_theta * sin(theta)^2 ...
    + I_phi * cos(phi)^2);
m45 = (I_theta - I_phi)* (cos(theta) * sin(phi) * cos(phi));
m54 = m45;
m77 = ml * l^2;
m88 = ml * l^2 * sin(alpha)^2;
M = [m11 0 0 0 0 0 m17 m18; ...
    0 m22 0 0 0 0 m27 m28; ...
    0 0 m33 0 0 0 m37 0; ...
    0 0 0 m44 m45 -I_psi * sin(alpha) 0 0;...
    0 0 0 m54 m55 0 0 0;...
    0 0 0 -I_psi * sin(theta) 0 I_psi 0 0;...
    m71 m72 m73 0 0 0 m77 0;...
    m81 m82 0 0 0 0 0 m88];
c17 = -ml * l * (cos(alpha) * sin(beta) * beta_d + sin(alpha)...
    * cos(beta) * alpha_d);
c18 = -ml * l *(cos(alpha) * sin(beta) * alpha_d + sin(alpha)...
    * cos(beta) * beta_d);
c27 = ml * l * (cos(alpha) * cos(beta) * beta_d - sin(alpha) ...
     * sin(beta) * alpha_d);
c28 = ml * l * (cos(alpha) * cos(beta) * alpha_d - sin(alpha) ...
     * sin(beta) * beta_d);
c44 = I_psi * theta_d * sin(theta) * cos(theta) - (I_theta + ...
    I_phi) * (theta_d * sin(theta) * cos(theta) * sin(theta)^2);
c45 = I_psi * psi_d * sin(theta) * cos(theta) - (I_theta + I_phi)...
    *(psi_d^2 * sin(theta) * cos(theta) * cos(theta)^2 - phi_d * ...
    cos(theta) * cos(phi)^2);
c46 = -(I_psi * theta_d * cos(theta));
c54 = psi_d * sin(theta) * cos(theta) * (-I_psi + I_theta * sin(theta)^2 ...
    + I_phi * cos(phi)^2);
c55 = 0;
c56 = I_psi * psi_d * cos(theta);
c64 = 0;
c65 = -I_psi * psi_d * cos(theta);
C = [0 0 0 0 0 0 c17 c18;...
    0 0 0 0 0 0 c27 c28; ...
    0 0 0 0 0 0 ml * l * cos(alpha) * alpha_d 0;...
    0 0 0 c44 c46 c46 0 0;...
    0 0 0 c54 c55 c56 0 0;...
    0 0 0 c64 c65 0 0 0;...
    0 0 0 0 0 0 0 -ml * l^2 * sin(alpha) * cos(beta) * beta_d;...
    0 0 0 0 0 0 ml * l^2 * sin(alpha) * cos(beta) * beta_d ml ...
    * l^2 * sin(alpha) * cos(beta) * alpha_d];
g = 9.8;
G = [0 0 (M_quad + ml) * g 0 0 0 ml * l * g * sin(alpha) 0]';
b11 = sin(alpha) * sin(psi) + cos(psi)* cos(phi) * sin(theta);
b21 = sin(psi)* cos(phi) * sin(theta) - cos(psi) * sin(phi);
b = [b11 0 0 0;...
    b21 0 0 0;...
    cos(theta) * cos(phi) 0 0 0;...
    0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 0 0 0;...
    0 0 0 0];
dydt(1) = xx(9);
dydt(2) = xx(10);
dydt(3) = xx(11);
dydt(4) = xx(12);
dydt(5) = xx(13);
dydt(6) = xx(14);
dydt(7) = xx(15);
dydt(8) = xx(16);
dydt(9:16) = -inv(M) * (C * xx(9:16) + G);
end