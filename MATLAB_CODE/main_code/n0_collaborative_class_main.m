clear all;
clc;
close all;

%% COLLABORATIVE CLASS
%   Matteo Mastromauro
%   Kevin Carrion
%   Alexandra Spertini
%   Alex Zagati

%% PARAMETERS

k = 0.8;
Jm = 4e-4;
Jl= 4e-4;
Bm = 0.015;
Bl = 0;
m = 0.3;
l = 0.3;
g = 9.81;
const = sqrt(2)*m*g*l/(2*Jl); 

%% 1. STATE SPACE REPRESENTATION OF THE NONLINEAR SYSTEM

% Simulink file "n1_nonlinear_system.slx"

%% 2. LINEARIZATION and EQUILIBRIUM
% linear tangent approximation at the point corresponding to θl =π/4

% 2.1 EQUILIBRIUM
x1_bar = 0;
x2_bar = pi/4;
x3_bar = 0;
x4_bar = m*g*l*cos(x2_bar)/k + x2_bar;
u_bar= m*g*l*cos(x2_bar);
x_bar = [x1_bar x2_bar x3_bar x4_bar]';

% SIMULINK MODEL INITILIZATION
dx_0 = [0;1;0;1]*0.01;           %integral initialization

% computation of the equilibrium position with MatLab
% 
% syms x1(t) x2(t) x3(t) x4(t) u(t) y(t) % variables
% syms  Bm Jl Jm m g l k %system parameters, comment if you want the numerical expression
% 
% x1_dot = -Bl/Jl*x1(t)-k/Jl*x2(t)-m*g*l/Jl*cos(x2(t))+k/Jl*x4(t);
% x2_dot = x1(t);
% x3_dot = k/Jm*x2(t)-Bm/Jm*x3(t)-k/Jm*x4(t)+u(t)/Jm;
% x4_dot = x3(t);
% 
% x1_dot = subs(x1_dot,x2(t),pi/4);
% x2_dot = subs(x2_dot,x2(t),pi/4);
% x3_dot = subs(x3_dot,x2(t),pi/4);
% x4_dot = subs(x4_dot,x2(t),pi/4);
% 
% eqn1 = 0 == x1_dot;
% eqn2 = 0 == x2_dot;
% eqn3 = 0 == x3_dot;
% eqn4 = 0 == x4_dot;
% 
% [x1_bar, x3_bar, x4_bar,u_bar] = solve([eqn1,eqn2,eqn3,eqn4],[x1(t),x3(t),x4(t),u(t)]) ;
% 
% x_bar = [x1_bar x2_bar x3_bar x4_bar] ;

% 2.2 LINEARIZATION

Alin= [ 0           -k/Jl + m*g*l*sin(x2_bar)/Jl             0           k/Jl;
        1                        0                           0             0;
        0                      k/Jm                       -Bm/Jm        -k/Jm;
        0                        0                           1              0];

Blin = [0; 0; 1/Jm; 0];

Clin= [0    1   0   0];            

Dlin = zeros(1);

% POLES OF THE SYSTEM

sys = ss(Alin, Blin, Clin, Dlin) ;
G = tf(sys) ;

pole(G); % computation of the poles of the linearised system 
zero(G); % computation of the zeros of the linearised system 

% 2.3 STABILIZING FEEDBACK CONTROL-LAW

% 2.3.1 POLE PLACEMENT

C0 = ctrb(Alin, Blin);          %controllability matrix
rank_C0 = rank(C0);             %rank of controllability matrix

% we decide to put the poles all in -50
rho = -50;                      % same rho for the two approaches to do some anaylsis 

p = rho + [-0, -1, -2, -3];     % poles choosen for pole placement
K = place(Alin, Blin, p);

p_pole_placement= eig(Alin-Blin*K);

% 2.3.2 POLE PLACEMENT WITH INTEGRATORS

A_ext= [Alin, zeros(4, 1);              %A extended matrix
       -Clin, zeros(1)];

B_ext = [Blin; zeros(1)];               %B extended matrix

C0_ext = ctrb(A_ext, B_ext);            %C extended matrix
rank_C0_ext = rank(C0_ext);

p_ext = rho + [0, -1, -2, -3, -4];      %poles of the extended system
K_ext = place(A_ext, B_ext, p_ext);     %pole placement gain

Kx = -K_ext(1, 1:4) ;                   %state feedback gain
Kv = -K_ext(1, 5) ;                     %additional integrator gain

p_pole_placement_ext = eig(A_ext-B_ext*K_ext); % poles of the extended closed loop system

% 2.3.3 OPTIMAL CONTROL : LQ 

% Tunable matrix intialization
Q = 0.01*eye(4);
R = 1*diag(1);

P_lq = are(Alin, Blin*R^(-1)*Blin', Q) ;
K_lq = R^(-1)*Blin'*P_lq;

p_LQ = eig(Alin-Blin*K_lq);

% 2.3.4 OPTIMAL CONTROL : LQ with the system enlargment

% now we want to enlarge the system in order to add an integrator 
A_tilde = [Alin zeros(4,1) ; -Clin zeros(1,1)] ;
B_tilde = [Blin; 0];

r =  1/0.001; 
q =  1;

%tunable matrixes 
Q_lq = q*eye(5) ; %1 
R_lq = r*eye(1) ; %1000

% check reachability and observability 
rank(ctrb(A_tilde, B_tilde));
C_q = sqrt(Q_lq) ;
rank(obsv(A_tilde, C_q));

P_lq_ext = are(A_tilde, B_tilde*R_lq^(-1)*B_tilde', Q_lq) ;
K_lq_ext = R_lq^(-1)*B_tilde'*P_lq_ext;

K_lqx = K_lq_ext(:,1:4) ;
K_lqeta = K_lq_ext(:,5) ;

p_LQ_ext = eig(Alin-Blin*K_lqx);

% 2.3.5 H2-CONTROL with the system enlargment

n = 5 ;             %n of states
p = 1 ;             %n of inputs
m = 1 ;             %number of outpuRts

q =       1;        %0.01
r =       1/0.008;        %10
qt =      0; 
rt =      0;

% definition of the matrixes 
A = [Alin zeros(4,1) ;
    -Clin   0 ];

B1 = [sqrt(qt)*eye(n) zeros(n,p)];
B2 = [Blin ; 0];

C1 = [sqrt(q)*eye(n) ; zeros(m,n)];
C2 = [Clin 0] ;

D11 = [zeros(n,n) zeros(n,p) ; zeros(m,n) zeros(m,p)];
D12 = [zeros(n,m) ; sqrt(r)*eye(m)] ;
D21 = [zeros(p,n)  sqrt(rt)*eye(p,p)] ;
D22 = [Dlin] ;

AA = A;
BB = [B1 B2] ;
CC = [C1 ; C2] ;
DD = [D11 D12 ; D21 D22] ;

P_H2 = are(A, B2*(D12'*D12)^(-1)*B2', C1'*C1);
K_H2 = (D12'*D12)^(-1)*B2'*P_H2;

K_H2_X = K_H2(:,1:4) ;
K_H2_V = K_H2(:, 5) ;
p_H2 = eig(Alin-Blin*K_H2_X);

%% 3. TEST THE LINEARIZED CONTROL-LAW

% Simulink file "n3_LTI_controller.slx"

%% 4 FEEDBACK LINEARIZATION

% 4.1 RELATIVE DEGREE OF NONLINEAR SYSTEM

% 4.1.1 Symbolic computation of the equilibrium
syms x1(t) x2(t) x3(t) x4(t) u(t) y(t) % variables
syms  Bm Jl Jm m g l k %system parameters, comment if you want the numerical expression

x1_dot = -Bl/Jl*x1(t)-k/Jl*x2(t)-m*g*l/Jl*cos(x2(t))+k/Jl*x4(t);
x2_dot = x1(t);
x3_dot = k/Jm*x2(t)-Bm/Jm*x3(t)-k/Jm*x4(t)+u(t)/Jm;
x4_dot = x3(t);

x1_dot = subs(x1_dot,x2(t),pi/4);
x2_dot = subs(x2_dot,x2(t),pi/4);
x3_dot = subs(x3_dot,x2(t),pi/4);
x4_dot = subs(x4_dot,x2(t),pi/4);

eqn1 = 0 == x1_dot;
eqn2 = 0 == x2_dot;
eqn3 = 0 == x3_dot;
eqn4 = 0 == x4_dot;

[x1_bar, x3_bar, x4_bar,u_bar] = solve([eqn1,eqn2,eqn3,eqn4],[x1(t),x3(t),x4(t),u(t)]);
 
% 4.1.2 computation of the relative degree "rd"
syms x1(t) x2(t) x3(t) x4(t) u(t) y(t) % variables
syms  Bm Jl Jm m g l k %system parameters, comment if you want the numerical expression
x1_dot = -Bl/Jl*x1(t)-k/Jl*x2(t)-m*g*l/Jl*cos(x2(t))+k/Jl*x4(t);
x2_dot = x1(t);
x3_dot = k/Jm*x2(t)-Bm/Jm*x3(t)-k/Jm*x4(t)+u(t)/Jm;
x4_dot = x3(t);
y(t) = pi/4-x2(t);

% checking the relative degree of the system
y_d = -x2_dot;
y_dd = -x1_dot;
y_ddd = diff(y_dd,t);
y_ddd  = (k*x2_dot)/Jl - (k*x4_dot)/Jl - (g*l*m*sin(x2(t))*x2_dot)/Jl;
y_dddd = diff(y_ddd,t);

%substituting diff(x,t) in order to include the exrpession in the simulink scheme

y_dddd = (k*x1_dot)/Jl - (k*x3_dot)/Jl - (g*l*m*sin(x2(t))*x1_dot)/Jl - (g*l*m*cos(x2(t))*x1(t)*x2_dot)/Jl;

%check the realtive degree of the system
y_d_sub = subs(y_d, [x1(t), cos(x2(t)), x2(t), x3(t), x4(t), u(t), k, Jl, Jm], [0, 0, 0, 0, 0, 1, 1, 1, 1]);
y_dd_sub = subs(y_dd, [x1(t), cos(x2(t)), x2(t), x3(t), x4(t), u(t), k, Jl, Jm], [0, 0, 0, 0, 0, 1, 1, 1, 1]);
y_ddd_sub = subs(y_ddd, [x1(t), cos(x2(t)), x2(t), x3(t), x4(t), u(t), k, Jl, Jm], [0, 0, 0, 0, 0, 1, 1, 1, 1]);
y_dddd_sub  = subs(y_dddd, [x1(t), cos(x2(t)), x2(t), x3(t), x4(t), u(t), k, Jl, Jm], [0, 0, 0, 0, 0, 1, 1, 1, 1]);

rd = zeros(1,1) ;

if y_d_sub ~= 0 
    rd =  1 ;
elseif y_dd_sub ~= 0 
        rd = 2 ;
elseif y_ddd_sub ~= 0 
        rd =  3 ;
elseif y_dddd_sub ~= 0 
        rd = 4;
end

% 4.2 STATE-FEEDBACK LINEARIZING CONTROL-LAW

% 4.2.1 Change of coordinates to bring equilibrium to [0, 0]
syms dx1(t) dx2(t) dx3(t) dx4(t)
dx1_dot = subs(x1_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx2_dot = subs(x2_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx3_dot = subs(x3_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx4_dot = subs(x4_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dy = subs(y(t), [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]); 

% 4.2.2 Diffeomorphism
x1_tilde = y;
x2_tilde = -dx1 ;
x3_tilde = dx1_dot ;
x4_tilde = (k*dx4_dot)/Jl - (k*dx2_dot)/Jl - (Bl*x1_dot)/Jl + (g*l*m*sin(pi/4 + dx2(t))*dx2_dot)/Jl;

% check if the diffeomorphism is well defined over the entire R^4: We need to compute the Jacobian and check
% if the Jacobian is full rank over the entire R^4
f = [x1_tilde ; x2_tilde ; x3_tilde ; x4_tilde] ;
J = jacobian(f, [dx1, dx2, dx3, dx4]);
rank(J);

% 4.2.3 Finding u(t)
syms x1 x2 x3 x4 u v_sim k Bm Jm Jl g l m 
v = (g*l*m*sin(x2)*((k*x2)/Jl - (k*x4)/Jl + (g*l*m*cos(x2))/Jl))/Jl - (k*(u/Jm - (Bm*x3)/Jm + (k*x2)/Jm - (k*x4)/Jm))/Jl - (k*((k*x2)/Jl - (k*x4)/Jl + (g*l*m*cos(x2))/Jl))/Jl - (g*l*m*cos(x2)*x1^2)/Jl;
eqn = v_sim == v;
U = solve(eqn, u);  % used in the block MatLab function "u computation"

% overwrite the symbolic parameter with the numerical values 
k = 0.8;
Jm = 4e-4;
Jl= 4e-4;
Bm = 0.015;
Bl = 0;
m = 0.3;
l = 0.3;
g = 9.81;
x1_bar = 0;
x2_bar = pi/4;
x3_bar = 0;
x4_bar = m*g*l*cos(x2_bar)/k + x2_bar;
u_bar= m*g*l*cos(x2_bar);
x_bar = [x1_bar x2_bar x3_bar x4_bar] ;

% 4.3 STABILIZING CONTROL-LAW on the feedback linearized system

% 4.3.1 FEEDBACK GAIN TUNING Kv_tilde
% linear system when applied the state-feedback linearizing control-law
A_tilde = [0      1       0       0;
           0      0       1       0;
           0      0       0       1;
           0      0       0       0;];

B_tilde = [0   0   0   1]';

C_tilde = [1   0   0   0];

p_tilde = -50 + [-0, -1, -2, -3]; % poles choosen for pole placement

Kv_tilde = place(A_tilde, B_tilde, p_tilde);
p_tilde_pole_placement = eig(A_tilde-B_tilde*Kv_tilde);

% 4.3.2 PI over the pole placement Kv_tilde
s = tf('s');
G_v =C_tilde*(s*eye(4)-(A_tilde-B_tilde*Kv_tilde))^-1*B_tilde;
wc = 3; % rad/s
R_v = 1/s;
L_v = R_v*G_v;
kp_v = 1/abs(freqresp(L_v,wc));
R_v = kp_v/s;
L_v = R_v*G_v;
phi_c = angle(freqresp(L_v,wc));
phi_m = (pi + phi_c)*180/pi;

%% 5. TEST THE STATE-FEEDBACK LINEARIZING CONTROL-LAW with LTI-controller

% Simulink file "n5_feedback_linearization_LTI_controller.slx"

%% 6. FEEDBACK LINEARIZATION with VARIABLE STRUCTURE CONTROLLER 

% 6.1 VARIABLE STRUCTURE CONTROLLER
% 6.1.1 linear system when applied the state-feedback linearizing control-law
%       the new state is called x_tilde
A_tilde = [0      1       0       0;
           0      0       1       0;
           0      0       0       1;
           0      0       0       0;];

B_tilde = [0   0   0   1]';

C_tilde = [1   0   0   0];

% 6.1.2 computation of beta_prime and alpha_prime coefficients
p = +10;                             % -p are the choosen poles when system slides on the sliding surface s(x) = 0 
beta_prime= [p^3 3*p^2 3*p 1];       % coefficients of the polynomial (s+p)^3
beta3= beta_prime(1);
b4= C_tilde(1);
gamma= beta3/b4;
alfa_prime= beta_prime*A_tilde;

% parameters q, r
q = 20;
r = 0.7;

% hysteresis parameters

BI = 0.02;
M = 1;

%% 7. FEEDBACK LINEARIZATION with VARIABLE STRUCTURE CONTROLLER other aspects analysis

% 7.1 ENLARGED SYSTEM WITH INTEGRATOR (for the high input signal switching issue)

% 7.1.1 enlargement of the system. The new state is called x_tilde2 = [x_tilde; u]
F = [A_tilde           B_tilde;
    zeros(1,4)  0];
G = [zeros(4,1); 1];
H = [C_tilde   0];

% 7.1.2 change of state variables
% matrix computation for the change of state variables x_canonical = T*x_tilde2 to get a new state 
% space reppresentation in the controllable canonical form: (F_tilde,G_tilde;H_tilde)
Mr = ctrb(F, G);

a_contr= poly(F);
a= flip(-a_contr);
% controllable canonical form we want to get
A_contr = [0     1       0       0       0;
           0     0       1       0       0;
           0     0       0       1       0;
           0     0       0       0       1;
           a(1:5)];

B_contr  = [0 0 0 0 1]';
Mr_contr = ctrb(A_contr, B_contr);

T= Mr_contr*Mr^-1;      % we get an identity matrix, i.e. T = eye(5), because 
                        % the system is already in the controllable canonical form

% space reppresentation in the controllable canonical form:
% the new state is x_canonical
F_contr = T*F*T^-1;
G_contr = T*G;
H_contr = H*T^-1;

% 7.1.3 computation of beta_prime_norm and alpha_prime_norm coefficients
p_contr = -10; 
beta_prime_contr = flip(poly([p_contr p_contr p_contr p_contr]));
beta3_contr= beta_prime_contr(1);
b4_contr= H_contr(1);
gamma_contr= beta3_contr/b4_contr;
alfa_prime_contr= beta_prime_contr*F_contr;


% 7.3 LOAD DISTURBANCES ON VSC

W = 2;      % amplitude of the disturbance w(t)

% Simulink file "n7_feedback_linearization_VSC_controller.slx"


