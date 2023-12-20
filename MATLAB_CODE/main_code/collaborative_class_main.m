clear all;
clc;
close all;

%% COLLABORATIVE CLASS
%   Matteo Mastromauro
%   Kevin Ortega
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
const = sqrt(2)*m*g*l/(2*Jl); %% needed for simulink simulation

%% EQUILIBRIUM
x1_bar = 0;
x2_bar = pi/4;
x3_bar = 0;
x4_bar = m*g*l*cos(x2_bar)/k + x2_bar;
u_bar= m*g*l*cos(x2_bar);
x_bar = [x1_bar x2_bar x3_bar x4_bar] ;

% computation of the equilibrium position

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
% 
% eqn1 = 0 == x1_dot;
% eqn2 = 0 == x2_dot;
% eqn3 = 0 == x3_dot;
% eqn4 = 0 == x4_dot;
% 
% [x1_bar, x3_bar, x4_bar,u_bar] = solve([eqn1,eqn2,eqn3,eqn4],[x1(t),x3(t),x4(t),u(t)]) ;
% 
% x_bar = [x1_bar x2_bar x3_bar x4_bar] ;

%% LINEARIZATION

Alin= [0    -k/Jl + m*g*l*sin(x2_bar)/Jl                     0           k/Jl;
        1                        0                           0             0;
        0                      k/Jm                       -Bm/Jm      -k/Jm;
        0                        0                           1              0];

Blin = [0; 0; 1/Jm; 0];

Clin= [0    1   0   0];             %applichiamo a x2(theta l) perchè abbiamo solo 1 input(2 output)

Dlin = zeros(1);


%% POLES OF THE SYSTEM

sys = ss(Alin, Blin, Clin, Dlin) ;
G = tf(sys) ;

pole(G); % computation of the poles of the linearised system 
zero(G); % computation of the zeros of the linearised system 


%% POLE PLACEMENT (FOR STABILIZATION)

C0 = ctrb(Alin, Blin);          %controllability matrix
rank_C0 = rank(C0);  %rank of controllability matrix

%we decide to put the poles all in -50
rho = -50; % same rho for the two approaches to do some anaylsis 

p = rho + [-0, -0.1, -0.2, -0.3]; % poles choosen for pole placement
K = place(Alin, Blin, p);

p_closeloop= eig(Alin-Blin*K);

%% POLE PLACEMENT WITH INTEGRATORS

A_ext= [Alin, zeros(4, 1);            %A extended matrix
       -Clin, zeros(1)];

B_ext = [Blin; zeros(1)];               %B extended matrix

C0_ext = ctrb(A_ext, B_ext);      %C extended matrix
rank_C0_ext = rank(C0_ext);

p_ext = rho+[0, -0.1, -0.2, -0.3, -0.4];    %poles of the extended system
K_ext = place(A_ext, B_ext, p_ext);           %pole placement gain

p_closeloop_ext = eig(A_ext-B_ext*K_ext); % poles of the extended closed loop system

%% SIMULINK MODEL INITILIZATION

x_0 = x_bar+0.01;           %integral initialization

Kx = -K_ext(1, 1:4) ;       %state feedback gain
Kv = -K_ext(1, 5) ;          %additional integrator gain


%% RESPONSE WITH K

x_init= [0, 0, 0, 0]';

sys = ss(Alin - Blin*K, zeros(4,1), eye(4), Dlin);
[y, t]= initial(sys, x_init);

for j= 1:1:length(t)
    u(j,:)= -K*y(j,:)';
end

% uncomment to plot
% figure(1)
% plot(t, u);
% grid on


%% PID for the linearized system

s = tf('s');
G_pp = Clin*(s*eye(4)-(Alin-Blin*K))^-1*Blin;
wc = 1; % rad/s
R_pp = 1/s;
L_pp = R_pp*G_pp;
kp = 1/abs(freqresp(L_pp,wc));
R_pp = kp/s;
L_pp = R_pp*G_pp;
phi_c = angle(freqresp(L_pp,wc));
phi_m = (pi + phi_c)*180/pi;

% uncomment to plot
figure(1)
bode(L_pp); grid on; title('L');
figure(2)
bode(R_pp); grid on; title('R');

%% ------------------------------------------------
%% OPTIMAL CONTROL : LQ 

% Tunable matrix intialization
Q = 0.1*eye(4);
R = 1*diag(1);

[K_lq, S, P] = lqr(Alin ,Blin, Q, R);

% K_lq response

A_lq = Alin-Blin*K_lq ;
B_lq = zeros(4,1) ;
C_lq = Clin ;
D_lq = 0;

% plot the evolution of the state x2 (control variable)

% %linearised system
% dxdt1 = @(t,y)mysystemode(t,y, A_ext-B_ext*K_ext) ;
% tspan1 = [0 5] ;
% y01 = ones(5,1);
% [t1,y1] = ode45(dxdt1, tspan1, y01) ;
% 
% %system with Lq
% dxdt2 = @(t,y)mysystemode(t,y, A_lq) ;
% tspan2 = [0 5] ;
% y02 = ones(4,1);
% [t2,y2] = ode45(dxdt2, tspan2, y02) ;


% uncomment to plot
% figure(3)
% plot(t1, y1(:, 2)) ;
% grid on
% xlabel('Time');
% ylabel('x2 - linearised system');
% title('Evolution of state 2 - pole placement with integrators');
% 
% figure(4)
% plot(t2, y2(:, 2),'c');
% grid on
% xlabel('Time');
% ylabel('x2');
% title('Evolution of state 2 - LQ control');

%% LQ with the system enlargment

%now we want to enlarge the system in order to add an integrator 

A_tilde = [Alin zeros(4,1) ; -Clin zeros(1,1)] ;
B_tilde = [Blin; 0];

r =  1/0.024; 
q =  1;

disp(q/r);

%tunable matrixes 
Q_lq = q*eye(5) ; %1 
R_lq = r*eye(1) ; %1000

% check reachability and observability 
rank(ctrb(A_tilde, B_tilde));
C_q = sqrt(Q_lq) ;
rank(obsv(A_tilde, C_q));

% P_lq = are(A_tilde, B_tilde*R_lq*B_tilde', Q_lq) 
% k = R_lq^(-1)*B_tilde'*P_lq;
% k_x = k(:, 1:4) ;
% eig(Alin-Blin*k_x) 
% 
% return 

% computation of the gain for l1
[K_lqe, S, P] = lqr(A_tilde, B_tilde, Q_lq, R_lq) ;

K_lqx = K_lqe(:,1:4) ;
K_lqeta = K_lqe(:, 5) ;

eig(Alin-Blin*K_lqx);
 

% %system with Lq
% dxdt3 = @(t,y)mysystemode(t,y, A_tilde-B_tilde*K_lqe) ;
% tspan3 = [0 5] ;
% y03 = ones(5,1);
% [t3,y3] = ode45(dxdt3, tspan3, y03) ;
% 
% figure(5)
% plot(t3, y3(:, 2),'b');
% grid on
% xlabel('Time');
% ylabel('x2');
% title('Evolution of state 2 - LQ control with integrators');



%% H2 PROBLEM 

n = 5 ; %n of states
p = 1 ; %n of inputs
m = 1 ; %number of outpuRts

q =     1; %0.01
r =      100;  %10
qt =      0 ; 
rt =      0 ;

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

P_H2 = are(A, B2*D12'*D12*B2', C1'*C1);
K_H2 = (D12'*D12)*B2'*P_H2;

K_H2_X = K_H2(:,1:4) ;
K_H2_V = K_H2(:, 5) ;
eig(Alin-Blin*K_H2_X)

return

%% ----------------------------------------------------------------------------------------------
%% FEEDBACK LINEARIZATION

%% NONLINEAR SYSTEM

%0. Nonlinear system + eq. computation
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

[x1_bar, x3_bar, x4_bar,u_bar] = solve([eqn1,eqn2,eqn3,eqn4],[x1(t),x3(t),x4(t),u(t)])
 
% 1. system definition
syms x1(t) x2(t) x3(t) x4(t) u(t) y(t) % variables
syms  Bm Jl Jm m g l k %system parameters, comment if you want the numerical expression
x1_dot = -Bl/Jl*x1(t)-k/Jl*x2(t)-m*g*l/Jl*cos(x2(t))+k/Jl*x4(t);
x2_dot = x1(t);
x3_dot = k/Jm*x2(t)-Bm/Jm*x3(t)-k/Jm*x4(t)+u(t)/Jm;
x4_dot = x3(t);
y(t) = pi/4-x2(t);

%2. checking the relative degree of the system
y_d = -x2_dot;
y_dd = -x1_dot;
y_ddd = diff(y_dd,t);
y_ddd  = (k*x2_dot)/Jl - (k*x4_dot)/Jl - (g*l*m*sin(x2(t))*x2_dot)/Jl;
y_dddd = diff(y_ddd,t);

%substituting diff(x,t) in order to include it in the simulink scheme

y_dddd = (k*x1_dot)/Jl - (k*x3_dot)/Jl - (g*l*m*sin(x2(t))*x1_dot)/Jl - (g*l*m*cos(x2(t))*x1(t)*x2_dot)/Jl;

%%
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

%%
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
%%

%3. Change of coordinates to bring equilibrium to [0, 0]
syms dx1(t) dx2(t) dx3(t) dx4(t)
dx1_dot = subs(x1_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx2_dot = subs(x2_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx3_dot = subs(x3_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx4_dot = subs(x4_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dy = subs(y(t), [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);

%4. Diffeomorphism
x1_tilde = dy;
x2_tilde = -dx1 ;
%x3_tilde = diff(x2_tilde,t);
x3_tilde = dx1_dot ;
%x4_tilde = diff(x3_tilde, t);
x4_tilde = (k*dx4_dot)/Jl - (k*dx2_dot)/Jl - (Bl*x1_dot)/Jl + (g*l*m*sin(pi/4 + dx2(t))*dx2_dot)/Jl;
%%

% check if it's a diffeomorphism. We need to compute the Jacobean and check
% if the Jacobea is full rank
f = [x1_tilde ; x2_tilde ; x3_tilde ; x4_tilde] ;

% compute the Jacobean
J = jacobian(f, [dx1, dx2, dx3, dx4]);
rank(J);
%%
%6. inverse relationship phi(^-1)
syms dx1 dx2 dx3 dx4 dx1_t dx2_t dx3_t dx4_t
eqn1 = dx1_t == -dx2;
eqn2 = dx2_t == -dx1;
eqn3 = dx3_t == k*(dx2+pi/4)/Jl+m*g*l*cos(dx2+pi/4)/Jl-k*(dx4+x4_bar)/Jl;
eqn4 = dx4_t == k/Jl*dx1-m*g*l/Jl*sqrt(2)/2*(sin(dx2)+cos(dx2))*dx1-k/Jl*dx3;
sys    = [eqn1, eqn2, eqn3, eqn4];
[x1_t, x2_t, x3_t, x4_t]       = solve(sys,[dx1,dx2,dx3,dx4]);

%6b. verify that the equilibrium is in [0,0,0,0]
x1_t = subs(x1_t,[dx1_t,dx2_t,dx3_t,dx4_t], [0,0,0,0]);
x2_t = subs(x2_t,[dx1_t,dx2_t,dx3_t,dx4_t], [0,0,0,0]);
x3_t = subs(x3_t,[dx1_t,dx2_t,dx3_t,dx4_t], [0,0,0,0]);
x4_t = subs(x4_t,[dx1_t,dx2_t,dx3_t,dx4_t], [0,0,0,0]);

%7. system in normal canonical form
syms x2_t x3_t x4_t
x1_t_dot = x2_t;
x2_t_dot = x3_t;
x3_t_dot = x4_t;
x4_t_dot = y_dddd;

y_dddd = (k*x1_dot)/Jl - (k*x3_dot)/Jl - (g*l*m*sin(x2(t))*x1_dot)/Jl - (g*l*m*cos(x2(t))*x1(t)*x2_dot)/Jl;
v = y_dddd;

 %8. Finding u
syms x1 x2 x3 x4 u v_sim k Bm Jm Jl g l m 
v = (g*l*m*sin(x2)*((k*x2)/Jl - (k*x4)/Jl + (g*l*m*cos(x2))/Jl))/Jl - (k*(u/Jm - (Bm*x3)/Jm + (k*x2)/Jm - (k*x4)/Jm))/Jl - (k*((k*x2)/Jl - (k*x4)/Jl + (g*l*m*cos(x2))/Jl))/Jl - (g*l*m*cos(x2)*x1^2)/Jl;
eqn = v_sim == v;
U = solve(eqn, u);

 %% PARAMETERS FOR SIMULINK SIMULATION

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

%% FEEDBACK GAIN TUNING Kv_lin

A_tilde = [0 1 0 0;
           0,0,1,0;
           0,0,0,1;
           0,0,0,0];
B_tilde = [0,0,0,1]';
C_tilde = [1 0 0 0];
p_tilde = -50+[-0, -1, -2, -3]; % poles choosen for pole placement

Kv_tilde = place(A_tilde, B_tilde, p_tilde);
p_tilde_closeloop= eig(A_tilde-B_tilde*Kv_tilde);

%% PID for feedback linearization

s = tf('s');
G_v =C_tilde*(s*eye(4)-(A_tilde-B_tilde*Kv_tilde))^-1*B_tilde;
wc = 3; % rad/s
R_v = 1/s;
L_v = R_v*G_v;
kp_v = 1/abs(freqresp(L_v,wc));
R_v = kp_v/s;
L_v = R_v*G_v;
phi_c = angle(freqresp(L_v,wc));
phi_m = (pi + phi_c)*180/pi

% uncomment to plot
figure(1)
bode(L_v); grid on; title('L_v');
figure(2)
bode(R_v); grid on; title('R_v');
