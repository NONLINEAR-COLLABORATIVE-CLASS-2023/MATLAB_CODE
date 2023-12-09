clear all;
clc;
close all;
%% COLLABORATIVE CLASS

%% PARAMETERS

k = 0.8;
Jm = 4e-4;
Jl= 4e-4;
Bm = 0.015;
Bl = 0;
m = 0.3;
l = 0.3;
g = 9.81;

%% EQUILIBRIUM

x1_bar = 0;
x2_bar = pi/4;
x3_bar = 0;
x4_bar = m*g*l*cos(x2_bar)/k + x2_bar;
u_bar= m*g*l*cos(x2_bar);
x_bar = [x1_bar x2_bar x3_bar x4_bar] ;

%% LINEARIZATION

Alin= [0    -k/Jl + m*g*l*sin(x2_bar)/Jl        0            k/Jl;
        1                        0                           0             0;
        0                      k/Jm                     -Bm/Jm      -k/Jm;
        0                       0                           1              0];

Blin = [0; 0; 1/Jm; 0];

Clin= [0    1   0   0];             %applichiamo a x2(theta l) perchè abbiamo solo 1 input(2 output)

Dlin = zeros(1);

%% POLES OF THE SYSTEM

poles_lin= eig(Alin);

%% POLE PLACEMENT (FOR STABILIZATION)

C0 = ctrb(Alin, Blin);          %controllability matrix
rank_C0 = rank(C0);          %rank of controllability matrix

p = [-10, -10.1, -10.2, -10.3]; % poles choosen for pole placement
K = place(Alin, Blin, p);

p_closeloop= eig(Alin-Blin*K);

%% POLE PLACEMENT WITH INTEGRATORS

A_ext= [Alin, zeros(4, 1);            %A extended matrix
       -Clin, zeros(1)];

B_ext = [Blin; zeros(1)];               %B extended matrix

C0_ext = ctrb(A_ext, B_ext);      %C extended matrix
rank_C0_ext = rank(C0_ext);

p_ext = 2*[-10, -10.1, -10.2, -10.3, -10.4];    %poles of the extended system
K_ext = place(A_ext, B_ext, p_ext);           %pole placement gain

p_closeloop_ext = eig(A_ext-B_ext*K_ext); % poles of the extended closed loop system

%% RESPONSE WITH K

x_init= [0, 0, 0, 0]';

sys = ss(Alin - Blin*K, zeros(4,1), eye(4), Dlin);
[y, t]= initial(sys, x_init);

for j= 1:1:length(t)
    u(j,:)= -K*y(j,:)';
end

figure(1)
plot(t, u);
grid on


%% PID 

s = tf('s');
G =Clin*(s*eye(4)-(Alin-Blin*K))^-1*Blin;

kp=10;
R2= kp*G^(-1)/(s*(s/500+1)^4); 
L2=G*R2;

%[N2, D2]= tfdata(R2);

figure
%bode(L); grid on; title('L');
bode(L2); grid on; title('L2');

%% RESPONSE WITH K_ext

x_init_ext= [0.5 0.5 0.5 0.5 1] ;

sys_ext = ss(A_ext - B_ext*K_ext, zeros(5, 1), eye(5), Dlin);
[y_ext, t]= initial(sys_ext, x_init_ext);

for j= 1:1:length(t)
    u_ext(j,:)= -K_ext*y_ext(j,:)';
end

figure(2)
plot(t, u_ext);
grid on

%% SIMULINK MODEL INITILIZATION

x_0 = x_bar+0.01;           %integral initialization

Kx = -K_ext(1, 1:4) ;       %state feedback gain
Kv = -K_ext(1, 5) ;          %additional integrator gain

%% ------------------------------------------------
%% OPTIMAL CONTROL : LQ 

% Tunable matrix intialization
Q = 10*eye(4);
R = 0.1*diag(1);

[K_lq, S, P] = lqr(Alin ,Blin, Q, R);

% K_lq response

A_lq = Alin-Blin*K ;
B_lq = zeros(4,1) ;
C_lq = Clin ;
D_lq = 0;

% plot the evolution of the state x2 (control variable)

%linearised system
dxdt1 = @(t,y)mysystemode(t,y, A_ext-B_ext*K_ext) ;
tspan1 = [0 5] ;
y01 = ones(5,1);
[t1,y1] = ode45(dxdt1, tspan1, y01) ;

%system with Lq
dxdt2 = @(t,y)mysystemode(t,y, A_lq) ;
tspan2 = [0 5] ;
y02 = ones(4,1);
[t2,y2] = ode45(dxdt2, tspan2, y02) ;

figure(3)
plot(t1, y1(:, 2)) ;
grid on
xlabel('Time');
ylabel('x2 - linearised system');
title('Evolution of state 2 - pole placement with integrators');

figure(4)
plot(t2, y2(:, 2),'c');
grid on
xlabel('Time');
ylabel('x2');
title('Evolution of state 2 - LQ control');

%% LQ with the system enlargment

%now we want to enlarge the system in order to add an integrator 

A_tilde = [Alin zeros(4,1) ; -Clin zeros(1,1)] ;
B_tilde = [Blin ; 0 ] ;

%tunable matrixes 
Q_lq = eye(5) ;
R_lq = 10000*eye(1) ;

% check reachability and observability 
rank(ctrb(A_tilde, B_tilde));
C_q = sqrt(Q_lq) ;
rank(obsv(A_tilde, C_q));

% computation of the gain for l1
[K_lqe, S, P] = lqr(A_tilde, B_tilde, Q_lq, R_lq) ;

K_lqx = K_lqe(:,1:4) ;
K_lqeta = K_lqe(:, 5) ;

%% ----------------------------------------------------------------------------------------------
%% FEEDBACK LINEARIZATION

syms x1(t) x2(t) x3(t) x4(t) u(t) y(t) % variables
syms Bl Bm Jl Jm m g l k %system parameters, comment if you want the numerical expression

% 1. system definition
x1_dot = -Bl/Jl*x1(t)-k/Jl*x2(t)-m*g*l/Jl*cos(x2(t))+k/Jl*x4(t);
x2_dot = x1(t);
x3_dot = k/Jm*x2(t)-Bm/Jm*x3(t)-k/Jm*x4(t)+u/Jm;
x4_dot = x3(t);
y(t) = pi/4-x2(t);

%2. checking the relative degree of the system
y_d = -x2_dot;
y_dd = -x1_dot;
%y_ddd = diff(-x1_dot,t)
y_ddd  = (Bl*x1_dot)/Jl + (k*x2_dot)/Jl - (k*x4_dot)/Jl - (g*l*m*sin(x2(t))*x2_dot)/Jl; %% this code ha to be changed if sth changes in the system
y_dddd = diff(y_ddd,t);
y_dddd = (k*x1_dot)/Jl - (k*x3_dot)/Jl - (Bl*((Bl*x1_dot)/Jl + (k*x2_dot)/Jl - (k*x4_dot)/Jl - (g*l*m*sin(x2(t))*x2_dot)/Jl))/Jl - (g*l*m*sin(x2(t))*x1_dot)/Jl - (g*l*m*cos(x2(t))*x1(t)*x2_dot)/Jl; %% this code ha to be changed if sth changes in the system
y_dddd = simplify(y_dddd);
y_dddd = subs(y_dddd, Bl, 0);

%3. Change of coordinates to bring equilibrium to [0, 0]
syms dx1(t) dx2(t) dx3(t) dx4(t)
dx1_dot = subs(x1_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx2_dot = subs(x2_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx3_dot = subs(x3_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dx4_dot = subs(x4_dot, [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);
dy = subs(y(t), [x1(t),x2(t),x3(t),x4(t)], [dx1(t)+x1_bar, dx2(t)+x2_bar, dx3(t)+x3_bar, dx4(t)+x4_bar]);

%4. diffeomorphism
x1_tilde = dy;
x2_tilde = -dx1 ; %x2_tilde = diff(dy,t)
x3_tilde = diff(dy,t,2);
x4_tilde = diff(dy,t,3);

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
%%v = simplify(y_dddd,'steps',10)
 %8. Finding u


 %% REGULATOR
s = tf('s');
kp = 1;
G = 1/s^4;
R = kp*(s+0.001)^3/(s/50+1)^4;
L = G*R;
bode(L), grid on;

syms s
R = kp*(s+0.001)^3/(s/50+1)^4;
v = ilaplace(R);