clear all
clc
close all

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
        1           0                           0             0;
        0           k/Jm                     -Bm/Jm      -k/Jm;
        0           0                           1           0];

Blin = [0; 0; 1/Jm; 0];

Clin= [0    1   0   0];             %applichiamo a x2(theta l) perchè abbiamo solo 1 input(2 output)

Dlin = zeros(1);

%% POLES OF THE SYSTEM

p_Alin= eig(Alin);

%% POLE PLACEMENT (FOR STABILIZATION)

C0 = ctrb(Alin, Blin);
rank_C0 = rank(C0);         %rank of controllability matrix

p = [-10, -10.1, -10.2, -10.3];
K = place(Alin, Blin, p);

p_closeloop= eig(Alin-Blin*K);

%% POLE PLACEMENTE WITH INTEGRATORS

A_ext= [Alin, zeros(4, 1);
        -Clin, zeros(1)];

B_ext = [Blin; zeros(1)];

C0_ext = ctrb(A_ext, B_ext);
rank_C0_ext = rank(C0_ext);

p_ext = 2*[-10, -10.1, -10.2, -10.3, -10.4];
K_ext = place(A_ext, B_ext, p_ext);

p_closeloop_ext = eig(A_ext-B_ext*K_ext);

%% RESPONSE WITH K

x_init= [1, 2, 3, 4]';

sys = ss(Alin - Blin*K, zeros(4,1), eye(4), Dlin);
[y, t]= initial(sys, x_init);

for j= 1:1:length(t)
    u(j,:)= -K*y(j,:)';
end

% figure(1)
% plot(t, u);
% grid on

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

%% SIMULATION IN MATLAB

x_0 = x_bar+0.01;

Kx = -K_ext(1, 1:4) ;
Kv = -K_ext(1, 5) ;

%% ------------------------------------------------
%% OPTIMAL CONTROL : LQ 

Q = 10*eye(4);
R = 0.1*diag([1]);

[K_lq, S, P] = lqr(Alin ,Blin, Q, R);


% risposta con K_lq

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

%now we want to enlarge the system in order toput the integrator 

A_tilde = [Alin zeros(4,1) ; -Clin zeros(1,1)] ;
B_tilde = [Blin ; 0 ] ;

%now let's design the LQ regulator 

Q_lq = eye(5) ;
R_lq = 10000*eye(1) ;

% Before apply LQif control we need to ceck that 
% 1. (A_tilde, B_tilde) is reachable 

rank(ctrb(A_tilde, B_tilde))

%2. (A_tilde, C_q) is observable 

C_q = sqrt(Q_lq) ;

rank(obsv(A_tilde, C_q))

%so we verify the two conditions 

[K_lqe, S, P] = lqr(A_tilde, B_tilde, Q_lq, R_lq) ;

%now we can separate the matrix K_lq

K_lqx = K_lqe(:,1:4) ;
K_lqeta = K_lqe(:, 5) ;


%% 




