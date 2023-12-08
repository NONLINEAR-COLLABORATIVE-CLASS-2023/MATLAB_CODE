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

p_ext = [-10, -10.1, -10.2, -10.3, -10.4];
K_ext = place(A_ext, B_ext, p_ext);

p_closeloop_ext = eig(A_ext-B_ext*K_ext);

%% RESPONSE WITH K


x_init= [1, 2, 3, 4]';

sys = ss(Alin - Blin*K, zeros(4,1), eye(4), Dlin);
[y, t]= initial(sys, x_init);

for j= 1:1:length(t)
    u(j,:)= -K*y(j,:)';
end

figure(1)
plot(t, u);

%% RESPONSE WITH K_ext

x_init_ext= [1, 2, 3, 4, 5]';


sys_ext = ss(A_ext - B_ext*K_ext, zeros(5, 1), eye(5), Dlin);
[y_ext, t]= initial(sys_ext, x_init_ext);

for j= 1:1:length(t)
    u_ext(j,:)= -K_ext*y_ext(j,:)';
end

figure(2)
plot(t, u_ext);
