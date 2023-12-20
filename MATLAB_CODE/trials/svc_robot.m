% clc
% close all
% clear all

% structure variable controller done on the system linearized through
% feedback linearization
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

% 2.1 EQUILIBRIUM
x1_bar = 0;
x2_bar = pi/4;
x3_bar = 0;
x4_bar = m*g*l*cos(x2_bar)/k + x2_bar;
u_bar= m*g*l*cos(x2_bar);
x_bar = [x1_bar x2_bar x3_bar x4_bar]';
y_bar = -x2_bar;

dx_0 = [0;1;0;1]*0.1;  


%% 


A= [0      1       0       0;
    0      0       1       0;
    0      0       0       1;
    0      0       0       0;];

B= [0   0   0   1]';

C= [1   0   0   0];

p = +10;                             %poli
beta_prime= [p^3 3*p^2 3*p 1];
beta3= beta_prime(1);
b4= C(1);
gamma= beta3/b4;
alfa_prime= beta_prime*A;

%% parameters q, r

q = 20;
r = 0.7;

%% hysteresis parameters

BI= 0.02;
M= 1;

%% enlarged system(with integrator)

F= [A           B;
    zeros(1,4)  0];
G= [zeros(4,1); 1];
H= [C   0];


Mr = ctrb(F, G);

a_tilde= poly(F);
a= flip(-a_tilde);
A_tilde= [0     1       0       0       0;
          0     0       1       0       0;
          0     0       0       1       0;
          0     0       0       0       1;
          a(1:5)];

B_tilde= [0 0 0 0 1]';
Mr_tilde = ctrb(A_tilde, B_tilde);

T= Mr_tilde*Mr^-1;

F_tilde= T*F*T^-1;
G_tilde= T*G;
H_tilde= H*T^-1;


%%


p_t = -10;                             %poli
beta_prime_t= flip(poly([p_t p_t p_t p_t]));
beta3_t= beta_prime_t(1);
b4_t= H_tilde(1);
gamma_t= beta3_t/b4_t;
alfa_prime_t= beta_prime_t*F_tilde;


%% LOAD DISTURBANCES(done on system without integrator)

W = 2;
A_new= A - B*alfa_prime;




