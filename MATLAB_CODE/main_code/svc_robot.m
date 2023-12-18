% clc
% close all
% clear all

% structure variable controller done on the system linearized through
% feedback linearization

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

q = 5;
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

