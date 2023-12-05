%% COLLABORATIVE CLASS

%% DATA INITIALIZATION
k = 0.8;
Jm = 4e-04;
Jl = Jm;
Bm = 0.015;
Bl = 0;
m = 0.3;
l = 0.3;
g = 9.81;
u_bar = sqrt(2)/2*m*g*l;
x_bar = [0 pi/4 0 u_bar/k+pi/4];
x_0 = x_bar+0.01;
Kx = -[-0.00024 0.758 0.0054 -0.5596 ];
Kv = -[-0.0221];