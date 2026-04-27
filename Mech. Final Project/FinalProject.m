%% System Definition

syms x1 x2 x3 x4 N Tm
g = 9.81;
%lp = 0.019;
lp = 0.025;
mp = 0.068 + 0.047;
mw = 0.034;
Rw = 0.022;
Iw = (1/2)*mw*Rw^2;
%Ip = 0;
Ip = (1/12)*mp*lp^2 + mp*(0.08-lp)^2;
Rm = 13;
Kb = 0.155;
Kt = 0.1;

A_sys = [mw     0   -1  0   0   1;
    0       0   0   -1  -1  0;
    -Iw/Rw  0   0   0   0   Rw;
    mp  -mp*lp*cos(x3) 1 0 0 0;
    0   -mp*lp*sin(x3) 0 1 0 0;
    0   -Ip      -lp*cos(x3) -lp*sin(x3) 0 0];

B_sys = [0; 
    -mw*g; 
    Tm; 
    -mp*lp*x4^2*sin(x3);
    mp*lp*x4^2*cos(x3)-mp*g;
    Tm];

z = A_sys\B_sys;

simplify(z(1))
simplify(z(2))

%% Get A and B

%Tm = ((Kt/Rm)*(V-Kb*(-x(2)/Rw-x(4))));


% define state derivatives
xdot = @(x,u) [
    x(2);  
    double(subs(z(1),[x3,x4,Tm],[x(3),x(4),((Kt/Rm)*(u-Kb*(-x(2)/Rw-x(4))))]));
    x(4);
    double(subs(z(2),[x3,x4,Tm],[x(3),x(4),((Kt/Rm)*(u-Kb*(-x(2)/Rw-x(4))))]));
];

% define initial conditions
x0 = [0; 0; 0; 0];
u0 = 0;

% small change for central difference calculations
delta = 1e-6;

% find A matrix
A = zeros(4, 4);
for i = 1:4
    dx = zeros(4,1);
    dx(i) = delta;
    A(:, i) = (xdot(x0 + dx, u0) - xdot(x0 - dx, u0)) / (2 * delta);
end

% find B matrix
B = (xdot(x0, u0 + delta) - xdot(x0, u0 - delta)) / (2 * delta);

A
B

%% lqr calculation
R = 1;
Q = [1000 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

K = lqr(A, B, Q, R)