clear;
close all;
clc;

syms t s
syms LL RL real positive
assume(RL > 0 & LL > 0);

% Define symbolic variables for states and inputs
syms x1 x2 x3 x4 u1 u2 real

R1 = 30;
L1 = 0.1;
C1 = (1*10^-6)/2;

R2 = 30;
L2 = 0.2;
C2 = (3*10^-6)/2;

LL = 1;
RL = 100;

E = 10000;
f = 50;
w = 2*pi*f;

% State vector
x = [x1; x2; x3; x4]; % x1 = current through line 1, x2 = current through line 2, ...
u = [u1; u2];          % Inputs

% System dynamics (derived from earlier state-space equations)
dx1 = (-R1/L1)*x1 - (1/L1)*x4 + (1/L1)*u1;
dx2 = (-R2/L2)*x2 - (1/L2)*x4 + (1/L2)*u2;
dx3 = (-RL/LL)*x3 + (1/LL)*x4;
dx4 = (1/(C1 + C2)) * (x1 + x2 - x3);

% Output equation
y = x4; % Voltage across the load

% State equations
f = [dx1; dx2; dx3; dx4]; % f(x, u)
g = y;                    % Output equation

% Compute Jacobians for linearization
A_pre = jacobian(f, x); % Partial derivatives of f(x, u) w.r.t. x
B_pre = jacobian(f, u); % Partial derivatives of f(x, u) w.r.t. u
C_pre = jacobian(g, x); % Partial derivatives of g(x, u) w.r.t. x
D_pre = jacobian(g, u); % Partial derivatives of g(x, u) w.r.t. u

% Evaluate the Jacobians at the operating point
x0 = [0; 0; 0; 0]; % Assume equilibrium point for states
u0 = [0; 0];       % Assume equilibrium point for inputs
A = subs(A_pre, [x; u], [x0; u0]);
B = subs(B_pre, [x; u], [x0; u0]);
C = subs(C_pre, [x; u], [x0; u0]);
D = subs(D_pre, [x; u], [x0; u0]);


controllability_matrix = [B, A*B, A^2*B, A^3*B]; % Controllability matrix
observeability_matrix = [C; C*A; C*A^2; C*A^3];

RREF_controllability_matrix = rref(controllability_matrix);
RREF_observeability_matrix = rref(observeability_matrix);

is_controllable = rank(RREF_controllability_matrix) == size(A, 1);
is_observable = rank(RREF_observeability_matrix) == size(A, 1);

% Display results
disp(['Controllable: ', num2str(is_controllable)]);
disp(['Observable: ', num2str(is_observable)]);

