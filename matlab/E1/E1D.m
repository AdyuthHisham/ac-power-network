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
A = jacobian(f, x); 
B = jacobian(f, u); 
C = jacobian(g, x); 
D = jacobian(g, u); 

% Evaluate the Jacobians at the operating point
x0 = [0; 0; 0; 0]; % Assume equilibrium point for states
u0 = [0; 0];       % Assume equilibrium point for inputs
A_l = subs(A, [x; u], [x0; u0]);
B_l = subs(B, [x; u], [x0; u0]);
C_l = subs(C, [x; u], [x0; u0]);
D_l = subs(D, [x; u], [x0; u0]);

% Convert symbolic expressions to numeric values
A_l = double(A_l);
B_l = double(B_l);
C_l = double(C_l);
D_l = double(D_l);

% Create the linearized state-space model
sys_l = ss(A_l, B_l, C_l, D_l);

% Time vector for simulation
t = linspace(0, 0.1, 1000); % 0 to 0.5 seconds

% Define a step input for both u1 and u2
u1_step = ones(size(t)); % Step input for u1
u2_step = ones(size(t));
u_step = [u1_step; u2_step]'; % Combine inputs for simulation

% Simulate the response of the linearized system
[y, t_out] = lsim(sys_l, u_step, t);

% Plot the output response of the linearized system
figure;
plot(t_out, y, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Output Response (V_{load})');
title('Response of the Linearized System to Step Inputs');
grid on;
