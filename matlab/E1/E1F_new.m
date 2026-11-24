clear;
close all;
clc;

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
A = jacobian(f, x); % Partial derivatives of f(x, u) w.r.t. x
B = jacobian(f, u); % Partial derivatives of f(x, u) w.r.t. u
C = jacobian(g, x); % Partial derivatives of g(x, u) w.r.t. x
D = jacobian(g, u); % Partial derivatives of g(x, u) w.r.t. u

% Evaluate the Jacobians at the operating point
x0 = [0; 0; 0; 0]; % Assume equilibrium point for states
u0 = [0; 0];       % Assume equilibrium point for inputs
A_lin = subs(A, [x; u], [x0; u0]);
B_lin = subs(B, [x; u], [x0; u0]);
C_lin = subs(C, [x; u], [x0; u0]);
D_lin = subs(D, [x; u], [x0; u0]);

% Convert symbolic expressions to numeric values
A_lin = double(A_lin);
B_lin = double(B_lin);
C_lin = double(C_lin);
D_lin = double(D_lin);

% Create the linearized state-space system
sys_lin = ss(A_lin, B_lin, C_lin, D_lin);
G = tf(sys_lin);

% % Poles and transmission zeros
% poles = pole(G);
% tzeros = zero(G);
% 
% % Display poles and zeros
% disp('Poles:');
% disp(poles);
% disp('Transmission zeros:');
% disp(tzeros);

% Check for minimum-phase system (poles and zeros in the left half-plane)
% is_minimum_phase = all(real(poles) < 0) && all(real(tzeros) < 0);
% if is_minimum_phase
%     disp('The system is minimum-phase.');
% else
%     disp('The system is not minimum-phase.');
% end

% Plot the pole-zero map
% figure;
% pzmap(sys_lin);
% title('Pole-Zero Map');

% Plot the zeros on the complex plane
[zeros_, poles_, ~] = zpkdata(G)

% % Plot the zeros on the complex plane
% plot(real(zeros_), imag(zeros_), 'o');
% hold on;
% 
% % Plot the poles on the complex plane (optional)
% plot(real(poles_), imag(poles_), 'x');
% 
% % Add axis labels
% xlabel('Real Axis (Seconds^{-1})');
% ylabel('Imaginary Axis (Seconds^{-1})');
% 
% title('Pole-Zero Map');
% hold off;

% Check for pole-zero cancellation
sys_tf_simplified = minreal(sys_lin);
if isequal(sys_lin, sys_tf_simplified)
    disp('No pole-zero cancellation occurred.');
else
    disp('Pole-zero cancellation occurred.');
end


