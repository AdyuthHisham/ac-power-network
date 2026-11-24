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

% Define parameters for the input signal
DeltaA = 1000;       % Amplitude (1 kV)
Deltaf = 5;          % Frequency difference (Hz)
Deltaomega = 2 * pi * Deltaf; % Angular frequency (rad/s)

% Time vector
t = linspace(0, 0.25, 5000); % Simulate for 1 second

% Define input signal
u1 = DeltaA * sin(Deltaomega * t); % Input 1
u2 = DeltaA * cos(Deltaomega * t); % Input 2
u = [u1; u2]'; % Combine inputs for simulation

% Define initial conditions for the states (varying x4(0))
initial_conditions = [-1000, -750, 0, 250, 1000]; % Different initial Vload values
% Number of initial conditions
num_conditions = length(initial_conditions);

% Preallocate matrices for time and output (assuming equal length for all time vectors)
% Note: You might need to adapt this if the time vectors have different lengths
responses_time = zeros(num_conditions, length(t));  % Preallocate time matrix
responses_output = zeros(num_conditions, length(t)); % Preallocate output matrix

% Create a figure to hold the plots
figure;
hold on;

% Loop through each initial condition
for i = 1:num_conditions
    x4_0 = initial_conditions(i); % Get current initial condition
    
    % Set the initial state vector (only x4 non-zero)
    x0 = [0; 0; 0; x4_0];

    % Simulate the system response to the inputs
    [y, t_out, x] = lsim(sys_lin, u, t, x0);
    
    % Store the results in the preallocated arrays
    responses_time(i, :) = t_out;  % Store time vector
    responses_output(i, :) = y;   % Store output (load voltage)
    
    % Plot the output response for the current initial condition
    plot(t_out, y, 'DisplayName', ['x4(0) = ' num2str(x4_0) ' V']);
end

% Formatting the plot
xlabel('Time (s)');
ylabel('Load Voltage Response (V_{load})');
title('Load Voltage Response for Varying Initial Conditions');
legend('show'); % Show the legend
grid on;