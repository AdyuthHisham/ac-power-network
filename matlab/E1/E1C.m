syms t s

%Transmission line 1
r1 = 0.3;     
l1 = 1e-3;        
c1 = 10e-9;    
d1 = 100;

%Transmission line 2
r2 = 0.15;     
l2 = 1e-3;        
c2 = 15e-9;        
d2 = 200;    

%Transmission Line 1, times distance
R1 = r1*d1;
L1 = l1 *d1;
C1 = c1*d1/2;    %Divided by 2, since split between two

%Transmission Line 2, times distance
R2 = r2*d2;
L2 = l2*d2;
C2 = c2*d2/2;

%Variables given in question
RL = 100;
LL = 1;
E = 10000;
f = 50;

% Define the system matrices from the state-space model
A = [-R1/L1 0 0 -1/L1;
    0 -R2/L2 0 -1/L2;
    0 0 -RL/LL 1/LL;
    1/(C1+C2) 1/(C1+C2) -1/(C1+C2) 0];

B = [1/L1 0;
    0 1/L2;
    0 0;
    0 0];

C = [0 0 0 1];
D = [0 0];

% Create the state-space system
sys = ss(A, B, C, D);
trans = tf(sys);

% Define the frequency of the input signal (50 Hz)
w = 2*pi*f;

% Laplace Transforms of the inputs
U1 = (sqrt(2)*E * w) / (s^2 + w^2); % Laplace of u1(t)
U2 = (sqrt(2)*E * s) / (s^2 + w^2); % Laplace of u2(t)

% Solve the system in the S-domain using the state-space system
%X(s) = (SI - A)^-1*Bu(s) --> Y(s) = C*X(s)
Y1_s = C * inv(s*eye(4) - A) * B * U1;   
Y2_s = C * inv(s*eye(4) - A) * B * U2;

% Combine the outputs for the total system response
Y_s = Y1_s + Y2_s;

% Inverse Laplace Transform to get the output in the time domain
Y_t = ilaplace(Y_s, s, t);

% Plot the result
fplot(Y_t, [0, 0.05]); % Plot from t=0 to 0.05 seconds
xlabel('Time (s)');
ylabel('Output Response');

% Define the input signal as a function of time
t = linspace(0, 0.08, 1000); % Time vector
u1_t = sqrt(2)*E * sin(w*t); % Input u1(t)
u2_t = sqrt(2)*E * cos(w*t); % Input u2(t)

% Combine the inputs
u = [u1_t; u2_t]; 

% Simulate the system response using lsim
[y, t_out] = lsim(sys, u', t);  

% Plot the output response
plot(t_out, y);
xlabel('Time (s)');
ylabel('Output Response');