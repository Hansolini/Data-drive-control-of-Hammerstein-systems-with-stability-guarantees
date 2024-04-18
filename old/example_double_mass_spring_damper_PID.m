%% Descrpition
% This script is ment to demonstrate a case where performing VRFT on a
% Hammerstein system may result in instability due to underestimation of
% the real gain of the system.

% Questions:
% - Why are some results negative once in a while?

%% Some general settings
clc; close;
rng(2);

%% Parameters
Ts = 0.1;
T = 0:Ts:500;
noise_var = 0.025;

%% Nonlinearity
alpha = 0.5;
sat = @(u) min(max(u,-alpha),alpha);
f = @(u) sat(u);

%% System
% Transfer function open loop
m1 = 1; m2 = 0.5;
c1 = 0.2; c2 = 0.5;
k1 = 1; k2 = 0.5;

s = tf('s');
G = minreal((m1*s^2+(c1+c2)*s+(k1+k2))/((m1*s^2+(c1+c2)*s+k1+k2)*(m2*s^2+c2*s+k2)-(k2+c2*s)^2));

% State space representation (augmented to include PID)
sys = ss(G);
A = sys.A; B = sys.B; C = sys.C;

% Funcitons for simulation
dx   = @(x,u) A*x + B*u;
dx_f = @(x,u) A*x + B*f(u);

% VRFT settings
% Reference controller (and model)
kp_ideal = 0.03;
ki_ideal = 0.12;
kd_ideal = 0.12 ;
K = (kp_ideal + ki_ideal/s + kd_ideal*s);

Gm = margin(K*G);

alpha = 1;%1/Gm;
M = feedback(1/alpha*K*G,1);

clc;
disp(evalfr(M,0))
disp(eig(M))
disp(alpha)

% Filter
W = c2d(tf(1,[0.3, 1]),Ts);

% Control strucuter
Bf = [
    tf(1,1,Ts, 'variable','z^-1');
    tf(Ts*[1 1],2*[1 -1],Ts,'variable','z^-1');
    tf([2 -2],Ts*[3 -1],Ts,'variable','z^-1')
];


figure(1);
clf; grid on; hold on;
bode(G)
margin(minreal(alpha*K*G))

figure(2);
clf; grid on; hold on;
step(M)
step(feedback((0.0683 + 0.1105/s + 0.2307*s)*G,1))
%% Sample open loop data
% Additional simulation parameters
x0 = zeros(size(A,1),1);

% Input data
U = 1-2*(1-(square(2 * pi * T / 120) + 1) / 2)';

% Noise free sample data
[~,x_sample_data_linear] = ode45(@(t,x) dx(x,interp1(T,U,t)),T,x0);
y_sample_data_linear = (C*x_sample_data_linear')';

[~,x_sample_data] = ode45(@(t,x) dx_f(x,interp1(T,U,t)),T,x0);
y_sample_data = (C*x_sample_data')';

%%
syms x(t) t a
dx = diff(x,t);
ddx = diff(dx,t);

eq = ddx + a*dx + x == 0;

solve(eq,x(t))

%%

%% VRFT convergence comparison
NUM_DATA = 1000:200:length(T);
n_batches = length(NUM_DATA);
n_samples = 1;

THETA = zeros(3,n_batches,n_samples);
THETA_LINEAR = zeros(3,n_batches,n_samples);
IS_STABLE = zeros(1,n_batches,n_samples);
IS_STABLE_LINEAR = zeros(1,n_batches,n_samples);

% Plot to show progress
figure(1);
clf; grid on; hold on;
xlabel('Data')
ylabel('Percentage unstable')

for i = 1:n_batches
    for j = 1:n_samples
        % Noise
        noise = sqrt(noise_var) * randn(size(T))';
 
        y = y_sample_data + noise;
        y_linear = y_sample_data_linear + noise;
    
        [~, theta]        = VRFT1_ry(U(1:NUM_DATA(i)), y(1:NUM_DATA(i)),        c2d(M,Ts), Bf, W, 4, []);
        [~, theta_linear] = VRFT1_ry(U(1:NUM_DATA(i)), y_linear(1:NUM_DATA(i)), c2d(M,Ts), Bf, W, 4, []);
    
        THETA(:,i,j) = theta;
        THETA_LINEAR(:,i,j) = theta_linear;

        % Compute closed loop
        M_vrft        = feedback((THETA(1,i,j) + THETA(2,i,j)/s + THETA(3,i,j)*s)*G,1);
        M_vrft_linear = feedback((THETA_LINEAR(1,i,j) + THETA_LINEAR(2,i,j)/s + THETA_LINEAR(3,i,j)*s)*G,1);

        IS_STABLE(1,i,j) = all(real(eig(M_vrft)) < 0);
        IS_STABLE_LINEAR(1,i,j) = all(real(eig(M_vrft_linear)) < 0);

        % Progress update
        percentage = ((i - 1)*n_samples + j)/(n_batches*n_samples)*100;
        
        clc;
        fprintf('TP | CP: %3.2f%% | %3.2f%% \n', percentage, j/n_samples*100)
    end

    % Plotting
    clf; hold on;
    bar(NUM_DATA, [sum(~IS_STABLE_LINEAR,3); sum(~IS_STABLE,3)]/(2*n_samples), 'stacked')
    plot(U, 'k--');
    legend({'Linear','Hammerstein'},'location','best')
    drawnow;

    % Saving data
    save('data\convergence.mat');
end



