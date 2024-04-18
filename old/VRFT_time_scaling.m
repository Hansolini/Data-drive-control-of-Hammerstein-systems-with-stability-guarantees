addpath('VRFT_toolbox_v1\*')

clc;
clear all;

rng(1);

%% Real system
m1 = 1;
m2 = 0.5;
c1 = 0.2;
c2 = 0.5;
k1 = 1;
k2 = 0.5;

s = tf('s');
P = minreal((m1*s^2+(c1+c2)*s+(k1+k2))/((m1*s^2+(c1+c2)*s+k1+k2)*(m2*s^2+c2*s+k2)-(k2+c2*s)^2));

P = 1/(s^2+s+1);
%% Sample data
% Parameters
Ts = 0.1;
t = 0:Ts:500;
var = 2*0.025;

% Data
noise = sqrt(var) * randn(size(t))';
u = (1-(square(2 * pi * t / 120) + 1) / 2)';
y = lsim(P,u,t) + noise;

% Figure
figure(1);
clf; grid on; hold on; box on;
plot(t, u, 'k--', 'linewidth', 2)
plot(t, y, 'r', 'linewidth', 2)

xlim([0,500]);
xlabel('Time')

%% VRFT + time-scaling
% Reference model
Mr = tf(1, [3, 1]);

% Figure
figure(2);
clf; grid on; hold on;
plot(t,step(Mr,t), 'k-', 'LineWidth', 2)
xlim([0,50])

% Run experiments
ALPHA = 1:3:10;
for i = 1:length(ALPHA)
    alpha = ALPHA(i);

    Mr_alpha = tf(1, [3/alpha, 1]);
    Mr_alpha_inv = tf(1, [3*alpha, 1]);
    
    W = c2d(tf(1,[0.3, 1]),Ts);
    B = [
            tf(1, 1, Ts, 'variable', 'z^-1'); 
            tf(Ts*[1 1], 2*[1 -1], Ts, 'variable', 'z^-1');
            tf([2 -2], Ts*[3 -1], Ts, 'variable', 'z^-1');
        ];
    
    % Using the VRFT1_ry function, though its exact implementation is not provided in the snippet
    [~, theta] = VRFT1_ry(u, y, c2d(Mr_alpha_inv,Ts), B, [], 2, []);
    
    % Extract the PI controller parameters from theta
    Kp = theta(1); Ki = theta(2); Kd = theta(3);
    
    % Controller
    C = Kp + 1/s*Ki + s*Kd;
        
    % Computing "real" closed loop
    M_alpha_inv = minreal(feedback(P*C,1));

    % Rescaling closed loop
    K1 = alpha;
    K2 = 1-1/alpha;
    M = minreal(feedback(K1*M_alpha_inv,K2));
    
    % Figure
    figure(2);
    set(gca,'ColorOrderIndex',i)
    plot(t/alpha,step(M_alpha_inv,t), ':', 'LineWidth', 2)
    set(gca,'ColorOrderIndex',i)
    plot(t,step(M,t), '--', 'LineWidth', 2)
    drawnow;
end

% Info:
% - Three possible end results for the closed loop
% -- (1) Mr
% -- (2) "poor" estimate when time-scalind, but perfect rescaling of time
% -- (3) "poor" estimate and "poor" rescaling
% - This plots all cases agains each other, for different values of alpha
% 
% Desired steady-state is achieved in all cases (?), but not knowing the
% full dynamics is what causes problems.
% It also seems like all of the models become non-minimum phase in this
% example.