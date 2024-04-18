%% TODO:
% x Choose a nonlinearity that acts as the saturation given the input, but
% that results in an unstable system
% --> tanh and atanh with limited domain!
% x Simulation of Hammerstein systems in closed loop!

% x Implement the stability constraint for linear systems
% x Implement the Hammerstein constraint
% - Verify with the VRFT toolbox
% - Implement the IV correlation stuff
% -- Not necessary, as we want to show that the condition holds (could also
% just not use noise at all)
% - OBS! This code does not implement a prediction error approach! This is
% actually cruical... This optimization scheme is nonlinear in th
% parameters! We should actually use a PEM approach to make this
% minimization...

rng(4);
addpath(genpath('functions\'))
warning('off', 'all')

%% Initialize system
% Define simulation time and sampling time
t = linspace(0,500,1000)';
dt = max(diff(t));

% Linear model
m1 = 1; m2 = 0.5;
c1 = 0.2; c2 = 0.5;
k1 = 1; k2 = 0.5;

s = tf('s');
G = c2d(minreal((m1*s^2+(c1+c2)*s+(k1+k2))/((m1*s^2+(c1+c2)*s+k1+k2)*(m2*s^2+c2*s+k2)-(k2+c2*s)^2)),dt,'zoh');

% Nonlinearity
a = 1.75;
f = @(u) tanh(a*u)/a;%2/pi*atan(pi/2*u);%(u/3).^3;
fi_hat = @(v) atanh(v);%2/sqrt(pi)*erfinv(v);%tan(a*v)/a;%real(3^3*tanh(a*v)/a).^(1/3);
f_hat  = @(u) tanh(u);%erf(sqrt(pi)/2*u);

% Plotting nonlinearity
u_range = linspace(-2,2,100);
figure(1);
clf; grid on; hold on; box on;

plot(u_range,f(u_range),'b','LineWidth',2)
plot(fi_hat(u_range),u_range,'r--','LineWidth',2)

plot(u_range,f(fi_hat(u_range)),'g','LineWidth',2)
plot(u_range,u_range,'k:','LineWidth',2)

plot(u_range(1:end-1),diff(f(u_range))./diff(u_range), 'm-')
plot(u_range(1:end-1),diff(f_hat(u_range))./diff(u_range), 'k-')

plot(u_range(1:end-1),diff(f(fi_hat(u_range)))./diff(u_range), 'r-')
xline(-1,'linewidth',2)
xline(1,'linewidth',2)
xlim([-1,1])
ylim([-1,1])

%%
alpha = f(fi_hat(1));
%% Sample data
v_hat = square(0.05*t);
v = f(fi_hat(v_hat)); % Using the estimated inverse to sample

noise = 0.05*randn(size(t));
y_linear      = lsim(G,v_hat,t) + noise;
y_hammerstein = lsim(G,v,t) + noise;

% Initialize controller and reference model
beta = [
    tf(1,1,dt);                     % P
    dt/2*tf([1, 1],[1,-1],dt);      % I
    2/dt*tf([1, -1],[3, -1],dt);    % D
    ]; % Tustin transform

beta = c2d(d2c(beta,'tustin'),beta.Ts,'zoh'); % Transform from Tustin to zoh (more numerically stable)

theta_star = [0.1, 0.11, 0.08]';
Kr = transpose(beta)*theta_star;
Mr = minreal(feedback(Kr*G,1));


%%
r = lsim(inv(Kr),v_hat,t);
figure(1);
clf; grid on; hold on;
plot(t,r, 'linewidth',2)

%% Control identification
% Initial condition
theta0 = theta_star.*rand(3,1);

% Parameters for stability constraints
deltaN = 0.8;

% Optimization
options = optimoptions("fmincon","Algorithm","sqp","Display","none","StepTolerance",1e-6);
THETA_linear = fmin_model_free(theta0,beta,Mr,G,v_hat,y_linear,t,deltaN,1,options);
[THETA_hammerstein, exitflag, output] = fmin_model_free(theta0,beta,Mr,G,v_hat,y_hammerstein,t,deltaN,alpha,options);

clc;
if exitflag <= 0
    fprintf('Simulation #%.f failed with exitflag %.f. \n',find(sum(abs(THETA_hammerstein)) == 0), exitflag)
    disp('The following message was recorded:')
    disp('------------------------------------')
    disp(output.message)
    return;
end

%% Results
disp('Ratio of real and estimated parameters:')
disp(theta_star./THETA_hammerstein)

% Derived results
r = sim_virtual_reference(Mr,y_hammerstein,t);
e = r - y_hammerstein;

K0 = transpose(beta)*theta0;
M0 = feedback(K0*G,1);

K = cell(size(THETA_hammerstein,2),1);
M = cell(size(THETA_hammerstein,2),1);
for i = 1:size(THETA_hammerstein,2)
    K{i} = transpose(beta)*THETA_hammerstein(:,i);
    M{i} = feedback(K{i}*G,1);
end

% Figures 
figure(1);
clf; grid on; hold on; box on;

plot(t,y_hammerstein, '-',   'linewidth', 2, 'DisplayName','$y$')
plot(t,v_hat, '--', 'linewidth', 2, 'DisplayName','$u$')
plot(t,r, '-',  'linewidth', 2, 'DisplayName','$r$')
plot(t,e, ':',  'linewidth', 2, 'DisplayName','$e$')

xlabel('$t$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
legend('location','best','interpreter','latex')
title('Sample data')

figure(2);
clf; 

DISPLAYNAMES = {
    '$J_{mr}$';
    '$J_{est}$';
    '$J_{vrft}$';
    '$J_{vrft}$ + ideal con.';
    '$J_{vrft}$ + est. con.';
    '$J_{vrft}$ + est. Hamm. con.';
};

subplot(2,1,1);
grid on; hold on;box on;
xlim([0,50])
ylim([0,2])

subplot(2,1,2);
grid on; hold on;box on;
xlim([0,50])
ylim([0,2])

for i = 3:length(M)
    if i == 4
        continue;
    end
    % Simulate step response
    yl = step_closed_loop_hammerstein(G,K{i},@(x) x,t); % Using linear "nonlinearity" here
    yh = step_closed_loop_hammerstein(G,K{i},@(x) f(fi_hat(x)),t);
    
    % Plot response
    subplot(2,1,1);
    plot(t,yl, '-', 'linewidth', 2*length(M)+2-2*i,'DisplayName',DISPLAYNAMES{i})
    
    subplot(2,1,2);
    plot(t,yh, '-', 'linewidth', 2*length(M)+2-2*i,'DisplayName',DISPLAYNAMES{i})
end

subplot(2,1,1);
plot(t,step(Mr,t), 'k--', 'linewidth', 3, 'displayname', '$M_r$')
xlabel('$t$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
legend('location','best','interpreter','latex')
title('Linear step response');

subplot(2,1,2);
plot(t,step(Mr,t), 'k--', 'linewidth', 3, 'displayname', '$M_r$')
xlabel('$t$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
legend('location','best','interpreter','latex')
title('Hammerstein step response');

save('data/constrained_vrft.mat','G','beta','theta_star','THETA_hammerstein','THETA_linear','Kr','Mr','K','t','v_hat','y_hammerstein','y_linear','f','fi_hat','-mat')
%% Functions
function [THETA, exitflag,output] = fmin_model_free(theta0,beta,Mr,G,u,y,t,deltaN,alpha,options)
    % Initialization
    THETA = zeros(length(beta),6);

    % Optimization
    [THETA(:,1),~,exitflag,output] = fmincon(@(theta) Jmr(theta,beta,Mr,G),theta0,[],[],[],[],zeros(size(theta0)),[],[],options);
    if exitflag <= 0
        return;
    end
    
    [THETA(:,2),~,exitflag,output] = fmincon(@(theta) Jest(theta,beta,Mr,G),theta0,[],[],[],[],zeros(size(theta0)),[],[],options);
    if exitflag <= 0
        return;
    end
    
    [THETA(:,3),~,exitflag,output] = fmincon(@(theta) Jvrft(theta,beta,Mr,u,y,t),theta0,[],[],[],[],zeros(size(theta0)),[],[],options);
    if exitflag <= 0
        return;
    end
    
    [THETA(:,4),~,exitflag,output] = fmincon(@(theta) Jvrft(theta,beta,Mr,u,y,t),theta0,[],[],[],[],zeros(size(theta0)),[],@(theta) ideal_stability_constraint(theta,beta,Mr,G,deltaN),options);
    if exitflag <= 0
        return;
    end
    
    [THETA(:,5),~,exitflag,output] = fmincon(@(theta) Jvrft(theta,beta,Mr,u,y,t),theta0,[],[],[],[],zeros(size(theta0)),[],@(theta) estimated_stability_constraint(theta,beta,Mr,u,y,t,deltaN,1,G),options);
    if exitflag <= 0
        return;
    end

    [THETA(:,6),~,exitflag,output] = fmincon(@(theta) Jvrft(theta,beta,Mr,u,y,t),theta0,[],[],[],[],zeros(size(theta0)),[],@(theta) estimated_stability_constraint(theta,beta,Mr,u,y,t,deltaN,alpha,G),options);
    if exitflag <= 0
        return;
    end
end

function y_pred = one_step_ahead_predictor_with_tf(u, y, B, A)
    % Inputs:
    %   u - Input signal
    %   y - Output signal
    %   B - Numerator coefficients of the transfer function K(z)
    %   A - Denominator coefficients of the transfer function K(z)

    % Number of data points
    N = length(u);

    % Initialize the predicted output
    y_pred = y;

    % Compute one-step-ahead prediction
    for t = max(length(A), length(B)) : N
        % Compute the output based on the transfer function
        y_pred(t) = -A(2:end) * y(t-1:-1:t-length(A)+1) + B * u(t:-1:t-length(B)+1);
    end
end