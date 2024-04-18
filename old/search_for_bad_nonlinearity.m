clc;
clear all;
%rng(2);

%%
syms a b s kp ki kd

Mr = 1/(s^2 + 2*s + 1);
G = 1/(s^2+a*s+b);
K = kp + s*kd;
M = simplify(G/(1+K*G));

pretty(M)

[num_M,den_M] = numden(M);

[num_Mr,den_Mr] = numden(Mr);

sol = solve(coeffs(expand(den_Mr),s) == coeffs(expand(den_M),s),[kp,kd]);

%% TODO
% - Parametrize an open loop system
% - Find a reference model that can be perfectly achieved with a PID
% controller
% - Find an expression of the control gains that lead to the perfect
% parametrization
% - Search for a nonlinearity AND system parameters that is unstable

%% Reference model
Mr = tf(1, [1, 2, 1]);

%% Closed loop simulation
% Simulation parameters
r = 1;
x0 = [0;0];
t = linspace(0,100,500);

n = 101;
u_domain = linspace(-50,50,n);

theta0 = cumtrapz(u_domain,rand(1,n));
theta0 = 0*(theta0 - theta0((n+1)/2));

theta0(end+1) = rand();
theta0(end+1) = rand();

% Constraint
% Function
Aeq =  zeros(1,n); Aeq((n+1)/2) = 1;
Beq = 0;

% Parameters
Aeq = [Aeq, zeros(size(Aeq,1),2)];

% Function
Aineq =  eye(n) - [zeros(n-1,1), eye(n-1); zeros(1,n)]; Aineq(end,end) = -1;
Bineq = -1e-1*ones(n,1);

alpha = 10;

Aineq = [Aineq; -Aineq(1:end-1,:)];
Bineq = [Bineq; alpha*ones(size(Bineq(1:end-1)))*max(diff(u_domain))];

% Parameters
Aineq = [
    Aineq, zeros(size(Aineq,1),2);
    zeros(2,n), -eye(2); % Stable open loop
    zeros(2,n), eye(2); % Positive control gains
    ];

Bineq = [
    Bineq;
    -1e-2;
    -1e-2;
    2;
    1;
    ];

% Solve
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
theta = fmincon(@(theta) J(theta,u_domain,r,t,x0),theta0,Aineq,Bineq,Aeq,Beq,[],[],[],options);

%% Results
t = linspace(0,200,500);
% Linear dynamics
A = [0,1;-theta(end-1),-theta(end)];
B = [0,1]';
C = [1,0];

% Controller
K = [1 - theta(end-1), 2 - theta(end)];

%%
% Nonlinear function
f = @(u) interp1(u_domain, theta(1:end-2), u);

% Differential equation
system_dynamics = @(time, x) A*x + B*f((r-K*x));

% Simulation
[~,x] = ode45(system_dynamics,t,x0);
y_M = (C*x');

% Differential equation
system_dynamics = @(time, x) A*x + B*(r-K*x);

% Simulation
[~,x] = ode45(system_dynamics,t,x0);
y_Mr = (C*x');

disp(theta(end-1))
disp(theta(end))

figure(1);
clf; grid on; hold on; box on;
plot(u_domain, theta0(1:end-2), 'k', 'linewidth', 2)
plot(u_domain, f(u_domain), 'b', 'linewidth', 2)
%plot(u_domain,alpha*u_domain, 'k--', 'linewidth', 2)

figure(2);
clf; grid on; hold on; box on;
plot(t, step(Mr,t), 'k', 'linewidth', 2)
plot(t, y_M, 'b', 'linewidth', 2)
plot(t, y_Mr, 'r--', 'linewidth', 2)

%% Functions
function val = J(theta,u_domain,r,t,x0)
    % Linear dynamics
    A = [0,1;-theta(end-1),-theta(end)];
    B = [0,1]';
    C = [1,0];
    
    % Controller
    K = [1 - theta(end-1), 2 - theta(end)];

    % Nonlinear function
    f = @(u) interp1(u_domain, theta(1:end-2), u);

    % Differential equation
    system_dynamics = @(time, x) A*x + B*f(K*(r-x));

    % Simulation
    [~,x] = ode45(system_dynamics,t,x0);
    y_M = (C*x');

    val = -max(abs(y_M));
end
