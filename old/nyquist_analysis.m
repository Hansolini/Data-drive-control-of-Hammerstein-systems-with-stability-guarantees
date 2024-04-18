%% Find a system that can be destabilized with a gain that is too large in closed loop
syms s

% Ideal control parameters
kp = 0.05;
ki = 0;
kd = 0;

K = kp + ki/s + kd*s;

% Desired closed loop
a = 1;
Mr = 1^3/(s+a)^4; % This can be destabilized by increasing the open loop gain

% Corresponding system
% M = K*G/(1 + K*G) -> M = K*G - M*K*G -> G = M/(K - M*K)
G = expand(simplify(expand(simplify(Mr/(K - Mr*K)))));

% Convert to transfer functions
pretty(G)
[num,den] = numden(G);
G = tf(flip(double(coeffs(num,s))),flip(double(coeffs(den,s))));

[num,den] = numden(Mr);
Mr = tf(flip(double(coeffs(num,s))),flip(double(coeffs(den,s))));

[num,den] = numden(K);
K = tf(flip(double(coeffs(num,s))),flip([0, double(coeffs(den,s))]));

figure(4);
clf; grid on; hold on;
nyquist(G)
ylim([-1,1])
xlim([-2,1])

%Mr
%minreal(feedback(K*G,1))

%% Constructing a suitable nonlinearity
% Derivative
alpha = 10;
sig = @(u) 1./(1+exp(-2*u));
df = @(u)  -100*(1/4 - (sig(0)*(1 - sig(0)) - sig(u).*(1-sig(u))));
dff = @(u) df(u) - df(0) + 0.1;

% Integrate to find nonlinearity
syms u
int_df = int(dff(u));
f = matlabFunction(int_df - subs(int_df,u,0));

% Plot
u_domain = linspace(-1,1,100);
figure(1);
clf, grid on; hold on;
plot(u_domain, f(u_domain),'b')
plot(u_domain, dff(u_domain),'r')
yline(mean(df(u_domain)),'r--')
yline(1,'k-')

disp(trapz(u_domain,df(u_domain)))

%% Sample data
% Parameters
Ts = 0.1;
t = 0:Ts:500;
var = 0;%.025;

% Input data
u = 2*(1-(square(2 * pi * t / 120) + 1) / 2)'-1;

% Noise
noise = sqrt(var) * randn(size(t))';

% System
sys = ss(G); 
A = sys.A;
B = sys.B;
C = sys.C;

system_dynamics = @(time, x) A*x + B*f(interp1(t,u,time));

% Sample data
x0 = zeros(size(A,1),1);
[t,x] = ode45(system_dynamics,t,x0);
y = (C*x' + noise')';

figure(2);
clf; grid on; hold on; box on;
plot(t, u, 'k--', 'linewidth', 2)
plot(t, y, 'r', 'linewidth', 2)

xlim([0,t(end)]);
xlabel('Time');

%% VRFT
W = c2d(tf(1,[0.3, 1]),Ts);
B = [
        c2d(tf(1), Ts);         % P
        %c2d(tf(1,[1,0]), Ts);   % I
    ];

[K, theta] = VRFT1_ry(u, y, c2d(Mr,Ts), B, W, 4, []);

disp(theta)
%%

% Closed loop
% Reference
r = 1;

% System
sys = ss(G); 
A = sys.A;
B = sys.B;
C = sys.C;

% Augment the system matrices for the integral action
Ai = [A, zeros(size(A, 1), 1); -C, 0];
Bi = [B; 0];
Ci = [C, 0];

system_dynamics = @(time, x) Ai*x + [zeros(size(Ai,1)-1,1);r] + Bi*f((theta(1)*(r - Ci*x) + ki*x(end)));

% Simulation
x0 = zeros(size(Ai,1),1);
t = linspace(0,100,500);
[t_sim,x] = ode45(system_dynamics,t,x0);
y_M = Ci*x';

%% Closed loop
M = minreal(feedback(G*(kp + ki/tf('s')),1));

figure(3);
clf; grid on; hold on; box on;
plot(t, step(Mr,t), 'k', 'linewidth', 2, 'DisplayName', 'Real step response (from transfer functions)')
plot(t_sim, y_M, 'b--', 'linewidth', 2, 'DisplayName', 'Step response from state space')

xlim([0,t(end)]);
xlabel('Time')
%ylim([0,2])

legend('location', 'best')