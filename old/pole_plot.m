clear all;

% System
G = tf([2,1,3], [2,1,1,2]);
K = tf([1, 2],1);

% Parameters
gain_range = linspace(0,10,100);

% Looop
figure(1);
clf; grid on; hold on; box on;
xline(0)
yline(0)

for i = 1:length(gain_range)    
    M = feedback(gain_range(i)*K*G,1);

    l = pole(M);

    plot(real(l), imag(l), 'r.', 'markersize',10)
    drawnow;
end
