x = linspace(-100,100,10000);

% Logistic function
L = 1.5;
Lhat = 0.5;
a = 0.01;
k = 5;
khat = k/(Lhat/L);
f = @(x) L./(1 + exp(-k*x)) - L/2 + a*x;
fhat = @(x) Lhat./(1 + exp(-khat*x)) - Lhat/2 + a*x;
fhati = @(y) interp1(fhat(x),x,y);

figure(1);
clf; grid on; hold on;

u = linspace(-10,10,10000);

plot(u,min(max(u,-0.5),0.5),'k--','linewidth',2,'displayname','sat')
plot(u,f(u),'linewidth',2,'displayname','$f(u)$')

plot(u,fhat(u),'linewidth',2,'displayname','$\widehat{f}(u)$')
%plot(u,fhati(u),'linewidth',2,'displayname','$\widehat{f}^{-1}(u)$')

plot(u,f(fhati(u)),'linewidth',2,'displayname','$f(\widehat{f}^{-1}(u))$')

xlim([-2,2])
ylim([-2,2])

legend('location','best')

%ylim([-1,1])
