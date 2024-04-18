function [y,t] = step_closed_loop_hammerstein(G,K,f,t)
    % System
    system = ss(G);
    As = system.A;
    Bs = system.B;
    Cs = system.C;

    if abs(system.D) > 1e-5
        disp(system.D)
        disp('There is direct feedthrough in the system! Hammerstein simulation will be inaccurate. Aborted.')
        y = -ones(size(t));
        return;
    end

    ns = size(As,1);

    % Controller
    controller = ss(K);
    Ac = controller.A;
    Bc = controller.B;
    Cc = controller.C;
    Dc = controller.D;

    nc = size(Ac,1);
    
    % Reference (step)
    r = @(t) 1;

    % Closed loop system
    A = [As, zeros(ns,nc); -Bc*Cs, Ac];
    B = blkdiag(Bs, Bc);
    C = [Cs, zeros(size(Cc))];

    phi = @(t,x) A*x + B*[f(Cc*x(ns+1:end) - Dc*Cs*x(1:ns) + Dc*r(t)); r(t)];

    % Simulation
    x = zeros(length(t),ns+nc);
    for i = 2:length(t)
        x(i,:) = phi(t(i),x(i-1,:)')';
    end

    y = x*C';
end