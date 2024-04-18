function [c,ceq] = estimated_stability_constraint(theta,beta,Mr,u,y,t,deltaN,bDelta,G)
    % K. Van Heudsen "Data-driven model reference control with asymptotically
    % guaranteed stability"
    
    % Construct controller
    K = transpose(beta)*theta;

    delta_1 = bDelta;
    delta_2 = 1;

    % Estimate virtual reference & error signal
    epsilon_1 = lsim(Mr,u,t) - delta_2/delta_1*lsim(K*(1-Mr),y,t);
    epsilon_2 = lsim(Mr,u,t) - delta_1/delta_2*lsim(K*(1-Mr),y,t);

    % Estimate H-inf norm
    l2 = 200;
    Pr = cpsd(u, u, rectwin(2*l2));
    Pre_1 = cpsd(epsilon_1, u, rectwin(2*l2));
    Pre_2 = cpsd(epsilon_2, u, rectwin(2*l2));

    delta_1 = max(abs(Pre_1./Pr));
    delta_2 = max(abs(Pre_2./Pr));

    % Calculate constraints
    ceq = 0;
    c = max(delta_1, delta_2) - deltaN;
    
end
