function [c,ceq] = ideal_stability_constraint(theta,beta,Mr,G,deltaN)
    % K. Van Heudsen "Data-driven model reference control with asymptotically
    % guaranteed stability"
    
    % Construct controller
    K = transpose(beta)*theta;

    % Calculate constraints
    ceq = 0;
    c = norm(minreal(zpk(Mr) - minreal((1-zpk(Mr))*zpk(K)*zpk(G))),inf) - deltaN;
end
