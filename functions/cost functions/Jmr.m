function v = Jmr(theta,beta,Mr,G)
    % Construct close loop transfer function
    K = transpose(beta)*theta;
    M = feedback(K*G,1);

    % Evaluate objective (zpk and minreal used to avoid numerical issues)
    v = norm(minreal(zpk(Mr) - zpk(M)),2)^2;
end
