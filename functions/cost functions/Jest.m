function v = Jest(theta,beta,Mr,G)
    % Construct controller
    K = transpose(beta)*theta;
    M = feedback(K*G,1);

    % Evaluate objective
    v = norm(minreal(zpk(Mr) - minreal((1-zpk(M))*zpk(K)*zpk(G))),2)^2;
end
