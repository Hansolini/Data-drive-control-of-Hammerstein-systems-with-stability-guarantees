function v = Jvrft(theta,beta,Mr,u,y,t)
    % Construct virtual reference & virtual error
    r = sim_virtual_reference(Mr,y,t);
    e = r - y;
    
    % Simulate controller output
    x = lsim(beta,e,t);

    % Compute cost
    v = norm(x*theta - u)^2;
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
    for i = max(length(A), length(B)):N
        y_pred(i) = -A(2:end) * y(i-1:-1:i-length(A)+1) + B * u(i:-1:i-length(B)+1);
    end
end
