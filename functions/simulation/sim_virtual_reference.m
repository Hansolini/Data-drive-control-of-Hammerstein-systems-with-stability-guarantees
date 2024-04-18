function r = sim_virtual_reference(Mr,y,t)
    % Find order of numerator and denominator
    [num, den] = tfdata(Mr, 'v');
    n_num = length(num) - find(num, 1, 'first');
    n_den = length(den) - find(den, 1, 'first');
    
    % Construct time-delay filter for strictly proper reference models
    % (https://proceedings.mlr.press/v120/breschi20a/breschi20a.pdf)
    L = tf(1,[1 zeros(1,n_den - n_num)],Mr.Ts);

    % Simulate the inverse
    r = lsim(minreal(L/Mr),y,t);

    % Shift the virtual reference to compensate for time-delay
    r = circshift(r,n_num - n_den);
end
