function t = conv_time(time, x, xhat, tolerance)

% Also works when replacing (x, xhat) by (Tx, z)
% Also works if tolerance is a vector

% We look for the last time such that the error is greater than "tolerance"

    t = zeros(size(tolerance));

    for i_tol = 1:length(tolerance)
    
        ind_i_tol = find(sqrt(sum((x - xhat).^2, 2)) > tolerance(i_tol), 1, 'last');
        if isempty(ind_i_tol)
            ind_i_tol = 1;
        elseif ind_i_tol < length(time)
            ind_i_tol = ind_i_tol + 1;
        end
        t(i_tol) = time(ind_i_tol);

    end

end