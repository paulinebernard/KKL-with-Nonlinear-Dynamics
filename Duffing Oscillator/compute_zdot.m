function zdot = compute_zdot(z,y,Lambda)



dim_z = length(Lambda);
a=5;  % strong slope far from 0 for speed
b=0.5;   % weak slope close to 0 for robustness



% KKL linear + non linear Lucas's style
% zdot = [Lambda.*(z(1:dim_z)-y);Lambda.*(z(dim_z+1:end)-atan(0.5*z(dim_z+1:end)+y))];
% 
% % KKL linear + non linear Laurent's style
%  zdot = [Lambda.*(z(1:dim_z)-y);Lambda.*(1+abs(z(dim_z+1:end)).^3).*(z(dim_z+1:end)-y)];

% zdot = [a * (Lambda.*(z(1:dim_z)-y)); Lambda.*(a*(z(dim_z+1:end)-y)+(b-a)*tanh(z(dim_z+1:end)-y))];  % linear versus nonlinear !
% factor (b-a), we want b small and a big

% KKL linear a + linear b + nonlinear 
zdot = [a*(Lambda.*(z(1:dim_z)-y)); ...
        b*(Lambda.*(z(dim_z+1:2*dim_z)-y)); ...
        Lambda.*(a*(z(2*dim_z+1:end)-y)+(b-a)*tanh(z(2*dim_z+1:end)-y))];
end