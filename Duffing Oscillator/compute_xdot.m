function xdot = compute_xdot(x)

% % % linear oscillator
% xdot = [x(2);-x(1)];

% nonlinear oscillator
%xdot = [x(2)^3;-x(1) ];
xdot = [x(2);-0.2*x(1)-x(1)^3 ];


end