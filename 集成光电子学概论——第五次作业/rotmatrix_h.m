function R = rotmatrix_h(t)
% Generate the rotating matrix
%   Detailed explanation goes here
R=[cosh(t) -sinh(t);-sinh(t) cosh(t)];
end