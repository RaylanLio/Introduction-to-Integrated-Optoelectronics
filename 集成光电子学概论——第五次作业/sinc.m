function y = sinc(x)
%sinc Summary of this function goes here
%   Detailed explanation goes here
y = zeros(1,length(x));
for i = 1:length(x)
    if x(i) == 0
        y(i) = 1;
    else
        y(i) = sin(x(i))./x(i);
    end
end