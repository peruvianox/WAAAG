function [ROM] = ROM (input)
% Ricky Pimentel
% 2017
% Range of motion function
Max = max(input);
Min  = min(input);

ROM = abs(Max-Min);

end