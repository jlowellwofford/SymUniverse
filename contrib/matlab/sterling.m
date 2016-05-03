function res = stirling( val )
%STIRLING Stirling approximation (2nd order) for log(n!)
%   input: value (double)
%   output: res (double)

if val == 0
    res = 0;
else
    res = val*log(val) - val + 0.5*log(2*pi()*val);
end

end

