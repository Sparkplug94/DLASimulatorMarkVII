function [ f ] = getGaussValue(x, mu, s)

    %extract value from gaussian distribution where peak value is 1
    
    %gets value of gaussian distribution at position x relative to mean mu,
    %with standard deviation s
    
    arg = ((x - mu)/(sqrt(2)*s));
    f = exp(-(arg.^2));

end
