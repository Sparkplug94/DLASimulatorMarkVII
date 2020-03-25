function [ beta ] = gamma2beta( gamma )
    %gives gamma from beta
    
    beta = sqrt(1-1./gamma.^2);
    
    
end

