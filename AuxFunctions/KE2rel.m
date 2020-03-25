function [ beta, gamma ] = KE2rel( T )
    % gives gamma and beta
    % takes argument of kinetic energy in eV
    me = 511e3;
    gamma = 1+T/me;
    beta = sqrt(1-1/gamma^2);
    
end

