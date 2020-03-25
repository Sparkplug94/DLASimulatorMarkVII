function [ phaseSpace ] = remove( phaseSpace, opt, param )
    
    %removes particle from phase space distribution if |opt| > param
    
    switch opt
        case 'x'
            bool = abs(phaseSpace.dist(1,:)) > param;
        case 'y'
            bool = abs(phaseSpace.dist(2,:)) > param;
        case 's'
            bool = abs(phaseSpace.dist(3,:)) > param;
        case 'xp'
            bool = abs(phaseSpace.dist(4,:)) > param;
        case 'yp'
            bool = abs(phaseSpace.dist(5,:)) > param;
        case 'delta'
            bool = abs(phaseSpace.dist(6,:)) > param;
    end
    
    
    phaseSpace.dist(:,bool == 1) = []; %remove particles
    phaseSpace.N = length(phaseSpace.dist); %update number of particles
    
    
end

