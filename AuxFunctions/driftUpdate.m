function [ phaseSpace ] = driftUpdate( phaseSpace, L )
    %Applies appropriate pushes for a drift length L
    
    %x - gaussian, described by mu_x and sigma_x, in m
    %y - gaussian, described by mu_y and sigma_y, in m
    %s - gaussian, described by sigma_s, s = z_synchronous - z
    %xp - gaussian, described by mu_xp and sigma_xp, in radians
    %yp - gaussian, described by mu_yp and sigma_yp, in radians
    %delta - gaussian - described by mu_T0, sigma_deltaW
    
    %extract particle dynamical variables from phaseSpace.dist
    x = phaseSpace.dist(1,:);
    y = phaseSpace.dist(2,:);
    s = phaseSpace.dist(3,:);
    xp = phaseSpace.dist(4,:);
    yp = phaseSpace.dist(5,:);
    delta = phaseSpace.dist(6,:);
    
    %Extract synchronous particle params
    phi_s = phaseSpace.phi_s;
    gamma_s = phaseSpace.gamma_s;
    beta_s = gamma2beta(gamma_s);
    
    %calculate pushes
    Dx = L*xp;
    Dy = L*yp;
    Ds = L*delta/(beta_s^2 * gamma_s^2);
    
    %update position variables
    phaseSpace.dist(1,:) = x + Dx;
    phaseSpace.dist(2,:) = y + Dy;
    phaseSpace.dist(3,:) = s + Ds;
    
    %update synchronous variables
    %no beta or gamma change
    phaseSpace.phi_s = mod(phi_s + 2*pi*L/(beta_s*phaseSpace.lam0), 2*pi); %I think this should be a plus, not a minus
    %logic - when particle drifts, it skips some of the laser oscillation -
    %skips time - skips phase - ADDS phase to oscillation??
    
end

