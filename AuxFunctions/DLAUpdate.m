function [ phaseSpace ] = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_las, Lambda )
    
    %Applies appropriate kicks for a DLA with field |eps|exp(i*arg(eps))
    %then applies drift
    
    %arguments:
        %phaseSpace struct, 
        %eps (a complex number which gives the phasematched fourier harmonic and its phase,
        %theta_r, the dual drive laser phase
        %rn, the r_n coefficient of the pillars
        %Lambda, the length of the DLA cell
            %the Lambda argument is optional. IF it is given, phasematching of the
            %synchronous particle is NOT assumed and the kicks are calculated with
            %the given value of Lambda. if Lambda is NOT provided, it is calculated
            %from the synchronous beta and laser wavelength - perfect phasematching
            %is assumed
    %sigma_tau_las, the standard deviation of the laser field gaussian
        %envelope. FIELD, not intensity. The center of the gaussian is always
        %aligned with the synchronous particle
        %IF sigma_tau_las = 0, use plane wave illumination
    
    %x - gaussian, described by mu_x and sigma_x, in m
    %y - gaussian, described by mu_y and sigma_y, in m
    %s - gaussian, described by sigma_s, s = z - z_synchronous
    %xp - gaussian, described by mu_xp and sigma_xp, in radians
    %yp - gaussian, described by mu_yp and sigma_yp, in radians
    %delta - gaussian - described by mu_T0, sigma_deltaW
    
    
    %constants
    me = 511e3; %electron mass, eV;
    c = 1; %speed of light = 1 for this, because i think it works out ok.
    q = 1; %elementary charge is 1, epsilon is in units of 
    c_SI = 299792458; %speed of light, m/s, used for conversion of laser pulse length to spatial coordinates
    
    %Synchronous particle params
    phi_s = phaseSpace.phi_s;
    gamma_s = phaseSpace.gamma_s;
    beta_s = gamma2beta(gamma_s);
    pz_s = me*c*beta_s*gamma_s; %electron momentum in eV/c
    W_s = me*c^2*gamma_s; %synchronous energy in eV
    
    if ~exist('Lambda','var')
        Lambda = beta_s*phaseSpace.lam0;
    else
        pmCheck = beta_s*phaseSpace.lam0/Lambda-1;
        disp(['Phasematching Check: ' num2str(pmCheck)])
    end
    
    
    %extract particle dynamical variables from phaseSpace.dist
    x = phaseSpace.dist(1,:);
    y = phaseSpace.dist(2,:);
    s = phaseSpace.dist(3,:);
    xp = phaseSpace.dist(4,:);
    yp = phaseSpace.dist(5,:);
    delta = phaseSpace.dist(6,:);
    
    
    %find particle momentum
    gamma = gamma_s*(1 + delta); 
    pz = me*c*sqrt(gamma.^2-1); %eV/c
    
    %calculate transverse dependence of kicks
    anp = (1-rn)*(exp(1i*theta_r)+1)/2;
    anm = (1+rn)*(exp(1i*theta_r)-1)/2;
    Gamma_n = 2*pi/(gamma_s*Lambda);
    xi_n = (1-delta)/gamma_s;
    trans_y = xi_n .* ( anm.*cosh(Gamma_n.*y) + anp.*sinh(Gamma_n.*y) );
    trans_z = 1i   .* ( anp.*cosh(Gamma_n.*y) + anm.*sinh(Gamma_n.*y) );
    
    %calculate transverse dependence of synchronous kick
    trans_z_s = 1i * anp;
    
    %calculate gaussian envelope of laser beam
    if sigma_tau_las == 0
        eps_env = eps; %if sigma_tau_las is not provided, assume plane wave illumination
    else
        sigma_space_las = sigma_tau_las * c_SI;
        eps_env = eps*getGaussValue(s,0,sigma_space_las);
    end
    
    %calculate other parts of kicks
    phi = mod(phi_s - 2*pi*s/Lambda, 2*pi); %should this be reversed? maybe? I'm going to try it this way
    g = q .* eps_env .* Lambda .* exp(-1i.*phi) .* ( 1 + 1i*pi.*delta./(beta_s.^2 * gamma_s.^2) ) / (beta_s * c);
    
    %calculate other parts of synchronous kick
    g_s = q * eps * Lambda * exp(-1i*phi_s) / (beta_s * c);
    
    %calculate momentum kicks
    deltaPx = 0; %kick in px
    deltaPy = imag(g .* trans_y); %kick in py
    deltaPz = imag(g .* trans_z); %kick in pz
    deltaPz_s = imag(g_s * trans_z_s); %kick in pz for synchronous particle
    
    %calculate energy gains
    deltaW = beta_s*c*deltaPz; %energy gain
    deltaW_s = beta_s*c*deltaPz_s; %energy gain of synchronous particle
    
    %calculate kicks for tracked variables
    Dxp = deltaPx/pz_s;
    Dyp = deltaPy/pz_s;
    Ddelta = (deltaW - deltaW_s)/W_s;
    
    %calculate adiabatic damping
    A = 1 + deltaPz./pz;
    
    %update phase Space
    phaseSpace.dist(4,:) = A.*xp + Dxp;
    phaseSpace.dist(5,:) = A.*yp + Dyp;
    phaseSpace.dist(6,:) = delta + Ddelta;
    
    %update other (synchronous parameters)
    %phi_s doesn't change
    phaseSpace.gamma_s = gamma_s + deltaW_s/(me*c^2);
    
    %drift for 1 DLA cell
    phaseSpace = driftUpdate( phaseSpace, Lambda );
    
    
end

