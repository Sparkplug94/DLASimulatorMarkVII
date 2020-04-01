%Demonstration file

%I want to simulate two DLA stages, each with a 100 MeV/m acceleration
%gradient, in the cosh mode. I want a drift length between them of length
%14.8 um. 

    %x - in meters
    %y - in meters
    %s - s = z - z_synchronous, in meters
    %xp - angle in x coordinate, radians
    %yp - angle in y coordinate, radians
    %delta - energy deviation (Energy - Energy of Ref Particle)/Energy of
    %Ref particle


%Standard header
clear all;
close all force;
clc;
addpath('AuxFunctions')
disp('Running Script...')

%% Constants
c_SI = 299792458; %speed of light, m/s
me = 511e3; %electron mass, eV

%% Laser Parameters
lam0 = 2e-6; %wavelength, m
las_FWHM = 0e-15; %laser intensity FWHM, s, %SET TO ZERO FOR PLANE WAVE ILLUMINATION
sigma_tau_LAS = sqrt(2)*las_FWHM/2.355; %the laser FIELD standard deviation, also in time

%% Electron Beam Parameters

%beam energy
T0 = 91.6e3; %injection energy, eV
[beta0, gamma0] = KE2rel(T0); %injection beta, gamma
deltaE = 1; %beam initial energy spread (stdev), eV

%beam size
mu_x = 0; %beam centroid relative to channel, m
sigma_x = 200e-9; %standard deviation, m
mu_y = 0; %beam centroid relative to channel, m
sigma_y = 200e-9; %standard deviation, m

%beam divergence
mu_xp = 0; %beam centroid relative to channel, m
sigma_xp = 0.5e-3; %standard deviation, m
mu_yp = 0; %beam centroid relative to channel, m
sigma_yp = 0.5e-3; %standard deviation, m

%beam length (time)
sigma_tau_BEAM = 10e-15; %macrobunch length (stdev), s
sigma_s = beta0*c_SI*sigma_tau_BEAM; %macrobunch length (stdev), m

%number of particles
N = 1e4;

%Reference particle phase RELATIVE TO LASER PHASE
phi_s0 = 0; %phi_s0 = 0 indicates that the reference particle is on-crest - i.e. experiences the maximum acceleration gradient


%% Initialize Gaussian Beam

%init tracked distribution
phaseSpace.dist = makeGaussBeam( mu_x, sigma_x,...
                            mu_y, sigma_y,...
                            mu_xp, sigma_xp,...
                            mu_yp, sigma_yp,...
                            sigma_s, T0, deltaE, N);

%save original distribution (optional)
phaseSpace.distOrig = phaseSpace.dist;
                        
%init additional tracked parameters
phaseSpace.gamma_s = gamma0; %synchronous gamma
phaseSpace.phi_s = phi_s0; %laser phase/synchronous particle phase RELATIVE TO PHASE OF GRADIENT/LASER

%init global parameters
phaseSpace.lam0 = lam0; %wavelength
phaseSpace.gamma0 = gamma0; %injection energy


%% Plot initial Gaussian Beam

%let's say I want to view the phase space 
figure('units','normalized','outerposition',[0 0 1 1])

subplot(4,3,1)
plotPhaseSpace(phaseSpace.dist, 'y')
title('y phase space, initial')

subplot(4,3,2)
plotPhaseSpace(phaseSpace.dist, 's')
title('s phase space, initial')

%I also want to plot the longitudinal particle density (i.e. bunching)
subplot(4,3,3)
bins = 100;
longHist(phaseSpace.dist,bins)
title('longitudinal density, initial')

%% 1st DLA segment

n1 = 20; %20 DLA periods
eps_1 = 100e6; %Peak Acceleration Gradient, in MeV/m 
theta_1 = 0; %relative laser phase
rn = 0;

%run n1 DLA cells
for i = 1:n1
    phaseSpace = DLAUpdate( phaseSpace, eps_1, theta_1, rn, sigma_tau_LAS);
    
    %myLambda = 1e-6;
    %phaseSpace = DLAUpdate( phaseSpace, eps_1, theta_1, rn, sigma_tau_LAS, myLambda);
        %Optionally, you may add an argument to the end of the DLAUpdate
        %function specifying the specific periodicity. If you do not, the code
        %will assume it is phasematched with the reference particle, and derive
        %period length from phaseSpace.lam0*beta_reference
end


%% Remove Outlier particles

%It is usually smart to remove outlier particles. Every now and then a
%particle will end up waaaay outside where it should be.

ymax = 250e-9; %for example, let's say my channel width is 500 nm
phaseSpace = remove(phaseSpace,'y',ymax); %Remove particles with abs(y) > ymax

deltamax = 0.05; %Remove particles with >5% energy deviation from reference particle
phaseSpace = remove(phaseSpace,'delta',deltamax); %Remove particles with abs(delta) > deltamax

%% Plot new phase space

subplot(4,3,4)
plotPhaseSpace(phaseSpace.dist, 'y')
title('y phase space, after DLA 1')

subplot(4,3,5)
plotPhaseSpace(phaseSpace.dist, 's')
title('s phase space, after DLA 1')

subplot(4,3,6)
bins = 100;
longHist(phaseSpace.dist,bins)
title('longitudinal density, after DLA 1')

%% Drift Segment

%drift for distance L
L = 13.7e-6;
phaseSpace = driftUpdate(phaseSpace, L);

%% Plot new phase space

subplot(4,3,7)
plotPhaseSpace(phaseSpace.dist, 'y')
title('y phase space, after drift')

subplot(4,3,8)
plotPhaseSpace(phaseSpace.dist, 's')
title('s phase space, after drift')

subplot(4,3,9)
bins = 100;
longHist(phaseSpace.dist,bins)
title('longitudinal density, after drift')

%% 2nd DLA segment.

n2 = 20; %10 DLA periods
eps_2 = 100e6; %Peak Acceleration Gradient, in MeV/m 
theta_2 = 0; %relative laser phase
rn = 0;

%run n1 DLA cells
for i = 1:n1
    phaseSpace = DLAUpdate( phaseSpace, eps_2, theta_2, rn, sigma_tau_LAS);

end

%Remove Outlier Particles
ymax = 250e-9; %for example, let's say my channel width is 500 nm
phaseSpace = remove(phaseSpace,'y',ymax); %Remove particles with abs(y) > ymax

deltamax = 0.05; %Remove particles with >5% energy deviation from reference particle
phaseSpace = remove(phaseSpace,'delta',deltamax); %Remove particles with abs(delta) > deltamax

% ypmax = 1; %for example, let's say my channel width is 500 nm
% phaseSpace = remove(phaseSpace,'yp',ypmax); %Remove particles with abs(y) > ymax

%% Plot new phase space

subplot(4,3,10)
plotPhaseSpace(phaseSpace.dist, 'y')
title('y phase space, after DLA 2')

subplot(4,3,11)
plotPhaseSpace(phaseSpace.dist, 's')
title('s phase space, after DLA 2')

subplot(4,3,12)
bins = 100;
longHist(phaseSpace.dist,bins)
title('longitudinal density, after DLA 2')


%% Let's make plots of the actual beam energy

%the particle distribution is store in phaseSpace.dist
%the coordinates are as follows:
    %phaseSpace.dist(1,:) = x - in meters
    %phaseSpace.dist(2,:) = y - in meters
    %phaseSpace.dist(3,:) = s - s = z - z_synchronous, in meters
    %phaseSpace.dist(4,:) = xp - angle in x coordinate, radians
    %phaseSpace.dist(5,:) = yp - angle in y coordinate, radians
    %phaseSpace.dist(6,:) = delta - energy deviation

%% Final Energy
%extract parameters from phaseSpace object
gamma_s = phaseSpace.gamma_s; %extract energy of synchronous particle
delta = phaseSpace.dist(6,:); %extract delta of particles

%convert to more interpretable units
W_s = me*gamma_s; %Find energy of synchronous particle
W = me*gamma_s*(1+delta);%find energy of all particles in phaseSpace

%Convert to kinetic energy
T_s = W_s - me;
T = W - me;

%% Initial Energy

%do the same for the original distribution
gamma_s0 = gamma0; %extract energy of synchronous particle
delta0 = phaseSpace.distOrig(6,:); %extract delta of particles

%convert to more interpretable units
W_s0 = me*gamma_s0; %Find energy of synchronous particle
W0 = me*gamma_s0*(1+delta0);%find energy of all particles in phaseSpace


%% Plots

%plot
figure
hist(T0/1e3,100)
title('Initial Energy')

%plot
figure
hist(T/1e3,100)
title('Final Energy')
