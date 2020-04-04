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
sigma_yp = 0e-3; %standard deviation, m

%beam length (time)
sigma_tau_BEAM = 10e-15; %macrobunch length (stdev), s
sigma_s = beta0*c_SI*sigma_tau_BEAM; %macrobunch length (stdev), m

%number of particles
N = 1e5;

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


%% 1st DLA segment

n1 = 100; %20 DLA periods
eps_1 = 100e6; %Peak Acceleration Gradient, in MeV/m 
theta_1 = 0; %relative laser phase
rn = 0;

%run n1 DLA cells
for i = 1:n1
    phaseSpace = DLAUpdate( phaseSpace, eps_1, theta_1, rn, sigma_tau_LAS);
   
end
