clear all;
close all force;
clc;
addpath('AuxFunctions')
disp('Running Script...')
tic;

%% Description
%Next version of DLA simulator, using x,y,s,xp,yp,delta
%Plane wave approximation right now - can add laser envelope later
%q = 1, c = 1
%This code needs testing and benchmarking, but initial tests look like it
%works - APF increases particle survival and bunching occurs properly (?)


%% Global Constants
c_SI = 299792458; %speed of light, m/s
c = 1; %normalized speed of light, I think this is fine for everything except initializations
me = 511e3; %mass of electron, eV
lam0 = 2e-6; %wavelength, m
T0 = 91.6e3; %injection energy, eV
[beta0, gamma0] = KE2rel(T0); %injection beta, gamma

%% Initial Beam Parameters

%beam energy
T0 = 91.6e3; %injection energy, eV
[beta0, gamma0] = KE2rel(T0); %injection beta, gamma
deltaE = 10; %energy spread (stdev), eV

%beam size
mu_x = 0; %beam centroid relative to channel, m
sigma_x = 100e-9; %standard deviation, m
mu_y = 0; %beam centroid relative to channel, m
sigma_y = 100e-9; %standard deviation, m

%beam divergence
mu_xp = 0; %beam centroid relative to channel, m
sigma_xp = 0.5e-3; %standard deviation, m
mu_yp = 0; %beam centroid relative to channel, m
sigma_yp = 0.5e-3; %standard deviation, m

%beam length (time)
sigma_tau_BEAM = 200e-15; %macrobunch length (stdev), s
sigma_s = beta0*c_SI*sigma_tau_BEAM; %macrobunch length (stdev), m

%number of particles
N = 1e5;

%synchronous phase
phi_s0 = pi/3;

%laser pulse length
%las_FWHM = 300e-15; %laser intensity FWHM, in time
las_FWHM = 0; %laser intensity FWHM, in time %SET TO ZERO FOR PLANE WAVE ILLUMINATION
sigma_tau_LAS = sqrt(2)*las_FWHM/2.355; %the laser FIELD standard deviation, also in time

%% Initialize Phase Space
    %6xN array 
    %x - gaussian, described by mu_x and sigma_x, in m
    %y - gaussian, described by mu_y and sigma_y, in m
    %s - gaussian, described by sigma_s, s = z - z_synchronous
    %xp - gaussian, described by mu_xp and sigma_xp, in radians
    %yp - gaussian, described by mu_yp and sigma_yp, in radians
    %delta - gaussian - described by mu_T0, sigma_deltaW

%init distribution
phaseSpace.distOrig = makeGaussBeam( mu_x, sigma_x,...
                            mu_y, sigma_y,...
                            mu_xp, sigma_xp,...
                            mu_yp, sigma_yp,...
                            sigma_s, T0, deltaE, N);

%init additional tracked parameters
phaseSpace.gamma_s = gamma0; %synchronous gamma
phaseSpace.phi_s = phi_s0; %laser phase/synchronous particle phase relative to laser - update laser phase during drifts
                        
%init global parameters
phaseSpace.lam0 = lam0; %wavelength
phaseSpace.gamma0 = gamma0; %injection energy
phaseSpace.N0 = N; %original particle number
phaseSpace.N = phaseSpace.N0; %Current particle number

%init tracked distribution
phaseSpace.dist = phaseSpace.distOrig;
                        


%% Test APF

%normalize number of particles to particles that initially enter
phaseSpace = remove(phaseSpace,'y',200e-9);
phaseSpace = remove(phaseSpace,'x',1e-6);
phaseSpace.N0 = phaseSpace.N;

%cosh mode kick
Lambda0 = beta0*lam0; %period length
phi_eps = 0; %phase of epsilon
eps = 1e8*exp(1i*phi_eps); %eV/m
rn = 0; %rn
theta_r = 0; %relative phase
M = 15; %#DLA cells

for i = 1:M
    phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS );
    phaseSpace = remove(phaseSpace,'y',200e-9);
    phaseSpace = remove(phaseSpace,'x',1e-6);
    Ws(i) = me*phaseSpace.gamma_s;
    Lambda(i) = gamma2beta(phaseSpace.gamma_s)*lam0;
    particles(i) = phaseSpace.N;
end
   
    %APF drift?? maybe??
    beta_s = gamma2beta(phaseSpace.gamma_s);
    APFdrift = beta_s*lam0*(1-phi_s0/pi);
    %phaseSpace = driftUpdate(phaseSpace, APFdrift);

for i = 1:M
    phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS );
    phaseSpace = remove(phaseSpace,'y',200e-9);
    phaseSpace = remove(phaseSpace,'x',1e-6);
    Ws(i+M) = me*phaseSpace.gamma_s;
    Lambda(i+M) = gamma2beta(phaseSpace.gamma_s)*lam0;
    particles(i+M) = phaseSpace.N;
end

    %APF drift?? maybe??
    beta_s = gamma2beta(phaseSpace.gamma_s);
    APFdrift2 = beta_s*lam0*(phi_s0/pi);
    %phaseSpace = driftUpdate(phaseSpace, APFdrift2);

for i = 1:M
    phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS );
    phaseSpace = remove(phaseSpace,'y',200e-9);
    phaseSpace = remove(phaseSpace,'x',1e-6);
    Ws(i+2*M) = me*phaseSpace.gamma_s;
    beta_s = gamma2beta(phaseSpace.gamma_s);
    Lambda(i+2*M) = gamma2beta(phaseSpace.gamma_s)*lam0;
    particles(i+2*M) = phaseSpace.N;
end

    %APF drift?? maybe??
    beta_s = gamma2beta(phaseSpace.gamma_s);
    APFdrift = beta_s*lam0*(1-phi_s0/pi);
    %phaseSpace = driftUpdate(phaseSpace, APFdrift);

for i = 1:M
    phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS );
    phaseSpace = remove(phaseSpace,'y',200e-9);
    phaseSpace = remove(phaseSpace,'x',1e-6);
    Ws(i+3*M) = me*phaseSpace.gamma_s;
    Lambda(i+3*M) = gamma2beta(phaseSpace.gamma_s)*lam0;
    particles(i+3*M) = phaseSpace.N;
end

    %APF drift?? maybe??
    beta_s = gamma2beta(phaseSpace.gamma_s);
    APFdrift2 = beta_s*lam0*(phi_s0/pi);
    %phaseSpace = driftUpdate(phaseSpace, APFdrift2);

for i = 1:M
    phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS );
    phaseSpace = remove(phaseSpace,'y',200e-9);
    phaseSpace = remove(phaseSpace,'x',1e-6);
    Ws(i+4*M) = me*phaseSpace.gamma_s;
    Lambda(i+4*M) = gamma2beta(phaseSpace.gamma_s)*lam0;
    particles(i+4*M) = phaseSpace.N;
end


    phaseSpace = remove(phaseSpace,'s',1e-3);
    phaseSpace = remove(phaseSpace,'xp',30e-3);
    phaseSpace = remove(phaseSpace,'yp',30e-3);

%% Complete
disp('Complete')
toc;
disp('Making Plots...')
tic;
%% Plotting

%validate phase space to remove NaN values (why do these appear??)
phaseSpace = validate(phaseSpace);

bins = 100;
figure('units','normalized','outerposition',[0 0 1 1])
plotFull(phaseSpace.dist, bins)

%survival percentage
survival = phaseSpace.N/phaseSpace.N0;
figure,plot(1:length(Ws),(Ws-me)/1e3),xlabel('DLA cell'),ylabel('Synchronous Particle Energy (keV)')
title(['Beam Energy, Survival: ' num2str(100*survival) '%'])
disp(['Survival: ' num2str(100*survival) '%'])

%DLA period length
figure,plot(1:length(Lambda),Lambda*1e6),xlabel('DLA cell'),ylabel('Synchronous Period Length (um)')
title('DLA period Length')

figure,plot(1:length(particles),100*particles/phaseSpace.N0),xlabel('DLA cell'),ylabel('Particle Survival %')
title('Particle Survival')

%%
disp('Complete')
toc;