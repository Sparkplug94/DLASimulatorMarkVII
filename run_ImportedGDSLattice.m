clear all;
close all force;
clc;
addpath('AuxFunctions')
disp('Running Script...')
tic;

%% load explicit lattice from gds file
%load('PeakGrad100_AvgGrad50_ExplicitLattice.mat')
load('C:\Users\dylan\Documents\MATLAB\DLASimulatorMarkVII\GDS_readingWritingTests\run_v5_consistency.mat')
%% Global Constants
c_SI = 299792458; %speed of light, m/s
c = 1; %normalized speed of light, I think this is fine for everything except initializations
me = 511e3; %mass of electron, eV
%lam0 = 1982e-6; %wavelength, m
gamma0 = 1/sqrt(1-beta0^2);
T0 = me*(gamma0-1);

%% Initial Beam Parameters

%beam energy
deltaE = 1; %energy spread (stdev), eV

%beam size
mu_x = 0; %beam centroid relative to channel, m
sigma_x = 100e-9; %standard deviation, m
mu_y = 0; %beam centroid relative to channel, m
sigma_y = 100e-9; %standard deviation, m

%beam divergence
mu_xp = 0; %beam centroid relative to channel, m
sigma_xp = 0.1e-3; %standard deviation, m
mu_yp = 0; %beam centroid relative to channel, m
sigma_yp = 0.1e-3; %standard deviation, m

%beam length (time)
sigma_tau_BEAM = 100e-15; %macrobunch length (stdev), s
sigma_s = beta0*c_SI*sigma_tau_BEAM; %macrobunch length (stdev), m

%number of particles
N = 1e5;

%synchronous phase RELATIVE TO LASER PHASE
phi_s0 = pi/3;

disp(['Normalized Emittance: ' num2str(sigma_x*sigma_xp*1e12*beta0*gamma0) ' pm'])

%% laser params
%las_FWHM = 300e-15; %laser intensity FWHM, in time
las_FWHM = 0e-15; %laser intensity FWHM, in time %SET TO ZERO FOR PLANE WAVE ILLUMINATION
sigma_tau_LAS = sqrt(2)*las_FWHM/2.355; %the laser FIELD standard deviation, also in time
theta_r = 0; %relative dual drive phase
rn = 0; %rn coefficient
cw = 200e-9; %channel width/2, used for removing particles
ch = 500e-9; %channel height/2, used for removing particles


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
phaseSpace.phi_s = phi_s0; %laser phase/synchronous particle phase RELATIVE TO PHASE OF GRADIENT/LASER

%HOW IS PHI DEFINED?? 

%init global parameters
phaseSpace.lam0 = lam0; %wavelength
phaseSpace.gamma0 = gamma0; %injection energy

%init tracked distribution
phaseSpace.dist = phaseSpace.distOrig;
                        
%normalize number of particles to particles that initially enter
phaseSpace = remove(phaseSpace,'y',cw);
phaseSpace = remove(phaseSpace,'x',ch);
phaseSpace.N0 = N; %original particle number
phaseSpace.N = phaseSpace.N0; %Current particle number

%% APF 
for i = 1:cells
    
    if lattice_Explicit(i,2)==1;
        opt = 'Drift';
    else
        opt = 'DLA';
    end
    
    switch opt
        case 'DLA'
            Lambda(i) = lattice_Explicit(i,1);
            %DLA kick
            eps = eps_avg(i)/cos(phi_s0);
            %disp(eps)
            phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS, Lambda(i));
            phaseSpace = remove(phaseSpace,'y',cw);
            phaseSpace = remove(phaseSpace,'x',ch);
            
            %Store variables to plot
            Ts(i) = me*(phaseSpace.gamma_s-1);
            particles(i) = phaseSpace.N;
            phi_s(i) = phaseSpace.phi_s;
            %Lambda(i) = gamma2beta(phaseSpace.gamma_s)*lam0; %DLA cell length
            beta(i) = gamma2beta(phaseSpace.gamma_s);
            latticeDesign{i,1} = 'DLA';
            latticeDesign{i,2} = Lambda(i);
            latticeDesign{i,3} = beta(i);
            
        case 'Drift'

            APFDrift = lattice_Explicit(i,1);
            phaseSpace = driftUpdate(phaseSpace, APFDrift);
            
            %store variables to plot
            Ts(i) = me*(phaseSpace.gamma_s-1);
            phi_s(i) = phaseSpace.phi_s;
            particles(i) = phaseSpace.N;
            beta(i) = gamma2beta(phaseSpace.gamma_s);
            Lambda(i) = APFDrift;
            latticeDesign{i,1} = 'Drift';
            latticeDesign{i,2} = Lambda(i);
            latticeDesign{i,3} = beta(i);
    end
    
    %Remove particles that hit the wall
    phaseSpace = remove(phaseSpace,'y',cw);
    phaseSpace = remove(phaseSpace,'x',ch);
    phaseSpace = remove(phaseSpace,'delta',.01);
    savedPhaseSpace{i} = phaseSpace;
end


%% remove particles
phaseSpace = remove(phaseSpace,'s',1e-3);
phaseSpace = remove(phaseSpace,'xp',30e-3);
phaseSpace = remove(phaseSpace,'yp',30e-3);

%% Complete
disp('Complete')
toc;

%% Plotting

disp('Making Plots...')
tic;

%validate phase space to remove NaN values (why do these appear??)
phaseSpace = validate(phaseSpace);

bins = 100;
%figure('units','normalized','outerposition',[0 0 1 1])
%plotFull(phaseSpace, bins)

%x axis for plotting (plot against distance)
xAxis = 1e6*linspace(0,sum(Lambda),length(Ts));


figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1)
        plot(xAxis,(Ts)/1e3)
        xlabel('z (um)')
        ylabel('Synchronous Particle Energy (keV)')
        title('Beam Energy')
    subplot(2,3,2)
        plot(xAxis,Lambda*1e6)
        xlabel('z (um)')
        ylabel('\Lambda (um)')
        title('Period (um)')
    subplot(2,3,3)
        plot(xAxis,100*particles/phaseSpace.N0)
        ylim([0 100])
        xlabel('z (um)')
        ylabel('Particle Survival %')
        title(['Survival: ' num2str(100*phaseSpace.N/phaseSpace.N0) '%'])
    subplot(2,3,4)
        scatter(beta,100*particles/phaseSpace.N0)
        ylim([0 100])
        xlabel('Particle \beta')
        ylabel('Particle Survival %')
        title(['Particles vs. \beta'])
    subplot(2,3,5)
        plot(xAxis,phi_s/(pi)) 
        ylim([0,2])
        xlabel('z (um)')
        ylabel('synchronous phase/\pi')
        title('Synchronous Phase / \pi')
    subplot(2,3,6)
        plot(xAxis,beta) 
        ylim([0,1])
        xlabel('z (um)')
        ylabel('\beta')
        title('beta')

plotEnergyHistory(savedPhaseSpace);
title('Particle Traces')
%% Complete 
disp('Complete')
toc;

disp(['Survival: ' num2str(100*phaseSpace.N/phaseSpace.N0) '%'])
disp(['Final Beam Energy: ' num2str(Ts(end)/1e3) ' keV'])
disp(['Avg Gradient: ' num2str((Ts(end)-T0)/sum(Lambda)/1e6) ' MeV/m'])
disp(['Avg Energy Gain per Cell: ' num2str((Ts(end)-T0)/dlaPeriods) ' eV'])
%% Display Lattice Used
%latticeDesign

figure, plot(100*(Lambda'-lattice_Explicit(:,1))./Lambda')
ylabel('% Difference')
xlabel('Cell')
title('Difference between perfectly phasematched lattice and input lattice')

writeLattice(latticeDesign,'run_ImportedGDSLattice')