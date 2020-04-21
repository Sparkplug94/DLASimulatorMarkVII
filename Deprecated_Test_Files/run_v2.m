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

%peak gradient
eps = 200e6; %can be complex number

%synchronous phase
phi_s0 = pi/3;

%% laser params
%las_FWHM = 300e-15; %laser intensity FWHM, in time
las_FWHM = 0; %laser intensity FWHM, in time %SET TO ZERO FOR PLANE WAVE ILLUMINATION
sigma_tau_LAS = sqrt(2)*las_FWHM/2.355; %the laser FIELD standard deviation, also in time
theta_r = 0; %relative dual drive phase
rn = 0; %rn coefficient
cw = 200e-9; %channel width/2, used for removing particles
ch = 1e-6; %channel height/2, used for removing particles


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

%init tracked distribution
phaseSpace.dist = phaseSpace.distOrig;
                        

%normalize number of particles to particles that initially enter
phaseSpace = remove(phaseSpace,'y',cw);
phaseSpace = remove(phaseSpace,'x',ch);
phaseSpace.N0 = N; %original particle number
phaseSpace.N = phaseSpace.N0; %Current particle number

%% APF Section 1

%Synchronous Phase Corrections are applied every time an APF drift is.

dlaPeriods = 250; %number of DLA cells

M = 5; %apf cell length


APF = 'on';

switch APF
    case 'on'
        flip = true;
        disp('APF ON')
        drifts = floor(dlaPeriods/M);
    case 'reverse'
        flip = false;
        disp('APF ON & Reversed')
        drifts = floor(dlaPeriods/M);
    case 'off'
        flip = NaN;
        disp('APF Off')
        drifts = 0;
end

cells = dlaPeriods+drifts;


for i = 1:cells
    
    if mod(i,M)==0;
        opt = 'Drift';
    else
        opt = 'DLA';
    end
    
    switch opt
        case 'DLA'
            
            disp('DLA')
            %DLA kick
            phaseSpace = DLAUpdate( phaseSpace, eps, theta_r, rn, sigma_tau_LAS);
            phaseSpace = remove(phaseSpace,'y',cw);
            phaseSpace = remove(phaseSpace,'x',ch);
            
            %Store variables to plot
            Ts(i) = me*(phaseSpace.gamma_s-1);
            particles(i) = phaseSpace.N;
            phi_s(i) = phaseSpace.phi_s;
            Lambda(i) = gamma2beta(phaseSpace.gamma_s)*lam0; %DLA cell length
            beta(i) = gamma2beta(phaseSpace.gamma_s);
        case 'Drift'
            disp('Drift')
            %previous Lambda
            if i-M >= 1
                prevLam = Lambda(i-M);
            elseif i-M == 0
                prevLam = beta0*lam0;
            end
            %current Lambda
            curLam = Lambda(i-1);
            %chirp 
            extraLam = curLam-prevLam;
            %phase correction drift
            phaseCorrDrift = cot(phaseSpace.phi_s)*extraLam/(2*pi);
            %phaseCorrDrift = 0;
            %get current speed
            beta_s = gamma2beta(phaseSpace.gamma_s);
           
            %apply one of the two APF Drifts plus the phase correction
            %drift
            if flip == true
                APFDrift = beta_s*lam0*(1-phi_s0/pi) + phaseCorrDrift;
                phaseSpace = driftUpdate(phaseSpace, APFDrift);
                Lambda(i) = APFDrift;
                flip = false;
            elseif flip == false
                APFDrift2 = beta_s*lam0*(phi_s0/pi) + phaseCorrDrift;
                phaseSpace = driftUpdate(phaseSpace, APFDrift2);
                Lambda(i) = APFDrift2;
                flip = true;
            else
            end
            
            %disp('APF Drift')
            
            %store variables to plot
            Ts(i) = me*(phaseSpace.gamma_s-1);
            phi_s(i) = phaseSpace.phi_s;
            particles(i) = phaseSpace.N;
            beta(i) = gamma2beta(phaseSpace.gamma_s);
            
    end
    
    %Remove particles that hit the wall
    phaseSpace = remove(phaseSpace,'y',cw);
    phaseSpace = remove(phaseSpace,'x',ch);

end


%% remove particles
phaseSpace = remove(phaseSpace,'s',1e-3);
phaseSpace = remove(phaseSpace,'xp',30e-3);
phaseSpace = remove(phaseSpace,'yp',30e-3);

%% Complete
disp('Complete')
toc;

%% Plotting

%validate phase space to remove NaN values (why do these appear??)
phaseSpace = validate(phaseSpace);

bins = 100;
%figure('units','normalized','outerposition',[0 0 1 1])
%plotFull(phaseSpace.dist, bins)

%x axis for plotting (plot against distance)
xAxis = 1e6*linspace(0,sum(Lambda),length(Ts));

disp('Making Plots...')
tic;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
plot(xAxis,(Ts)/1e3)
xlabel('z (um)')
ylabel('Synchronous Particle Energy (keV)')
title('Beam Energy')

subplot(2,3,3)
plot(xAxis,100*particles/phaseSpace.N0)
xlabel('z (um)')
ylabel('Particle Survival %')
title(['Survival: ' num2str(100*phaseSpace.N/phaseSpace.N0) '%'])

subplot(2,3,5)
plot(xAxis,phi_s/(pi)) 
ylim([0,2])
xlabel('z (um)')
ylabel('synchronous phase/\pi')
title('Synchronous Phase / \pi')


subplot(2,3,2)
plot(xAxis,Lambda*1e6)
xlabel('z (um)')
ylabel('\Lambda (um)')
title('Period (um)')

subplot(2,3,6)
plot(xAxis,beta) 
ylim([0,1])
xlabel('z (um)')
ylabel('\beta')
title('beta')

%% Complete 
disp('Complete')
toc;

disp(['Survival: ' num2str(100*phaseSpace.N/phaseSpace.N0) '%'])
disp(['Final Beam Energy: ' num2str(Ts(end)/1e3) ' keV'])
disp(['Avg Gradient: ' num2str((Ts(end)-T0)/sum(Lambda)/1e6) ' MeV/m'])
disp(['Avg Energy Gain per Cell: ' num2str((Ts(end)-T0)/dlaPeriods) ' eV'])
%%