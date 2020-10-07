clear all
close all
clc
importname = 'C:\Users\dylan\Documents\MATLAB\DLASimulatorMarkVII\GDS_readingWritingTests\run_v5_consistency.csv';
imported = csvread(importname);
exportname = strrep(importname,'.csv','.mat'); 

%What is the wavelength?
lambda0 = 2; %in um
me = 511e3;

%units of um
centers = [];
for i = 1:length(imported)
    centers = [centers; imported(1,i) imported(2,i)];
%     pillar{i}.x = imported(1,i);
%     pillar{i}.y = imported(2,i);
%     pillar{i}.xdiam = imported(3,i);
%     pillar{i}.ydiam = imported(4,i);
end

%remove all pillars with x < 0
bool = centers(:,1) <= 0;
centers(bool == 1,:) = [];

%sort list by y value
sorted = sort(centers,1);
sorted_centers = sorted(:,2);

%find period lengths
period = diff(sorted_centers);

%find beta_s
beta_s_spiked = period/lambda0;

%find beta0
if beta_s_spiked == 0
    beta0 = beta_s_spiked(2);
else
    beta0 = beta_s_spiked(1);
end

%remove spikes due to drifts
beta_s = medfilt1(beta_s_spiked,4);
%replace first entry (medfilt screws up at first entry
beta_s(1) = beta_s_spiked(1);

%plot beta
figure
plot(beta_s_spiked);
hold on
plot(beta_s);
title('Recovered Electron Speed')

%find average gradient
gamma_s = 1./sqrt(1-beta_s.^2);
T = me*(gamma_s-1);
figure, plot(T), ylabel('Energy (eV)')
title('Recovered Electron Energy')
hold on

cellfit = linspace(1,length(T),length(T));
p = polyfit(cellfit',T,2);
Tfit = polyval(p,cellfit);
plot(Tfit)

%find energy gain per unit cell
dE = diff(Tfit); %in eV
figure
plot(diff(T))
title('Energy Gain per Cell (eV)')
hold on
plot(dE)

disp(['Average Gradient: ' num2str(sum(dE)/sum(period)) ' MeV/m'])

%cut off last entry of all exported variables so that everything now
%matches
%variables for exporting
cells = length(period)-1; %truncate last entry
isDrift = beta_s_spiked - beta_s > 1e-2;
isDrift(end) = []; %truncate last entry
drifts = sum(isDrift);
dlaPeriods = cells-drifts;
period(end) = []; %truncate last entry

figure,plot(period)
xlabel('DLA Cell')
ylabel('Periodicity or Drift Length (um)')
hold on
plot(isDrift)

lattice_Explicit = [period/1e6 isDrift];

lam0 = lambda0/1e6;

eps_avg = 1e6*dE'./period;

save(exportname,'lattice_Explicit','beta0','cells','dlaPeriods','drifts','lam0','eps_avg')

