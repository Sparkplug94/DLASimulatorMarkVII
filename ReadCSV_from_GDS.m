clear all
close all
clc
imported = csvread('PeakGrad100_dims.csv');
%units of um
centers = [];
for i = 1:length(imported)
    centers = [centers; imported(1,i) imported(2,i)];
    pillar{i}.x = imported(1,i);
    pillar{i}.y = imported(1,i);
    pillar{i}.xdiam = imported(3,i);
    pillar{i}.ydiam = imported(4,i);
end

%remove all pillars with x < 0
bool = centers(:,1) < 0;
centers(bool == 1,:) = [];

%sort list by y value
sorted = sort(centers,1);
sorted_centers = sorted(:,2);

%find period lengths
period = diff(sorted_centers);

%find beta_s
lambda0 = 1.982; %in um
beta_s_spiked = period/lambda0;

%find beta0
beta0 = beta_s_spiked(1);

%remove spikes due to drifts
beta_s = medfilt1(beta_s_spiked,7);

%find average gradient
gamma_s = 1./sqrt(1-beta_s.^2);
deltaE = 511e3*(gamma_s(end)-gamma_s(1));
dE = 511e3*diff(gamma_s);
dE = smooth(medfilt1([dE(1); dE]));

g_avg = deltaE/sum(period); %eV/um = MeV/m;

disp(['Average Gradient: ' num2str(g_avg) ' MeV/m'])



cells = length(period);
isDrift = beta_s_spiked - beta_s > 1e-2;
drifts = sum(isDrift);
dlaPeriods = cells-drifts;


figure,plot(period)
xlabel('DLA Cell')
ylabel('Periodicity or Drift Length (um)')
hold on
plot(isDrift)

lattice_Explicit = [period/1e6 isDrift];

lam0 = lambda0/1e6;

eps_avg = 1e6*dE./period;

save('PeakGrad100_AvgGrad50_ExplicitLattice.mat','lattice_Explicit','beta0','cells','dlaPeriods','drifts','lam0','eps_avg')

