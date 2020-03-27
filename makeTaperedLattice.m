%make APF Lattice
clear all
clc

%ADD Plot Betafunction 

%make lattice

M = [12, 24, 25, 26]; %length of APF cell
num = [12, 200, 200, 200]; %number of periods

count = 0;
lattice = [];
for i = 1:length(M)
    M_tmp = M(i);
    for j = 1:num(i)
        lattice = [lattice, 0]; %0 = DLA
        count = count+1;
        
        if count == M_tmp
            lattice = [lattice, 1]; %1 = Drift
            count = 0;
        end
        
    end
end
cells = length(lattice);
dlaPeriods = cells - sum(lattice);
drifts = sum(lattice);
save('lattice.mat','lattice','cells','dlaPeriods','drifts')