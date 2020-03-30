function [ ] = plotEnergyHistory( savedPhaseSpace )
    
    
    figure,
    hold on
    for i = 1:length(savedPhaseSpace)
        x = savedPhaseSpace{i}.dist(1,:);
        y = savedPhaseSpace{i}.dist(2,:);
        s = savedPhaseSpace{i}.dist(3,:);
        xp = savedPhaseSpace{i}.dist(4,:);
        yp = savedPhaseSpace{i}.dist(5,:);
        delta = savedPhaseSpace{i}.dist(6,:);
    
        gamma_s(i) = savedPhaseSpace{i}.gamma_s;
        T_s(i) = 511*(gamma_s(i)-1);
        T{i} = 511*(gamma_s(i)*(1+delta)-1);
        scatter(i*ones(1,length(T{i})),T{i},'k.');
    end
        xlabel('Cell')
        ylabel('Kinetic Energy (keV)')
        plot(T_s,'r.')
        
    
end

