function [] = longHist( dist, bins )
    %Creates histogram of longitudinal coordinate
    s = dist(3,:);
    [N,bins] = hist(s,bins);
    plot(bins*1e6,N/max(N));
    xlabel('s (um)')
    ylabel('Number of Particles (Normalized)')
    set(gca,'FontSize',14)
    
end

