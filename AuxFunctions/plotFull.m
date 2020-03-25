function [ output_args ] = plotFull( dist, bins )
    
    
    subplot(2,3,1)
    plotPhaseSpace(dist,'x')
    subplot(2,3,2)
    plotPhaseSpace(dist,'y')
    subplot(2,3,3)
    plotPhaseSpace(dist,'s')
    subplot(2,3,5)
    plotPhaseSpace(dist,'3D')
    subplot(2,3,6)
    longHist(dist, bins)
    subplot(2,3,4)
    plotPhaseSpace(dist,'foc')
    
end

