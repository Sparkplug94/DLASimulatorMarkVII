function [ ] = plotFull( phaseSpace, bins )
    
    
    subplot(2,3,1)
    plotPhaseSpace(phaseSpace.dist,'x')
    subplot(2,3,2)
    plotPhaseSpace(phaseSpace.dist,'y')
    subplot(2,3,3)
    plotPhaseSpace(phaseSpace.dist,'s')
    subplot(2,3,5)
    plotPhaseSpace(phaseSpace.dist,'3D')
    subplot(2,3,6)
    longHist(phaseSpace.dist, bins)
    subplot(2,3,4)
    plotPhaseSpace(phaseSpace.dist,'foc')
    
end

