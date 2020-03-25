function [ phaseSpace ] = validate( phaseSpace )
    
    %removes particles with NaN attributes (not sure why these appear...)
    boolNaN = sum(isnan(phaseSpace.dist),1); %find all columns where NaN appears
    boolNaN(boolNaN>1) = 1;
    
    phaseSpace.dist(:,boolNaN == 1) = [];
    
    %removes particles with complex values (also not sure why these
    %appear...)
    boolImag = sum(abs(imag(phaseSpace.dist)),1);
    boolImag(boolImag>1) = 1;
    
    phaseSpace.dist(:,boolImag == 1) = [];
    
end

