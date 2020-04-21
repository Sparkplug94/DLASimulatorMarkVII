function [ ] = writeLattice( latticeDesign, filename )
    
    outputCSV = [filename '.csv'];    %must end in csv
    writetable( cell2table(latticeDesign), outputCSV, 'writevariablenames', false, 'quotestrings', true);

end

