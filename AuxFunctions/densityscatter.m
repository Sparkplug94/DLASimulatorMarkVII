function densityscatter(PCx,PCy,bins,Msize)
% PCx: Column vector
% PCy: Column vector
% bins: Number of bins of the histogram
% Msize: Markersize
%
% Example:
% densityscatter(randn(1000,1),randn(1000,1),50,15)
%
% Copyright:
% José Manuel Amigo
[N C] = hist3([PCx,PCy],[bins bins]);
CX = C{1};
CY = C{2};
N2 = N;
N2(N2 == 0) = [];
Nunique = unique(N2);
colors = jet(length(Nunique));
for i = 1:length(PCx)
    if isnan(PCx(i))
        PCxnew(i,1) = NaN;
        PCynew(i,1) = NaN;
        J(i,1) = NaN;
    else
        whichoneX = find(min(abs(CX - PCx(i))) == abs(CX - PCx(i)));
        PCxnew(i,1) = CX(whichoneX(1));
        whichoneY = find(min(abs(CY - PCy(i))) == abs(CY - PCy(i)));
        PCynew(i,1) = CY(whichoneY(1));
        J(i,1) = sub2ind([bins,bins],whichoneX(1),whichoneY(1)); 
    end
end
for i = 1:bins
    for j = 1:bins
        temp = sub2ind([bins,bins],i,j);
        Jthese = find(J == temp);
        if ~isempty(Jthese)
            Ntemp = N(temp);
            Nthis = find(Nunique == Ntemp);
            plot(PCx(Jthese),PCy(Jthese),'.','color',colors(Nthis,:),'Markersize',Msize);
            hold on;
        end  
    end
end
hold off;
        
