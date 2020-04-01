function [] = plotPhaseSpace( dist , opt)
    
    %dist format
    %x - gaussian, described by mu_x and sigma_x
    %y - gaussian, described by mu_y and sigma_y
    %s - 
    %xp - gaussian, described by mu_xp and sigma_xp
    %yp - gaussian, described by mu_yp and sigma_yp
    %delta - gaussian - described by mu_T0, sigma_deltaW
    
    %extract particle dynamical variables from dist
    x = dist(1,:);
    y = dist(2,:);
    s = dist(3,:);
    xp = dist(4,:);
    yp = dist(5,:);
    delta = dist(6,:);
    
    bins = 100;
    
    switch opt
        
        %plots trace space in x
        case 'x'
        %figure
        scatter(x*1e9,xp*1e3,'.')
        xlabel('x (nm)')
        ylabel('x'' (mrad)')
        
        case 'x_heat'
        densityscatter(x'*1e9,xp'*1e3,bins,1)
        xlabel('x (nm)')
        ylabel('x'' (mrad)')
        
        %plots trace space in y
        case 'y'
        %figure
        scatter(y*1e9,yp*1e3,'.')
        xlabel('y (nm)')
        ylabel('y'' (mrad)')
        
        case 'y_heat'
        densityscatter(y'*1e9,yp'*1e3,bins,1)
        xlabel('y (nm)')
        ylabel('y'' (mrad)')
        
        %plots whatever in t
        case 's'
        %figure
        scatter(s*1e6,delta,'.')
        xlabel('s (\mum)')
        ylabel('\delta (W/Ws-1)')
        
        case 's_heat'
        densityscatter(s'*1e6,delta',100,1)
        xlabel('s (\mum)')
        ylabel('\delta (W/Ws-1)')
        
        case '3D'
        P = [x'*1e9 y'*1e9 s'*1e6];
        k = boundary(P);
        scatter3(P(:,1),P(:,2),P(:,3),'.')
        hold on
        trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','red','FaceAlpha',0.1)
        xlabel('x (nm)')
        ylabel('y (nm)')
        zlabel('s (\mum)')
        
        case 'foc'
        %figure
        scatter(s*1e6,yp*1e3,'.')
        xlabel('s (\mum)')
        ylabel('y'' (mrad)')
        
        case 'xy'
        scatter(x*1e6,y*1e6,'.')
        xlabel('x (\mum)')
        ylabel('y (\mum)')  
            
    end
    set(gca,'FontSize',14)
    
    
end

