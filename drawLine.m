
function [xyz] = drawLine(p1, p2, draw, color)
    %eval(['denseQ' num2str(d(i)) '= dense;']);
    u = (p2-p1)/norm(p2-p1);   % unit vector, p1 to p2
    d_ = norm(p2 - p1);
    d = linspace(0, d_, 1000)';
    %d = (0:.1:d)';            % displacement from p1, along u
    xyz = p1 + d.*u; 
    if(draw)
        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'-', 'color', color, 'LineWidth', 1);
    end
end