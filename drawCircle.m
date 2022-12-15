%A circle in 3D is parameterized by six numbers: two for the orientation of its unit normal vector, one for the radius, and three for the circle center.
function [P] = drawCircle(rad, pos, n, color, LineWidth, draw)
    %https://demonstrations.wolfram.com/ParametricEquationOfACircleIn3D/
    %draws a 3D circle at position pos with radius rad, normal to the
    %circle n, and color color.
    phi = atan2(n(2),n(1)); %azimuth angle, in [-pi, pi]
    theta = atan2(sqrt(n(1)^2 + n(2)^2) ,n(3));% zenith angle, in [0,pi]    
    t = 0:pi/32:2*pi;
    %t = linspace(0, 2*pi, floor(2*pi*rad/10));
    x = pos(1)- rad*( cos(t)*sin(phi) + sin(t)*cos(theta)*cos(phi) );
    y = pos(2)+ rad*( cos(t)*cos(phi) - sin(t)*cos(theta)*sin(phi) );
    z = pos(3)+ rad*sin(t)*sin(theta);
    P = [x; y; z]';
    %plot3(x,y,z,color, 'LineWidth', LineWidth)
    %fill3(x,y,z,color, 'LineWidth', LineWidth, 'FaceAlpha', .2) % was 0.35