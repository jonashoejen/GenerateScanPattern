function area_now = calculateArea(R_outer, R_inner, r, resolution_mm, I, draw)

% R: the radius of the survey area circle, e.g. 50 mm for 10 cm area
% r: the radius of the target detection, e.g. 5 mm for a 1 cm target
% resoluion_mm: the amount of cells in the 

N = 2*R_outer/resolution_mm; 
N = N + 1;

canvas = zeros(N, N);
cen_canvas = floor(N/2 + 1); % center of the circle in the matrix (both x and y)
Rcanvas = floor(N/2); % radius of the circle in the matrix
Rcanvas_inner = floor(R_inner/resolution_mm);
[x, y] = meshgrid(1:N);
canvas((x - cen_canvas).^2 + (y - cen_canvas).^2 < Rcanvas^2) = 1;
canvas((x - cen_canvas).^2 + (y - cen_canvas).^2 < Rcanvas_inner^2) = 0;
area_circ = sum(canvas, 'all');

for i = 1:size(I, 1)
    canvas((x - (N*(R_outer + I(i, 1))/(2*R_outer))).^2 + (y - (N*(R_outer + I(i, 2))/(2*R_outer))).^2 < (N/(2*R_outer/r))^2) = 0;
end

if(draw)
    imshow(flip(canvas), 'InitialMagnification', 'fit')
    g = gcf;
    g.WindowState = 'maximized';
end

area_now = (1 - sum(canvas, 'all')/area_circ) * 100;

