[O_sparse, unit_vec_sparse, O_dense, unit_vec_dense] = analyseData();

%%
% If 70 mm is the mean, then we can 
% either go a cm up or down and see the results. 
% Positive is further away, negative is closer

%--------------------------

%Nref = sparseN70;
%Nref = sparseN(13, 1:3);
%D = double(sparseD70);
D = 105;
dist2plane = norm([0, 0, D] - O_sparse);
%on_plane = O_sparse + dist2plane.*Nref;

a = 25.5;
v = 18/180*pi;
% distance to target plane (tp) along the boresight (z-axis)
h_tp = 56; % mm

% y-coordinate
a_medium = h_tp*sin(v) + 1.60033;
h_y = a_medium/cos(v);
p1 = [0, -h_y, 0];
h_small = 1.60033/sin(v);
p2 = [0, 0, h_tp + h_small];

% Unit vec from p2 to p1 % end - start
uv_xray = (p1 - p2)./norm(p1 - p2);
h_xray = 25.5 + 1.60033/tan(v);
o_xray = h_xray*sin(v);
a_xray = h_xray*cos(v);

% x-ray location
p_xray = [0, -o_xray, p2(3) - a_xray];
h_gnd_xray = D - p2(3);
len_gnd_xray = h_gnd_xray*cos(v);

pivot_point_dist = 190 + h_small; % - norm(p_xray - p2); %%len_gnd_xray;
pivot_point = uv_xray*pivot_point_dist + p2;

%pos = -(h_xray + len_gnd_xray).*uv_xray + p_xray;
plane_dist = [+10, -10]; 
pos = -(h_xray + len_gnd_xray + plane_dist(1)).*uv_xray + p_xray;
pos2 = -(h_xray + len_gnd_xray + plane_dist(2)).*uv_xray + p_xray;
pos3 = -(h_xray + len_gnd_xray).*uv_xray + p_xray; %  + 1

Nref = -uv_xray;



%% 

% Draw circle
rad = 50; % 100mm = 10cm
r = 5; % radius of detectable object
%CP = drawCircle(rad,pos,Nref,'b', 2);
CP4 = drawCircle(rad,pos3,Nref,'r', 2);
CP = drawCircle(rad,pos,Nref,'r', 2);
CP2 = drawCircle(rad,pos2,Nref,'g', 2);
CP3 = [CP; CP2];
k = convhull(CP3(1:end-1, 1), CP3(1:end-1, 2));
CP3 = unique(CP3(k, 1:2), 'rows');
area_total = polyarea(CP4(:, 1), CP4(:, 2));

%Step 1: Find the unweighted mean of the vertices:
cx = mean(CP(:, 1));
cy = mean(CP(:, 2));

%Step 2: Find the angles:
angles = atan2(CP(:, 2) - cy, CP(:, 1) - cx);

%Step 3: Find the correct sorted order:
[~, order] = sort(angles);

%Step 4: Reorder the coordinates:
CP3 = CP(order, 1:2);
%%



%function analyseData()

%clearvars
close all

folder = 'C:\Users\jonas\Desktop\30220 Synthesis in Earth and Space Physics Fall 22\Jonas_scan_pattern_01\';
dir = {'pic1', 'pic2', 'pic3', 'pic4', 'pic5', 'pic6', 'pic7', 'pic8', 'pic9', 'pic10', 'pic11', 'pic12', 'pic13', 'pic14', 'pic15', 'pic16', 'pic17', 'pic18', 'pic19', 'pic20', 'pic21', 'pic22', 'pic23', 'pic24', 'pic25', 'pic26', 'pic27', 'pic28', 'pic29', 'pic30'};
showImg = 0;
fig = 0;
fig_add = 0;
denseQ = [];
sparseQ = [];
denseN = [];
sparseN = [];
denseD = [];
sparseD = [];

rot = [...
    5.00, -10.0, 0; ... % 1
    -6.8, -10.0, 0; ... % 2
    0.40, -8.10, 0; ... % 3
    8.70, -7.00, 0; ... % 4
    1.80, -6.70, 0; ... % 5
    -1.8, -5.20, 0; ... % 6
    -6.6, -4.70, 0; ... % 7
    -2.5, -2.80, 0; ... % 8
    -5.5, -2.80, 0; ... % 9
    -3.9, -2.00, 0; ... % 10
    -1.3, -1.20, 0; ... % 11
    10.7, -0.10, 0; ... % 12
    0.000, 0.00, 0; ... % 13
    -2.20, 0.00, 0; ... % 14
    -8.70, 0.80, 0; ... % 15
    -8.50, 1.40, 0; ... % 16
    -7.20, 2.80, 0; ... % 17
    7.600, 3.30, 0; ... % 18
    3.100, 3.60, 0; ... % 19
    4.700, 4.90, 0; ... % 20
    -1.80, 5.10, 0; ... % 21
    -4.10, 5.10, 0; ... % 22
    8.400, 5.10, 0; ... % 23
    9.300, 6.00, 0; ... % 24
    -5.70, 7.40, 0; ... % 25
    -8.20, 8.40, 0; ... % 26
    2.70, 8.900, 0; ... % 27
    -0.80, 9.40, 0; ... % 28
    -7.20, 10.0, 0; ... % 29
    6.50, 10.90, 0];    % 30

%rot = rot.*[1, -1, 1];

for i = 1:length(dir)
    [fig_add, im1, im2] = readFiles(strcat(folder, dir{i}), fig, showImg);
    dense = [im1.sliData.Q(im1.sliData.sliID ~= 255, :)/1e3, im1.sliData.sliID(im1.sliData.sliID ~= 255, :)];
    sparse = [im2.sliData.Q(im2.sliData.sliID ~= 255, :)/1e3, im2.sliData.sliID(im2.sliData.sliID ~= 255, :)];
    
    N_d = im1.sliHeader.planeN;
    D_d = im1.sliHeader.planeD/1e3;
    N_s = im2.sliHeader.planeN;
    D_s = im2.sliHeader.planeD/1e3;
    
    denseQ = [denseQ; [dense, i*ones(size(dense, 1), 1)]];
    sparseQ = [sparseQ; [sparse, i*ones(size(sparse, 1), 1)]];
    denseN = [denseN; [N_d, i*ones(size(N_d, 1), 1)]];
    sparseN = [sparseN; [N_s, i*ones(size(N_s, 1), 1)]];
    denseD = [denseD; [D_d, i*ones(size(D_d, 1), 1)]];
    sparseD = [sparseD; [D_s, i*ones(size(D_s, 1), 1)]];
    
    %eval(['denseQ' num2str(i) '= dense;']);
    %eval(['sparseQ' num2str(i) '= sparse;']);
    %eval(['denseN' num2str(i) '= N_d;']);
    %eval(['sparseN' num2str(i) '= N_s;']);
    %eval(['denseD' num2str(i) '= D_d;']);
    %eval(['sparseD' num2str(i) '= D_s;']);
end
R = [];
N = [0, 0, 1];
e = cross(N, -uv_xray)/norm(cross(N, -uv_xray)); % normalised rotation axis; 
%ve = subspace(N, Nref'); % rotation angle
ve = atan2(norm(cross(-uv_xray,N)),dot(-uv_xray,N));
q = [cos(ve/2), e*sin(ve/2)]; % quaternion

Rc = [1 - 2*(q(3)^2 + q(4)^2), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(1)*q(3) + q(2)*q(4));
     2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(2)^2 + q(4)^2), 2*(q(3)*q(4) - q(1)*q(2));
     2*(q(2)*q(4) - q(1)*q(3)), 2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2)];

Rc = [Rc, zeros(3, 1); ones(1, 4).*[0, 0, 0, 1]]; 
rot_denseQ = [];
denseQ_ = [];
rdQ = [];
rot_sparseQ = [];
rsQ = [];
sparseQ_ = [];
D_sparse = [];
D_sparse2 = [];
I = [];
I2 = [];
pos_j = [];
rot_xy_ = [];
rot_xy = [];

for j = 1:30%29%1:30
    
    for i = j
        I_sparse = rotateSLI(Nref, uv_xray, O_sparse, unit_vec_sparse, pivot_point, dist2plane, pos3, p_xray, rot(i, :), [0.5, 0.5, 0.5], 0);
        I_dense = rotateSLI(Nref, uv_xray, O_dense, unit_vec_dense, pivot_point, dist2plane, pos3, p_xray, rot(i, :), [0.5, 0.5, 0.5], 0);
        I = [I; [[I_sparse, j*ones(size(I_sparse, 1), 1)]; [I_dense, j*ones(size(I_dense, 1), 1)]]];
    end
    
    for i = j
        I_sparse = rotateSLI_mistake(Nref, uv_xray, O_sparse, unit_vec_sparse, pivot_point, dist2plane, pos3, p_xray, rot(i, :), [0.5, 0.5, 0.5], 0);
        I_dense = rotateSLI_mistake(Nref, uv_xray, O_dense, unit_vec_dense, pivot_point, dist2plane, pos3, p_xray, rot(i, :), [0.5, 0.5, 0.5], 0);
        I2 = [I2; [[I_sparse, j*ones(size(I_sparse, 1), 1)]; [I_dense, j*ones(size(I_dense, 1), 1)]]];
    end
    
    % Translation matrix
    roll =  -rot(j, 1);  % x
    pitch = rot(j, 2); % y
    yaw =   rot(j, 3);   % z

    T = [1, 0, 0, pivot_point(1); 0, 1, 0, pivot_point(2); 0, 0, 1, pivot_point(3); 0, 0, 0, 1];

    Rroll = [1, 0, 0; 0, cos(roll/180*pi), -sin(roll/180*pi); 0, sin(roll/180*pi), cos(roll/180*pi)];
    Rpitch = [cos(pitch/180*pi), 0, sin(pitch/180*pi); 0, 1, 0; -sin(pitch/180*pi), 0, cos(pitch/180*pi)];
    Ryaw = [cos(yaw/180*pi), -sin(yaw/180*pi), 0; sin(yaw/180*pi), cos(yaw/180*pi), 0; 0, 0, 1];

    % Rotation matrix - XYZ
    R = Ryaw*Rpitch*Rroll;
    %R = Rroll*Rpitch*Ryaw;
    R = [R, zeros(3, 1); ones(1, 4).*[0, 0, 0, 1]];
    
    %--------------------------------------------------------------------
%     N = sparseN(sparseN(:, 4) == 13, 1:3);
%     N2 = sparseN(sparseN(:, 4) == j, 1:3);
%     e = cross(N, N2)/norm(cross(N, N2)); % normalised rotation axis; 
%     %ve = subspace(N, Nref'); % rotation angle
%     ve = atan2(norm(cross(N2,N)),dot(N2,N));
%     q = [cos(ve/2), e*sin(ve/2)]; % quaternion
% 
%     R = [1 - 2*(q(3)^2 + q(4)^2), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(1)*q(3) + q(2)*q(4));
%          2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(2)^2 + q(4)^2), 2*(q(3)*q(4) - q(1)*q(2));
%          2*(q(2)*q(4) - q(1)*q(3)), 2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2)];
% 
%     R = [R, zeros(3, 1); ones(1, 4).*[0, 0, 0, 1]];
    %--------------------------------------------------------------------
    
    figure(35)
    hold on
    
    %pos_j(j, :) = Rc*R^-1*Rc^-1*[pos, 1]';
    
    D_sparse_ = T*Rc*R*Rc^-1*T^-1*[0, 0, double(sparseD(j)), 1]';
    D_sparse(j, :) = D_sparse_(1:3);
    
    D_sparse_ = T*Rc^-1*R*Rc^1*T^-1*[0, 0, double(sparseD(j)), 1]';
    D_sparse2(j, :) = D_sparse_(1:3);
    %rotateSLI(Nref, uv_xray, O_sparse, [0, 0, double(sparseD(j))], pivot_point, dist2plane, pos, p_xray, rot(j, 1:3), [0.5, 0.5, 0.5], 0);
    %drawLine([0, 0, 0], D_sparse(j, :), 1, 'm');
    %plot3(D_sparse(j, 1), D_sparse(j, 2), D_sparse(j, 3), 'm.', 'MarkerSize', 15)
    %plot3(D_sparse2(j, 1), D_sparse2(j, 2), D_sparse2(j, 3), 'c.', 'MarkerSize', 15)
    %drawLine([0, 0, 0], D_sparse2(j, :), 1, 'c');
    %pos - D_sparse(j, :)
    
    % Finding the delta-x-y coordinates of the rotations
    
    O_dense_SLI = T*Rc*R*Rc^-1*T^-1*[O_dense, 1]';
    O_dense_SLI = O_dense_SLI(1:3)';
    
    O_sparse_SLI = T*Rc*R*Rc^-1*T^-1*[O_sparse, 1]';
    O_sparse_SLI = O_sparse_SLI(1:3)';
    
    p2_rot = T*Rc*R*Rc^-1*T^-1*[p2, 1]';
    p2_rot = O_dense_SLI(1:3);
    
    uv_xray_rot = Rc*R*Rc^-1*[uv_xray(1:3), 1]';
    uv_xray_rot = uv_xray_rot(1:3)';
    
    [rot_xy_, ~] = line_plane_intersection(uv_xray_rot, p2_rot, Nref, pos3);
    rot_xy = [rot_xy; [j, rot_xy_(1:2)]];
    
    denseQ_coor = denseQ(denseQ(:, 5) == j, :);
    sparseQ_coor = sparseQ(sparseQ(:, 5) == j, :);
    for k = 1:size(denseQ_coor, 1)
        denseQ_r = T*Rc*R*Rc^-1*T^-1*[denseQ_coor(k, 1:3), 1]';
        denseQ_(k, 1:3) = denseQ_r(1:3)';
        rdQ(k, :) = [denseQ_(k, 1:3), denseQ_coor(k, 4:5)];
        %rdQ_ = (denseQ_(k, 1:3) - O_dense_SLI)./norm(denseQ_(k, 1:3) - O_dense_SLI);
        %[rdQ_, ~] = line_plane_intersection(rdQ_, O_dense_SLI, Nref, pos3);
        %rdQ(k, :) = [rdQ_, denseQ_coor(k, 4:5)];
    end
    for k = 1:size(sparseQ_coor, 1)
        sparseQ_r = T*Rc*R*Rc^-1*T^-1*[sparseQ_coor(k, 1:3), 1]';
        sparseQ_(k, 1:3) = sparseQ_r(1:3)';
        rsQ(k, :) = [sparseQ_(k, 1:3), sparseQ_coor(k, 4:5)];
        %rsQ_ = (sparseQ_(k, 1:3) - O_sparse_SLI)./norm(sparseQ_(k, 1:3) - O_sparse_SLI);
        %[rsQ_, ~] = line_plane_intersection(rsQ_, O_sparse_SLI, Nref, pos3);
        %rsQ(k, :) = [rsQ_, sparseQ_coor(k, 4:5)];
    end
    
%     rot_denseQ = [rot_denseQ; denseQ_];
%    rot_sparseQ = [rot_sparseQ; sparseQ_];
    rot_denseQ = [rot_denseQ; rdQ];
    rot_sparseQ = [rot_sparseQ; rsQ];
    
   
    
    %drawLine(10*sparseN(j, 1:3) + D_sparse(j, 1:3), D_sparse(j, 1:3), 1, 'r')
end

%rotateSLI(Nref, uv_xray, O_dense, denseQ(denseQ(:, 5) == j, 1:3), pivot_point, dist2plane, pos, p_xray, rot(i, :), c(i, :), view(i))
%rotateSLI(Nref, uv_xray, O_sparse, denseQ(sparseQ(:, 5) == j, 1:3), pivot_point, dist2plane, pos, p_xray, rot(i, :), c(i, :), view(i))

figure(35)
%axis ij
hold on
xlabel('x-axis $\left[\mathrm{mm}\right]$','interpreter','latex', 'FontSize', 14)
ylabel('y-axis $\left[\mathrm{mm}\right]$','interpreter','latex', 'FontSize', 14)
zlabel('z-axis $\left[\mathrm{mm}\right]$','interpreter','latex', 'FontSize', 14)

%plot3(denseQ(:, 1), denseQ(:, 2), denseQ(:, 3), 'x', 'MarkerSize', 5)
%plot3(sparseQ(:, 1), sparseQ(:, 2), sparseQ(:, 3), 'x', 'MarkerSize', 5)

%plot3(denseQ(denseQ(:, 5) == j, 1), denseQ(denseQ(:, 5) == j, 2), denseQ(denseQ(:, 5) == j, 3), 'r.', 'MarkerSize', 15)
%plot3(sparseQ(sparseQ(:, 5) == j, 1), sparseQ(sparseQ(:, 5) == j, 2), sparseQ(sparseQ(:, 5) == j, 3), 'r.', 'MarkerSize', 15)




I_data = [rot_sparseQ; rot_denseQ];
plot3(I_data(:, 1), I_data(:, 2), I_data(:, 3), 'r.', 'MarkerSize', 12)
%plot3(rot_denseQ(:, 1), rot_denseQ(:, 2), rot_denseQ(:, 3), 'r.', 'MarkerSize', 15)
%plot3(rot_sparseQ(:, 1), rot_sparseQ(:, 2), rot_sparseQ(:, 3), 'r.', 'MarkerSize', 15)


% Data points from model
plot3(I(:, 1), I(:, 2), I(:, 3), 'bo', 'LineWidth', 1.5)
%plot3(I2(:, 1), I2(:, 2), I2(:, 3), 'yo', 'LineWidth', 1.2  )


%plot3(O_sparse(1), O_sparse(2), O_sparse(3), 'g.', 'MarkerSize', 15)
plot3(CP4(:, 1), CP4(:, 2), CP4(:, 3), '.-', 'color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

legend('Measured SLI points', 'Predicted SLI points', 'Target Area','Interpreter','latex');

vars = {'dense', 'sparse', 'dir', 'folder', 'im1', 'im2', 'N_d', 'N_s', 'vars'};
clear(vars{:})

%------------------------------------------------

TO_pred = delaunayTriangulation(I(:, 1), I(:, 2));           
within_vec = ones(size(TO_pred.ConnectivityList, 1), 1);

for i = 1:size(TO_pred.ConnectivityList, 1)
        tri = [TO_pred.Points(TO_pred.ConnectivityList(i, :), :), I(TO_pred.ConnectivityList(i, :), 3)];
        
    if(isDetected2(tri, r)) % If the distance from a vertice to the furthest point is less than five (mm), the object can be detected
        within_vec(i) = 1;
    else
        within_vec(i) = 0;
    end
end

[area_pred, within_vec_best_pred] = coveredArea(TO_pred, within_vec, TO_pred.Points, rad, r, pos3, CP, area_total, Nref);
area_pred
 
%%


figure(56)
%subplot(2, 1, 1)
%pbaspect([1 1.2 1])
hold on
xlabel('x-axis','interpreter','latex', 'FontSize', 14)
ylabel('y-axis','interpreter','latex', 'FontSize', 14)
zlabel('z-axis','interpreter','latex', 'FontSize', 14)
fill(CP(:, 1), CP(:, 2), '.-', 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', .2);
trisurf(TO_pred.ConnectivityList(within_vec_best_pred == 0, :), TO_pred.Points(:, 1), TO_pred.Points(:, 2), zeros(size(TO_pred.Points(:, 2))), 'EdgeColor', [0, 0, 0], 'FaceColor', [0.6350 0.0780 0.1840]) 
trisurf(TO_pred.ConnectivityList(within_vec_best_pred == 1, :), TO_pred.Points(:, 1), TO_pred.Points(:, 2), zeros(size(TO_pred.Points(:, 2))), 'EdgeColor', [0, 0, 0], 'FaceColor', [0.5, 1, 0.3])
%trisurf(TO_data.ConnectivityList, I_data(:, 1), I_data(:, 2), zeros(size(I_data(:, 3))), 'EdgeColor', [0, 0, 0], 'FaceColor', [0.5, 1, 0.3])
plot(CP(:, 1), CP(:, 2), '.-', 'color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

[~, rows] = unique(I_data(:, 1:3), 'rows');
I_data = I_data(rows, :);
I_data = I_data(I_data(:, 4) ~= 255, :);
TO_data = delaunayTriangulation(I_data(:, 1), I_data(:, 2));            
within_vec = ones(size(TO_data.ConnectivityList, 1), 1);
for i = 1:size(TO_data.ConnectivityList, 1)
        tri = [TO_data.Points(TO_data.ConnectivityList(i, :), :), I_data(TO_data.ConnectivityList(i, :), 3)];
        
    if(isDetected2(tri, r)) % If the distance from a vertice to the furthest point is less than five (mm), the object can be detected
        within_vec(i) = 1;
    else
        within_vec(i) = 0;
    end
end

[area_data, within_vec_best_data] = coveredArea(TO_data, within_vec, TO_data.Points, rad, r, pos3, CP, area_total, Nref);
area_data
   
figure(55)
%subplot(2,1,2)
%pbaspect([1 1.2 1])
hold on
xlabel('x-axis','interpreter','latex', 'FontSize', 14)
ylabel('y-axis','interpreter','latex', 'FontSize', 14)
zlabel('z-axis','interpreter','latex', 'FontSize', 14)
fill(CP(:, 1), CP(:, 2), '.-', 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', .2);
trisurf(TO_data.ConnectivityList(within_vec_best_data == 0, :), TO_data.Points(:, 1), TO_data.Points(:, 2), zeros(size(TO_data.Points(:, 2))), 'EdgeColor', [0, 0, 0], 'FaceColor', [0.6350 0.0780 0.1840]) 
trisurf(TO_data.ConnectivityList(within_vec_best_data == 1, :), TO_data.Points(:, 1), TO_data.Points(:, 2), zeros(size(TO_data.Points(:, 2))), 'EdgeColor', [0, 0, 0], 'FaceColor', [0.5, 1, 0.3])
%trisurf(TO_data.ConnectivityList, I_data(:, 1), I_data(:, 2), zeros(size(I_data(:, 3))), 'EdgeColor', [0, 0, 0], 'FaceColor', [0.5, 1, 0.3])
plot(CP(:, 1), CP(:, 2), '.-', 'color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

%plot3(I_data(:, 1), I_data(:, 2), I_data(:, 3), 'r.', 'MarkerSize', 15)
%plot3(rot_denseQ(:, 1), rot_denseQ(:, 2), rot_denseQ(:, 3), 'r.', 'MarkerSize', 15)
%plot3(rot_sparseQ(:, 1), rot_sparseQ(:, 2), rot_sparseQ(:, 3), 'r.', 'MarkerSize', 15)


% Data points from model
%plot3(I(:, 1), I(:, 2), I(:, 3), 'bo', 'LineWidth', 2)



(size(rot_denseQ, 1) + size(rot_sparseQ, 1)) / size(I, 1) * 100
%% - Find the pairwise distance between any measured point and it's equivalent predicted point

% distribution of the spread for the data

mean_img = zeros(size(rot, 1), 1);
std_img = zeros(size(rot, 1), 1);


figure(75)
hold on

set(gca,'FontSize',14)
set(gca,'LineWidth',1)
set(gcf, 'paperunits', 'centimeters', 'Paperposition', [0 0 20 13]);
set(gcf, 'PaperPositionMode', 'auto')
title('Empirical kernel','interpreter','latex', 'FontSize', 15)
xlabel('Spread of the data from its mean $\left[\sigma\right]$','interpreter','latex', 'FontSize', 14)
ylabel('Density','interpreter','latex', 'FontSize', 14)
grid on
box on

pair = [];
g = [];
z = 0;
for l = 1:size(rot, 1)
    pair_ = [];
    m_range = unique(I_data(I_data(:, 5) == l, 4))';
    for m = 1:size(m_range, 2)
        %pair_(m) = 1e3*(I_data(I_data(:, 4) == m_range(m) & I_data(:, 5) == l, 3) - I(I(:, 4) == m_range(m) & I(:, 5) == l, 3));
        pair_(m) = 1e3*norm(I_data(I_data(:, 4) == m_range(m) & I_data(:, 5) == l, 1:3) - I(I(:, 4) == m_range(m) & I(:, 5) == l, 1:3));
    end
    mean_img(l) = mean(pair_);
    std_img(l) = std(pair_);
    
    % z-score
    z = (pair_' - mean_img(l))/std_img(l);
    
    pd=fitdist(z,'Kernel');
    s=(-3:0.01:3);
    y=pdf(  pd,s);
    p=plot(s,y, 'LineWidth', 1.5);
    
    pair = [pair; pair_'];
    g = [g; l*ones(size(pair_', 1), 1)];

end
figure(70)
hold on

set(gca,'FontSize',14)
set(gca,'LineWidth',1)
set(gcf, 'paperunits', 'centimeters', 'Paperposition', [0 0 20 13]);
set(gcf, 'PaperPositionMode', 'auto')
title('Boxplot of spread between model and measured data','interpreter','latex', 'FontSize', 15)
xlabel('Picture number','interpreter','latex', 'FontSize', 14)
ylabel('Distance between model and measured data [$\mu$m]','interpreter','latex', 'FontSize', 14)
grid on
box on
%ylim([-0.025 0.375])
boxplot(pair', g)


%%
 %Nref = -uv_xray;
    %N = sparseN70;
    N = sparseN(sparseN(:, 4) == 13, 1:3);
    N2 = sparseN(sparseN(:, 4) == j, 1:3);
    e = cross(N, N2)/norm(cross(N, N2)); % normalised rotation axis; 
    %ve = subspace(N, Nref'); % rotation angle
    ve = atan2(norm(cross(N2,N)),dot(N2,N));
    q = [cos(ve/2), e*sin(ve/2)]; % quaternion

    Rnorm = [1 - 2*(q(3)^2 + q(4)^2), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(1)*q(3) + q(2)*q(4));
         2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(2)^2 + q(4)^2), 2*(q(3)*q(4) - q(1)*q(2));
         2*(q(2)*q(4) - q(1)*q(3)), 2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2)];

    Rnorm = [Rnorm, zeros(3, 1); ones(1, 4).*[0, 0, 0, 1]];

    Rnorm(1:3, 1:3);
    
    
    roll =  rot(j, 1);  % x
    pitch = rot(j, 2); % y
    yaw =   rot(j, 3);   % z

    Rroll = [1, 0, 0; 0, cos(roll/180*pi), -sin(roll/180*pi); 0, sin(roll/180*pi), cos(roll/180*pi)];
    Rpitch = [cos(pitch/180*pi), 0, sin(pitch/180*pi); 0, 1, 0; -sin(pitch/180*pi), 0, cos(pitch/180*pi)];
    Ryaw = [cos(yaw/180*pi), -sin(yaw/180*pi), 0; sin(yaw/180*pi), cos(yaw/180*pi), 0; 0, 0, 1];

    % Rotation matrix - XYZ
    R = Ryaw*Rpitch*Rroll;
    
    
    ve2 = 2*acos(q(1)); s = sqrt(1 - q(1)^2);
    ex = [q(2)/s, q(3)/s, q(4)/s];
    
    
    
    roll_x = atan2d(Rnorm(2, 3), Rnorm(3, 3)); 
    pitch_y = -asind(Rnorm(1, 3));
    yaw_z = atan2d(Rnorm(1, 2), Rnorm(1, 1));
    
    rot(j, :);
    -[roll_x, pitch_y, yaw_z];
    
%%
%close allx
%plot(rot(:, 1), -rot(:, 2), 'm.', 'MarkerSize', 15)
%axis ij
xlabel('x-axis','interpreter','latex', 'FontSize', 14)
ylabel('y-axis','interpreter','latex', 'FontSize', 14)

hold on
for i = 1:size(rot, 1)
    %drawLine([0, 0, 0], [rot(i, 1), -rot(i, 2), 0], 1, [0.5, 0.5, 0.5]);
end

%view([90 -90])

%%
roll_x = atan2d(Rc(2, 3), Rc(3, 3)); 
pitch_y = -asind(Rc(1, 3));
yaw_z = atan2d(Rc(1, 2), Rc(1, 1));

-[roll_x, pitch_y, yaw_z];
