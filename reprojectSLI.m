function Ix = reprojectSLI(O, uv, p2, h_xray, pivot_point, pos, rot, color, draw)
    roll =  -rot(1);  % x % was -rot(1)
    pitch = -rot(2); % y
    yaw =   -rot(3);   % z
    
    % Create rotation matrix
    R = createRotationMatrix(yaw, pitch, roll);
    R = [R, zeros(3, 1); ones(1, 4).*[0, 0, 0, 1]];
    
    % Rotation to align with base-frame
    best_angle_x = 18; %18.27; % 18
    best_angle_y = 0; %-0.38; % 0
    best_angle_z = 0; % 0.2700; % 0
    Rc = createRotationMatrix(best_angle_z, best_angle_y, best_angle_x);
    Rc = [Rc, zeros(3, 1); ones(1, 4).*[0, 0, 0, 1]];
    
    %T_rot = [1, 0, 0, hexa_pos(1, j); 0, 1, 0, hexa_pos(2, j); 0, 0, 1, hexa_pos(3, j); 0, 0, 0, 1];
    T_p2 = [1, 0, 0, p2(1); 0, 1, 0, p2(2); 0, 0, 1, p2(3); 0, 0, 0, 1];
    T_base = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, h_xray; 0, 0, 0, 1];
    %T_base^-1*T_p2*Rc*T_p2^-1
    %pivot_point
    pivot_point_r = (T_base^-1*T_p2*Rc*T_p2^-1*[pivot_point, 1]')'; 
    %pivot_point_r = (Rc*[pivot_point, 1]')'; 
    % c2x_y = 23%*10^3;        %um
    % x2c_z = 22.997;          %um
    T_pivot = [1, 0, 0, pivot_point_r(1); 0, 1, 0, pivot_point_r(2); 0, 0, 1, pivot_point_r(3); 0, 0, 0, 1];
    
    %RT = R*T_rot*Rc*T_pivot^-1;
    
    %pivot_point = 190; %um
    %c2x_y = 23;        %um
    %x2c_z = 22997*1e-3; %um
    %desired_distance = 30; %um
    %T_pivot = [1, 0, 0, 0; 0, 1, 0, c2x_y; 0, 0, 1, (pivot_point-desired_distance-x2c_z); 0, 0, 0, 1];
    %T = [1, 0, 0, hexa_pos(1, j); 0, 1, 0, hexa_pos(2, j); 0, 0, 1, hexa_pos(3, j); 0, 0, 0, 1];
    %RT = T_pivot*R'*T_pivot^-1*Rc;
    %pos = [0, 0, 235.1198];
    % changed from -> 
    
    RT = T_pivot*R*T_pivot^-1*T_base^-1*T_p2*Rc*T_p2^-1;
    
    pos_rot = RT*[pos, 1]';
    pos_rot = pos_rot(1:3)';
    
    pos = T_base^-1*T_p2*Rc*T_p2^-1*[pos, 1]';
    pos = pos(1:3)';
    
    
    O_SLI = RT*[O, 1]';
    O_SLI = O_SLI(1:3)';
    
    p2_rot = RT*[p2, 1]';
    p2_rot = p2_rot(1:3)';

    Ix = zeros(size(uv(:, :)));
    
% Pure rotation - without translation
%     for i = 1:size(uv, 1)    
%         % Rotate unit vectors (uv) first to the base-frame, then the
%         % specified angles
%         uv_ = R*Rc*[uv(i, 1 :3), 1]';
%         uv_ = uv_(1:3)';
%         % Find out where the rotated unit vectors intersect with the plane
%         [Ix_, ~] = line_plane_intersection(uv_, O_SLI, [0, 0, 1], pos_rot);
%         Ix(i, :) = [Ix_, uv(i, 4)];
%         % If specified with 'draw', the lines can be drawn.
%         drawLine(Ix(i, 1:3), O_SLI, draw, color);
%     end
    
% With translaltion
    for i = 1:size(uv, 1)    
        % Rotate unit vectors (uv) first to the base-frame, then the
        % specified angles
        uv_ = R*Rc*[uv(i, 1 :3), 1]';
        uv_ = uv_(1:3)';
        % Find out where the rotated unit vectors intersect with the plane
%         [Ix_, ~] = line_plane_intersection(uv_, O_SLI, [0, 0, 1], pos); %Nref.*dist2plane + O);
        [Ix_, ~] = line_plane_intersection(uv_, O_SLI, [0, 0, 1], pos_rot);
        N = (pos_rot - p2_rot)./norm(pos_rot - p2_rot);
        % Changed this one -> 
        %N = (pos_rot - O_SLI)./norm(pos_rot - O_SLI);
        [Ix_, ~] = line_plane_intersection(N, Ix_, [0, 0, 1], pos); %Nref.*dist2plane + O);
        Ix(i, :) = [Ix_, uv(i, 4)];
        % If specified with 'draw', the lines can be drawn.
        %drawLine(Ix(i, 1:3), O_SLI, draw, color);
    end
    
    
    % If specified with 'draw', the intersection points can be drawn.
    %if(draw) plot3(Ix(:, 1), Ix(:, 2), Ix(:, 3), '.', 'color', color, 'MarkerSize', 20); end
    