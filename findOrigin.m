function [O, unit_vec] = findOrigin(Q65, Q70, Q75, draw, color)    
    % y-axis distortion factor
    %y_df = (8.3/8.6);
    %Q65(:, 2) = Q65(:, 2).*y_df;
    uniqueIDs = unique(sort([Q65(:, 4);Q70(:, 4);Q75(:, 4)]));
    line = []; %zeros(length(0:0.1:100), length(uniqueIDs));
    unit_vec = [];
    if(draw)
        figure(75)
        hold on
        xlabel('x-axis')
        ylabel('y-axis')
        zlabel('z-axis')
    end
    for i = 1:length(uniqueIDs)
        p1 = ismember(uniqueIDs(i), Q65(:, 4));
        p2 = ismember(uniqueIDs(i), Q70(:, 4));
        p3 = ismember(uniqueIDs(i), Q75(:, 4));
        if(p1 + p2 + p3 > 1)
            if(p1)
                start = Q65(Q65(:, 4) == uniqueIDs(i), 1:3);
            else
                start = Q70(Q70(:, 4) == uniqueIDs(i), 1:3);
            end
            if(p3)
                ending = Q75(Q75(:, 4) == uniqueIDs(i), 1:3);
            else
                ending = Q70(Q70(:, 4) == uniqueIDs(i), 1:3);
            end
        %dL = drawLine(ending, start, dz);
        %line = [line, dL];
        line = [line, drawLine(ending, ending - 100*(ending - start)/norm(ending - start), draw, color)];
        %line = [line, drawLine(ending, start, 20, 1)];
        unit_vec = [unit_vec; (ending - start)/norm(ending - start), uniqueIDs(i)];
        end
    end

    n = size(line, 2)/3;
    coords = zeros((n)*(n+1)/2, 3);
    c_count = 1;
    for i = 1:n - 1
        for j = 1:n
            if(i ~= j && j > i)
                [cl_l1, cl_l2] = closestPoint(line(:,i*3 - 2:i*3), line(:,j*3-2:j*3));
                coords(2*c_count - 1:2*c_count, :) = [cl_l1;cl_l2];
                c_count = c_count + 1;
            end
        end
    end
    O = mean(coords, 1);
    
    unit_vec = unit_vec(ismember(unit_vec(:, 4), excludePoints()), :);
 
%     illPoints = uniqueIDs(~ismember(unique(sort([Q65(:, 4);Q70(:, 4);Q75(:, 4)])), unit_vec(:, 4)));
%     
%     for k = 1:length(illPoints)
%         illQ65 = Q65(Q65(:, 4) == illPoints(k), 1:3);
%         illQ70 = Q70(Q70(:, 4) == illPoints(k), 1:3);
%         illQ75 = Q75(Q75(:, 4) == illPoints(k), 1:3);
%         
%         if(~isempty(illQ65))
%             unit_vec = [unit_vec; (illQ65 - O)/norm(illQ65 - O), illPoints(k)];
%         end
%         if(~isempty(illQ70))
%             unit_vec = [unit_vec; (illQ70 - O)/norm(illQ70 - O), illPoints(k)];
%         end
%         if(~isempty(illQ75))
%             unit_vec = [unit_vec; (illQ75 - O)/norm(illQ75 - O), illPoints(k)];
%         end
%     end
    %unit_vec;
    end





