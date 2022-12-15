
function [O_line1, O_line2] = closestPoint(line1, line2)
    O_line1 = line1(find(vecnorm(abs(line1 - line2)') == min(vecnorm(abs(line1 - line2)'))), :);
    O_line2 = line2(find(vecnorm(abs(line1 - line2)') == min(vecnorm(abs(line1 - line2)'))), :);
end