function [O_sparse, unit_vec_sparse, O_dense, unit_vec_dense] = analyseData()
%clearvars
%close all

folder = 'C:\Users\jonas\OneDrive - Danmarks Tekniske Universitet\Skrivebord\Tidligere kurser & bøger\30220 Synthesis in Earth and Space Physics Fall 22\Repo\Jonas_Init_meas\';
dir_ = {'65\2022_09_08_14_37_12\download\', '70\2022_09_08_14_41_52\download', '75\2022_09_08_14_46_41\download'};
showImg = 0;
fig = 0;
fig_add = 0;

% create 3D point cloud
d = [65, 70, 75];
d_plot = {'m.', 'y.', 'c.'};
s_plot = {'r.', 'g.', 'b.'};
for i = 1:length(d)
    [fig_add, im1, im2] = readFiles(strcat(folder, dir_{i}), fig, showImg);
    dense = [im1.sliData.Q(im1.sliData.sliID ~= 255, :)/1e3, im1.sliData.sliID(im1.sliData.sliID ~= 255, :)];
    sparse = [im2.sliData.Q(im2.sliData.sliID ~= 255, :)/1e3, im2.sliData.sliID(im2.sliData.sliID ~= 255, :)];
    
    N_d = im1.sliHeader.planeN;
    D_d = im1.sliHeader.planeD/1e3;
    N_s = im2.sliHeader.planeN;
    D_s = im2.sliHeader.planeD/1e3;
    
    eval(['denseQ' num2str(d(i)) '= dense;']);
    eval(['sparseQ' num2str(d(i)) '= sparse;']);
    eval(['denseN' num2str(d(i)) '= N_d;']);
    eval(['sparseN' num2str(d(i)) '= N_s;']);
    eval(['denseD' num2str(d(i)) '= D_d;']);
    eval(['sparseD' num2str(d(i)) '= D_s;']);
end

% figure(75)
% set(gca,'FontSize',14)
% set(gca,'LineWidth',1)
% set(gcf, 'paperunits', 'centimeters', 'Paperposition', [0 0 20 13]);
% set(gcf, 'PaperPositionMode', 'auto')
% title('Traced paths of the SLI points as taken from a distance of 65, 70 and 75 mm','interpreter','latex', 'FontSize', 15)
% xlabel('x-axis','interpreter','latex', 'FontSize', 14)
% ylabel('y-axis','interpreter','latex', 'FontSize', 14)
% grid on
% box on

% Sparse array (RED)
[O_sparse, unit_vec_sparse] = findOrigin(sparseQ65, sparseQ70, sparseQ75, 0, [0.6350 0.0780 0.1840]); 
% Dense array (Green)
[O_dense, unit_vec_dense] = findOrigin(denseQ65, denseQ70, denseQ75, 0, [0.4660 0.6740 0.1880]); 
% Q65 = ImageCorrection(sparseQ65); Q70 = ImageCorrection(sparseQ70); Q75 = ImageCorrection(sparseQ75); 