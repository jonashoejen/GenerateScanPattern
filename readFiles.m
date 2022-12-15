%clearvars
%close all

function[fig, im_unk1, im_unk2] = readFiles(myDir, fig, showImg)

% myDir = uigetdir; %gets directory - forces you to decide yourself
% or just use the hardcoded path

myFiles_unc = dir(fullfile(myDir,'*.unc')); %gets all unc files in struct
myFiles_unk = dir(fullfile(myDir,'*.unk')); %gets all unk files in struct

T_unc = struct2table(myFiles_unc); % convert the struct array to a table
T_unc_sorted = sortrows(T_unc, 'date'); % sort the table by 'DOB'
myFiles_unc = table2struct(T_unc_sorted); 

T_unk = struct2table(myFiles_unk); % convert the struct array to a table
T_unk_sorted = sortrows(T_unk, 'date'); % sort the table by 'DOB'
myFiles_unk = table2struct(T_unk_sorted);

for k = 1:min([length(myFiles_unc), length(myFiles_unk)])
    baseFileName = myFiles_unc(k).name;
    fullFileName = fullfile(myDir, baseFileName);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    if (k < length(myFiles_unk) + 1)
        [header, im, metaData] = readAscImage(fullFileName);
        eval(['header_unc' num2str(k) '= header;']);
        eval(['im_unc' num2str(k) '= im;']);
        eval(['metaData_unc' num2str(k) '= metaData;']);
        if showImg 
            figure(k + fig)
            imshow(im);
        end
    end
end

% Only the first two are relevant in the first case
for k = 1:min([length(myFiles_unc), length(myFiles_unk)])
    baseFileName = myFiles_unk(k).name;
    fullFileName = fullfile(myDir, baseFileName);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    %figure(k)
    [header, im, metaData] = readAscImage(fullFileName);
    eval(['header_unk' num2str(k) '= header;']);
    eval(['im_unk' num2str(k) '= im;']);
    eval(['metaData_unk' num2str(k) '= metaData;']);
    
    if showImg 
        figure(k + fig)
        hold on
        for m = 1:eval(['im_unk' num2str(k) '.sliHeader.Ncen;'])
            cen = eval(['im_unk' num2str(k) '.sliData.cen(m, :);']);
            if (eval(['im_unk' num2str(k) '.sliData.sliID(m);']) ~= 255)
                plt = plot(cen(1), cen(2), 'wo'); 
                plt.MarkerSize = 5;
                plt.LineWidth = 1.25;
                txt = text(cen(1) + 7, cen(2), strcat('\leftarrow ', num2str(m), ', ', num2str(cen(3))));
                txt.Color = 'white';
                txt.FontSize = 10;
            else
                plt = plot(cen(1), cen(2), 'ro'); 
                plt.MarkerSize = 5;
                plt.LineWidth = 1.25;
                txt = text(cen(1) + 7, cen(2), strcat('\leftarrow ', num2str(m), ', ', num2str(cen(3))));
                txt.Color = 'red';
                txt.FontSize = 10;
            end
        end
    end
end

fig = length(myFiles_unk);

vars = {'header', 'im', 'metaData', 'k', 'myDir', 'baseFileName', 'fullFileName', 'myFiles_unc', 'myFiles_unk', 'vars'};
clear(vars{:});

end
