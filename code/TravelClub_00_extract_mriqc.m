% Extract MRIQC result
% Xueru Fan Jun 3 2025 @ BNU

clear; clc;

rootDir = 'E:\PhDproject\3RB2\QC';
cd(rootDir)

siteDirs = dir(fullfile(rootDir, '3RB2*'));
siteDirs = siteDirs([siteDirs.isdir]);

resultTable = table('Size', [0, 6], ...
                   'VariableTypes', {'string', 'string', 'string', 'double', 'double', 'double'}, ...
                   'VariableNames', {'sub', 'site', 'run', 'fd_mean', 'snr', 'tsnr'});

for s = 1:length(siteDirs)

    sitePath = fullfile(siteDirs(s).folder, siteDirs(s).name);

    % site
    siteName = extractAfter(siteDirs(s).name, '3RB2_');
    
    subDirs = dir(fullfile(sitePath, 'sub-*'));
    subDirs = subDirs([subDirs.isdir]);
    
    for sub = 1:length(subDirs)
        subPath = fullfile(subDirs(sub).folder, subDirs(sub).name);
        funcDir = fullfile(subPath, 'func');
      
        jsonFiles = dir(fullfile(funcDir, '*.json'));
        
        for f = 1:length(jsonFiles)

            jsonFile = fullfile(jsonFiles(f).folder, jsonFiles(f).name);
            jsonData = jsondecode(fileread(jsonFile));
            
            % run
            [~, filename] = fileparts(jsonFiles(f).name);
            match = regexp(filename, 'task-rest_(.*?)_bold$', 'tokens');
            runInfo = match{1}{1};

            % 提取所需数据
            newRow = {...
                string(jsonData.bids_meta.subject_id), ...
                siteName, ...
                runInfo, ...
                jsonData.fd_mean, ...
                jsonData.snr, ...
                jsonData.tsnr ...
            };
            
            resultTable = [resultTable; newRow];
        end
    end
end

writetable(resultTable, fullfile(rootDir, 'QC_result.xlsx'))

% 筛选并保存fd_mean > 0.2
excludeTable = resultTable(resultTable.fd_mean > 0.2, :);
writetable(excludeTable, fullfile(rootDir, 'QC_exclude.xlsx'));

disp("DONE!")