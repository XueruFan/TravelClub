function TravelClub_MSHBM_Prep(project_dir, num_sessions, num_subs, num_runs)
% This function prepares all the txt files for the input of MSHBM
% Input:
%   - project_dir: The directory containing the data
%   - num_sessions: The total number of sessions expected
%   - num_subs: The total number of subjects expected
%   - num_runs: The total number of runs expected
% Note: This code is designed for our in-house 3RB2 data storage format, site-sub-run
% Written by Xue-Ru Fan 2024-06-01 @ Beijing Normal University

    % Define the data and output folders
    dataFolder = fullfile(project_dir, 'data');
    outputFolder = fullfile(project_dir, 'output', ...
        'generate_profiles_and_ini_params', 'data_list', 'fMRI_list');
    
    % Create the output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Get all session folders
    sessionFolders = dir(dataFolder);
    sessionFolders = sessionFolders([sessionFolders.isdir]);
    sessionFolders = sessionFolders(~ismember({sessionFolders.name}, {'.', '..'}));
    
    % Check the number of sessions
    if length(sessionFolders) ~= num_sessions
        error('The number of session folders does not match the expected number of sessions.');
    end

    % Define subject names
    subs = {'sub-001', 'sub-002', 'sub-003', 'sub-004', 'sub-007'};
    
    % Iterate through each session folder
    for i = 1:num_sessions
        sessionFolderName = sessionFolders(i).name;
        sessionLetter = sessionFolderName(end); % Extract the last letter of the session folder name
        sessionFolderPath = fullfile(dataFolder, sessionFolderName);
        
        % Iterate through each subject
        for j = 1:num_subs
            subfolderName = subs{j};
            subfolderPath = fullfile(sessionFolderPath, subfolderName);
            
            % Initialize a cell array to store all run file paths for the subject
            niiFilePaths = cell(1, num_runs);
            
            if exist(subfolderPath, 'dir')
                for r = 1:num_runs
                    % Find folders that match the pattern *run%d
                    pattern = sprintf('*run%d', r);
                    folderInfo = dir(fullfile(subfolderPath, 'MNINonLinear', 'Results', pattern));
                    
                    if isempty(folderInfo)
                        % If no folder is found, use 'NONE'
                        niiFilePaths{r} = 'NONE';
                    else
                        % Use the first matching folder
                        runFolderName = folderInfo(1).name;
                        runFolderPath = fullfile(subfolderPath, 'MNINonLinear', 'Results', runFolderName);
                        
                        % Look for the specific file in the subfolder
                        niiFile = dir(fullfile(runFolderPath, '*.dtseries.nii'));
                        
                        if isempty(niiFile)
                            % If the file is not found, use 'NONE'
                            niiFilePaths{r} = 'NONE';
                        else
                            % Add the file path to the cell array
                            niiFilePaths{r} = fullfile(runFolderPath, niiFile(1).name);
                        end
                    end
                end
            else
                % If the subfolder does not exist, use 'NONE' for all runs
                niiFilePaths(:) = {'NONE'};
            end
            
            % Combine all run file paths for this session into a single string
            niiFilePathsStr = strjoin(niiFilePaths, ' ');
            
            % Create a txt file for the subject and write all nii file paths for all runs
            txtFileName = fullfile(outputFolder, sprintf('%s_sess%s.txt', subfolderName, sessionLetter));
            fileID = fopen(txtFileName, 'w');
            fprintf(fileID, '%s\n', niiFilePathsStr);
            fclose(fileID);
        end
    end
    
    disp('DONE!');
end
