function TravelClub_MSHBM_ProfilePrep(project_dir)
% This functional prepare all the txt files for the input of group priors
% estimation of MSHBM
% Input:
%   - project_dir:
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University
% Email:xueru@mail.bnu.edu.cn

    % Define the input and output folders
    profilesFolder = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params', 'profiles');
    outputFolder = fullfile(project_dir, 'output', 'estimate_group_priors', 'profile_list', 'training_set');
    
    % Create the output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Get all subject folders
    subfolders = dir(profilesFolder);
    subfolders = subfolders([subfolders.isdir]);
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
    
    % Filter out non-subject folders (only keep folders starting with 'sub')
    subfolders = subfolders(startsWith({subfolders.name}, 'sub'));
    
    % Get all session folders for the first subject to determine the number of sessions
    if isempty(subfolders)
        error('No subject folders found in the profiles directory.');
    end
    firstSubFolder = fullfile(profilesFolder, subfolders(1).name);
    sessionFolders = dir(firstSubFolder);
    sessionFolders = sessionFolders([sessionFolders.isdir]);
    sessionFolders = sessionFolders(~ismember({sessionFolders.name}, {'.', '..'}));
    
    % Iterate through each session folder
    for j = 1:length(sessionFolders)
        sessionFolderName = sessionFolders(j).name;
        
        % Initialize a cell array to store paths for each subject's session mat file
        matFilePaths = cell(length(subfolders), 1);
        
        % Iterate through each subject folder
        for i = 1:length(subfolders)
            subfolderName = subfolders(i).name;
            subfolderPath = fullfile(profilesFolder, subfolderName, sessionFolderName);
            
            % Check if the mat file exists
            matFiles = dir(fullfile(subfolderPath, '*.mat'));
            if isempty(matFiles)
                matFilePaths{i} = 'NONE';
            else
                matFilePaths{i} = fullfile(subfolderPath, matFiles(1).name);
            end
        end
        
        % Create a txt file for the session and write all mat file paths
        txtFileName = fullfile(outputFolder, [sessionFolderName, '.txt']);
        fileID = fopen(txtFileName, 'w');
        for k = 1:length(matFilePaths)
            fprintf(fileID, '%s\n', matFilePaths{k});
        end
        fclose(fileID);
    end
    
    disp('All txt files have been created and saved in the output/estimate_group_priors/profile_list/training_set folder.');
end
