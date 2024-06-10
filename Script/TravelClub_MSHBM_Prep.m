function TravelClub_MSHBM_Prep(project_dir)
% This functional prepare all the txt files for the input of MSHBM
% Input:
%   - project_dir:
% Written by Xue-Ru Fan 2024-06-01 @ Beijing Normal University
% Email:xueru@mail.bnu.edu.cn

    % Define the data and output folders
    dataFolder = fullfile(project_dir, 'data');
    outputFolder = fullfile(project_dir, 'output', ...
        'generate_profiles_and_ini_params', 'data_list', 'fMRI_list');
    
    % Create the output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Get all subject folders
    subfolders = dir(dataFolder);
    subfolders = subfolders([subfolders.isdir]);
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
    
    % Iterate through each subject folder
    for i = 1:length(subfolders)
        subfolderName = subfolders(i).name;
        subfolderPath = fullfile(dataFolder, subfolderName);
        
        % Get all session folders in each subject folder
        sessionFolders = dir(subfolderPath);
        sessionFolders = sessionFolders([sessionFolders.isdir]);
        sessionFolders = sessionFolders(~ismember({sessionFolders.name}, {'.', '..'}));
        
        % Iterate through each session folder
        for j = 1:length(sessionFolders)
            sessionFolderName = sessionFolders(j).name;
            sessionFolderPath = fullfile(subfolderPath, sessionFolderName);
            
            % Initialize a cell array to store all run file paths for the session
            niiFilePaths = {};
            
            % Get all run folders in each session folder
            runFolders = dir(sessionFolderPath);
            runFolders = runFolders([runFolders.isdir]);
            runFolders = runFolders(~ismember({runFolders.name}, {'.', '..'}));
            
            % Iterate through each run folder to get nii file paths
            for k = 1:length(runFolders)
                runFolderName = runFolders(k).name;
                runFolderPath = fullfile(sessionFolderPath, runFolderName);
                
                % Get nii files
                niiFiles = dir(fullfile(runFolderPath, '*.nii*'));
                
                % Add nii file paths to the cell array
                for m = 1:length(niiFiles)
                    niiFilePath = fullfile(runFolderPath, niiFiles(m).name);
                    niiFilePaths{end+1} = niiFilePath; %#ok<AGROW>
                end
            end
            
            % Create a txt file for the session and write all nii file paths
            txtFileName = fullfile(outputFolder, [subfolderName, '_', sessionFolderName, '.txt']);
            fileID = fopen(txtFileName, 'w');
            fprintf(fileID, '%s ', niiFilePaths{:});
            fclose(fileID);
        end
    end
    
    disp('All txt files have been created and saved in the output/data_list/fMRI_list folder.');
end
