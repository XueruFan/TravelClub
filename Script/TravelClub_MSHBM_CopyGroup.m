function TravelClub_MSHBM_CopyGroup(project_dir)
% This functional copy the group.m to the right folder
% Input:
%   - project_dir:
% Written by Xue-Ru Fan 2024-06-01 @ Beijing Normal University
% Email:xueru@mail.bnu.edu.cn

    % Define the source and destination directories
    sourceFolder = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params', 'group');
    destinationFolder = fullfile(project_dir, 'output', 'estimate_group_priors');
    
    % Ensure the source folder exists
    if ~exist(sourceFolder, 'dir')
        error('Source folder does not exist: %s', sourceFolder);
    end
    
    % Create the destination folder if it does not exist
    if ~exist(destinationFolder, 'dir')
        mkdir(destinationFolder);
    end
    
    % Define the full destination path including the group folder
    destinationGroupFolder = fullfile(destinationFolder, 'group');
    
    % Copy the source folder to the destination
    copyfile(sourceFolder, destinationGroupFolder);
    
    disp('Group folder and its contents have been copied successfully.');
end