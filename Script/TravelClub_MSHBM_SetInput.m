function TravelClub_MSHBM_SetInput()

% Input:
%     - seed_mesh: fs_LR_900
%     - targ_mesh: fs_LR_32k
%     - nsub: 3
%     - maxsess: 2
%     - split_flag: 0
%     - output_dir: E:\PhDproject\HCP\HCP_test4MSHBM\output
   
% Prompt the user for input parameters
    seed_mesh = input('Enter the resolution of the ROI mesh: ', 's');
    targ_mesh = input('Enter the surface mesh of input fMRI data: ', 's');
    nsub = input('Enter the total number of subjects: ', 's');
    maxsess = input('Enter the maxmium mumber of sessions: ', 's');
    split_flag = input('Enter the split flag: ', 's');
    output_dir = input('Enter the output directory: ', 's');

    % Create project_dir
    project_dir = fullfile(output_dir, 'generate_profiles_and_ini_params');
    
    % Display the input parameters
    fprintf('Input Parameters:\n');
    fprintf('seed_mesh: %s\n', seed_mesh);
    fprintf('targ_mesh: %s\n', targ_mesh);
    fprintf('nsub: %s\n', nsub);
    fprintf('maxsess: %s\n', maxsess);
    fprintf('split_flag: %s\n', split_flag);
    fprintf('project_dir: %s\n', project_dir);
    
    % Assign variables to the base workspace
    assignin('base', 'seed_mesh', seed_mesh);
    assignin('base', 'targ_mesh', targ_mesh);
    assignin('base', 'output_dir', output_dir);
    assignin('base', 'nsub', nsub);
    assignin('base', 'maxsess', maxsess);
    assignin('base', 'split_flag', split_flag);
    assignin('base', 'project_dir', project_dir);
end
