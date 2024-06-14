%% Pipeline of Xue-Ru's anaysis: Mapping individual networks
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University
%% Step 0
clear;clc
%—————————————————————————— 下面这里需要手动定义
project_dir = '/Users/xuerufan/Downloads/3RB2_test';
code_dir = '/Users/xuerufan/TravelClub/Script';
site = {'A', 'C', 'D'};
subid = [1 4 7]; % 注意这里列出要用的也就是数据全的被试编号

%—————————————————————————— 以上需要自己定义，下面不用管
mex -setup C
mex -setup C++
nsite = length(site);
subs = arrayfun(@(x) {sprintf('sub-%03d', x)}, subid);
subnums = arrayfun(@(x) {sprintf('-%03d', x)}, subid);
nsub = length(subid);
cd(code_dir)  
%% Step 1 prepare the fMRI list files
dataFolder = fullfile(project_dir, 'data');
outputFolder = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params', 'data_list', 'fMRI_list');

% Create the output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Get all session folders
sessionFolders = dir(dataFolder);
sessionFolders = sessionFolders([sessionFolders.isdir]);
sessionFolders = sessionFolders(~ismember({sessionFolders.name}, {'.', '..'}));

% Check the number of sessions
if length(sessionFolders) ~= nsite
    error('The number of session folders does not match the expected number of sessions.');
end

% Iterate through each session folder
for i = 1:nsite
    sessionFolderName = sessionFolders(i).name;
    sessionLetter = sessionFolderName(end); % Extract the last letter of the session folder name
    sessionFolderPath = fullfile(dataFolder, sessionFolderName);
    
    % Iterate through each subject
    for j = 1:nsub
        subfolderName = subs{j};
        subfolderPath = fullfile(sessionFolderPath, subfolderName);
        
        if exist(subfolderPath, 'dir')
            % Initialize a cell array to store all run file paths for the subject
            niiFilePaths = {};

            % Find all run folders that match the pattern *run*
            runFolders = dir(fullfile(subfolderPath, 'MNINonLinear', 'Results', '*run*'));

            for k = 1:length(runFolders)
                runFolderName = runFolders(k).name;
                runFolderPath = fullfile(subfolderPath, 'MNINonLinear', 'Results', runFolderName);
                
                % Look for the specific file in the subfolder
                niiFile = dir(fullfile(runFolderPath, '*.dtseries.nii'));
                
                if ~isempty(niiFile)
                    % Add the file path to the cell array
                    niiFilePaths{end+1} = fullfile(runFolderPath, niiFile(1).name);
                end
            end
            
            % Only create a txt file if there are valid nii files
            if ~isempty(niiFilePaths)
                % Create a txt file for the subject and write all nii file paths for all runs
                txtFileName = fullfile(outputFolder, sprintf('%s_sess%s.txt', subfolderName, sessionLetter));
                fileID = fopen(txtFileName, 'w');
                for n = 1:length(niiFilePaths)
                    fprintf(fileID, '%s\n', niiFilePaths{n});
                end
                fclose(fileID);
            end
        end
    end
end

disp('DONE!');


%% Step 2 set input
% see Kong2019_MSHBM README.md for input suggestions
seed_mesh = 'fs_LR_900';
targ_mesh = 'fs_LR_32k';
split_flag = '0';
out_dir = fullfile(project_dir, 'output/generate_profiles_and_ini_params');
%% Step 3 Generating profiles and initialization parameters
% make sure you have added the CBIG folder into your path
for su = 1:nsub
 for se = 1:nsite
    CBIG_MSHBM_generate_profiles(seed_mesh, targ_mesh, out_dir, ...
        subnums{su}, site{se}, split_flag);
 end
end
%% Step 4 To obtain the group averaged profiles
% modified from CBIG_MSHBM_avg_profiles.m

if(~exist(fullfile(out_dir,'profiles','avg_profile')))
    mkdir(fullfile(out_dir,'profiles','avg_profile'));
end
num_data = 0;
for su = 1:nsub
    out_profile_dir = fullfile(out_dir,'profiles',['sub' subnums{su}],['sess' site{se}]);       
    avg_profile_file = fullfile(out_dir,'profiles','avg_profile',[targ_mesh '_roi' seed_mesh '_avg_profile.mat']);
    profile_file = fullfile(out_profile_dir,['sub' subnums{su} '_sess' site{se} '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.mat']);      
    if(exist(profile_file))
        num_data = num_data + 1;
        profile_data = load(profile_file);
        if(num_data == 1)
           avg_profile = profile_data;
        else
           avg_profile.profile_mat = avg_profile.profile_mat + profile_data.profile_mat;
        end
    else
        fprintf('Skip: %s \n',profile_file);
    end
end
profile_mat = avg_profile.profile_mat./num_data;
save(avg_profile_file,'profile_mat','-v7.3');

%% Step 5 To run Yeo2011 clustering algorithm for generate or own group prior
% modified from CBIG_MSHBM_generate_ini_params.m
% 这里用的就是自己的数据
num_clusters = 15;
num_initialization = 2; % !!!!!!!!!!!!!use 1000 for real analysis
% 这一步Du用HCP的40个人的数据，从2~21个网络都run了一遍，各做了一个组水平的
% 先验，最终是结合后面基于ROI的连接图，确定了15个网络

if(~exist(fullfile(out_dir,'group')))
    mkdir(fullfile(out_dir,'group'));
end

output_file = fullfile(out_dir,'group','group.mat');
avg_profile_file = fullfile(out_dir,'profiles','avg_profile',[targ_mesh '_roi' seed_mesh '_avg_profile.mat']);
CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, '', num_clusters, output_file, ...
avg_profile_file, 'NONE', 0, num_initialization, 0, 100, 1);
    
% Reorganize output variables
if(exist(output_file))
    load(output_file);
    clustered.mtc = mtc;
    clustered.lowerbound = lowerbound;
    clustered.lambda = lambda;
    save(output_file,'lh_labels','rh_labels','clustered');
else
    error('could not find clustering results group.mat')
end

%% Step 6 把每个被试的功能连接结果路径按照要求整理成txt文件
% 这部分具体的真理要求参阅CBIG_MSHBM_estimate_group_priors.m里的描述
TravelClub_MSHBM_ProfilePrep(project_dir);
%% Step 7 把生成的group.mat文件和文件夹拷贝到MSHBM_Step2的文件夹里
TravelClub_MSHBM_CopyGroup(project_dir)
%% Step 8 估计模型的参数同时得出个体分区
% set up input
work_dir = fullfile(project_dir, 'output', 'estimate_group_priors');
mesh = 'fs_LR_32k';
num_sub = num2str(nsub);
num_sess = num2str(nsite);
num_clusters = '15';
% estimate group priors 这里要用雪如修改过的代码
Params = TravelClub_CBIG_MSHBM_estimate_group_priors(work_dir, mesh, num_sub, ...
    num_sess, site, num_clusters);

% mex -v -largeArrayDims mtimesx.c -lmwblas -lmwlapack 如果需要编译
% Add the result directory and its subfolders to the MATLAB path
addpath(genpath(project_dir));
savepath
%% Step 9 提取出个体分区的结果label
% indi_dir = fullfile(project_dir, 'output', 'generate_individual_parcellations');   
% % Create the output folder if it doesn't exist
% if ~exist(indi_dir, 'dir')
%     mkdir(indi_dir);
% end
% addpath(genpath(indi_dir));
% savepath
% 跳过Kong2019_Step3为新个体估计个体网络的部分
% 使用Xue2021的提取MSHBM结果的代码
work_dir = fullfile(project_dir, 'output', 'estimate_group_priors');
CBIG_IndCBM_extract_MSHBM_result(work_dir)
%% Step 10 可视化
project_dir = 'E:\PhDproject\HCP\HCP_test4MSHBM';
vis_dir = fullfile(project_dir, 'visual'); % dscalar模板文件的地址
data_dir = fullfile(project_dir, 'output\estimate_group_priors\ind_parcellation');
num_sub = '3';
num_clusters = '15';
TravelClub_MSHBM_VisNet(data_dir, vis_dir, num_sub, num_clusters)