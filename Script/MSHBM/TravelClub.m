%% Pipeline of Xue-Ru's anaysis: Mapping individual networks
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 0 Set up parameters

clear;clc

% setenv('CBIG_CODE_DIR', '/Users/xuerufan/matlab-toolbox/CBIG');
% setenv('FREESURFER_HOME', '/Applications/freesurfer');
% project_dir = '/Users/xuerufan/Downloads/3RB2_test';
% code_dir = '/Users/xuerufan/TravelClub/Script/MSHBM'; 
% temp_dir = '/Users/xuerufan/TravelClub/Script/Templet';
% addpath('/Applications/freesurfer');
% addpath('/Applications/workbench/bin_macosx64/');

%——————————
setenv('CBIG_CODE_DIR', '/media/ubuntu/Zuolab_Xueru/Script/CBIG');
setenv('FREESURFER_HOME', '/home/ubuntu/Softwares/FREESURFER/freesurfer');
project_dir = '/media/ubuntu/Zuolab_Xueru/3RB2_test';
code_dir = '/media/ubuntu/Zuolab_Xueru/Script/MSHBM'; 
temp_dir = '/media/ubuntu/Zuolab_Xueru/Script/Templet';
addpath('/home/ubuntu/Softwares/FREESURFER/freesurfer');
addpath('/home/ubuntu/Softwares/workbench/bin_linux64/');


site = {'3RB2_A', '3RB2_B', '3RB2_C', '3RB2_D', '3RB2_E', '3RB2_G'};
% subid = [1 7 15]; % 注意只有一个site扫描的人去掉
subid = [2 100 200]; % test
seed_mesh = 'fs_LR_900';
targ_mesh = 'fs_LR_32k';
num_clusters = 15;
threshold = 0.1;
niter = 2; % Step4 正式设置为1000


mex -setup C
mex -setup C++
nsite = length(site);
subs = arrayfun(@(x) {sprintf('sub-%03d', x)}, subid);
% subs = arrayfun(@(x) {sprintf('%d', x)}, subid); % test
nsub = length(subid);
cd(code_dir)


disp('Step 0 DONE!');

%% Step 1 prepare the fMRI list files

dataFolder = fullfile(project_dir, 'data');
outputFolder = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params', ...
    'data_list', 'fMRI_list');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

siteFolders = dir(dataFolder);
siteFolders = siteFolders([siteFolders.isdir]);
siteFolders = siteFolders(~ismember({siteFolders.name}, {'.', '..'}));

if length(siteFolders) ~= nsite
    error('Number of site folders doesnt match the expected number of sites. Please check!');
end

for i = 1:nsite
    siteFolderName = siteFolders(i).name;
    siteLetter = siteFolderName(end);
    siteFolderPath = fullfile(dataFolder, siteFolderName);
    
    for j = 1:nsub
        subFolderName = subs{j};
        subfolderPath = fullfile(siteFolderPath, subFolderName);
        
        if exist(subfolderPath, 'dir')
            niiFilePaths = {};

%             runFolders = dir(fullfile(subfolderPath, 'MNINonLinear', 'Results', '*PA_run*'));
            runFolders = dir(fullfile(subfolderPath, 'MNINonLinear', 'Results', 'rfMRI*'));

            for k = 1:length(runFolders)
                runFolderName = runFolders(k).name;
                runFolderPath = fullfile(subfolderPath, 'MNINonLinear', 'Results', ...
                    runFolderName);
                
                niiFile = dir(fullfile(runFolderPath, '*.dtseries.nii'));
                
                if isempty(niiFile)
                    error('This sub doesnt have this data. Please check!');
                end

                niiFilePaths{end+1} = fullfile(runFolderPath, niiFile(1).name);
            end
            
            if ~isempty(niiFilePaths)
                txtFileName = fullfile(outputFolder, sprintf('%s_site-%s.txt', subFolderName, ...
                    siteLetter));
                fileID = fopen(txtFileName, 'w');
                for n = 1:length(niiFilePaths)
                    if n == length(niiFilePaths)
                        fprintf(fileID, '%s', niiFilePaths{n});
                    else
                        fprintf(fileID, '%s ', niiFilePaths{n});
                    end
                end
                fclose(fileID);
            end
        end
    end
end

disp('Step 1 DONE!');

%% Step 2 Generating profiles and initialization parameters
% modified from CBIG_ComputeCorrelationProfile.m
% 由于多回波ICA融合时会进行降噪，这里需要小心设置阈值

out_dir = fullfile(project_dir, 'output/generate_profiles_and_ini_params');

for j = 1:nsub
 for i = 1:nsite

    out_profile_dir = fullfile(out_dir,'profiles', subs{j}, ['site-' site{i}]);

    if ~exist(out_profile_dir)
        mkdir(out_profile_dir);
    end

    profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_' targ_mesh ...
        '_roi' seed_mesh '.surf2surf_profile.mat']);
    fMRI_list = fullfile(out_dir,'data_list','fMRI_list',[subs{j} '_site-' site{i} '.txt']);

    if exist(fMRI_list)
        TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, targ_mesh, profile_file, 'NONE', num2str(threshold), ...
            fMRI_list, 'NONE', 'NONE', '0');
        disp([subs{j} ' site-' site{i} ' DONE!' ]);
    else
        disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ]);
    end
 end
end

disp('Step 2 DONE!');

%% Step 3 To obtain the group averaged profiles
% modified from CBIG_MSHBM_avg_profiles.m

out_dir = fullfile(project_dir, 'output/generate_profiles_and_ini_params');

if(~exist(fullfile(out_dir,'profiles','avg_profile')))
    mkdir(fullfile(out_dir,'profiles','avg_profile'));
end

num_data = 0;
for j = 1:nsub
    for i = 1:nsite

        out_profile_dir = fullfile(out_dir,'profiles', subs{j}, ['site-' site{i}]);       
        avg_profile_file = fullfile(out_dir,'profiles','avg_profile',[targ_mesh '_roi' seed_mesh '_avg_profile.mat']);
        profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_' targ_mesh ...
            '_roi' seed_mesh '.surf2surf_profile.mat']);
    
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
            disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ])
        end
    end
end

profile_mat = avg_profile.profile_mat./num_data;
save(avg_profile_file,'profile_mat','-v7.3');

disp('Step 3 DONE!');

%% Step 4 To run Yeo2011 clustering algorithm for generate or own group prior
% modified from CBIG_MSHBM_generate_ini_params.m

gro_dir = fullfile(project_dir, 'output', 'estimate_group_priors', 'group');

if(~exist(fullfile(gro_dir)))
    mkdir(fullfile(gro_dir));
end

output_file = fullfile(gro_dir,'group.mat');
avg_profile_file = fullfile(out_dir,'profiles','avg_profile',[targ_mesh '_roi' seed_mesh '_avg_profile.mat']);

CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, '', num_clusters, output_file, ...
    avg_profile_file, 'NONE', 0, niter, 0, 'max_iter', 1);

load(output_file);
clustered.mtc = mtc;
clustered.lowerbound = lowerbound;
clustered.lambda = lambda;
save(output_file,'lh_labels','rh_labels','clustered');

disp('Step 4 DONE!');

%% Step 5 把每个被试的功能连接结果路径按照要求整理成txt文件
% 这部分具体的真理要求参阅CBIG_MSHBM_estimate_group_priors.m里的描述

profilesFolder = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params', 'profiles');
outputFolder = fullfile(project_dir, 'output', 'estimate_group_priors', 'profile_list', 'training_set');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for i = 1:nsite
    
    siteFolderName = ['site-' site{i}];
 
    matFilePaths = {};

    txtFileName = fullfile(outputFolder, [siteFolderName, '.txt']);

    for j = 1:nsub
        
        subFolderName = subs{j};
        subfolderPath = fullfile(profilesFolder, subFolderName, siteFolderName);
        
        matFiles = dir(fullfile(subfolderPath, '*.mat'));
        if isempty(matFiles)
            matFilePaths{end+1}  = 'NONE';
        else
            matFilePaths{end+1}  = fullfile(subfolderPath, matFiles(1).name);
        end

        fileID = fopen(txtFileName, 'w');
        for n = 1:length(matFilePaths)
            fprintf(fileID, '%s\n', matFilePaths{n});
        end
        fclose(fileID);
    end  
end

disp('Step 5 DONE!');

%% Step 6 估计模型的参数同时得出个体分区

work_dir = fullfile(project_dir, 'output', 'estimate_group_priors');

% estimate group priors 这里要用雪如修改过的代码
Params = TravelClub_CBIG_MSHBM_estimate_group_priors(work_dir, targ_mesh, num2str(nsub), num2str(nsite), ...
    site, num2str(num_clusters)) % Du 设置了 max_iter 5

% mex -v -largeArrayDims mtimesx.c -lmwblas -lmwlapack 如果需要编译
% Add the result directory and its subfolders to the MATLAB path
addpath(genpath(project_dir));
savepath

disp('Step 6 DONE!');

%% Step 7 提取出个体分区的结果label
% modified from CBIG_IndCBM_extract_MSHBM_result.m

final = load(fullfile(work_dir, 'priors', 'Params_Final.mat'));
vertex_num = size(final.Params.s_lambda, 1) / 2;

output_dir = fullfile(project_dir, 'output', 'ind_parcellation');
if(~exist(output_dir))
    mkdir(output_dir);
end

for j = 1:nsub
    if(isempty(final.Params.s_lambda))
        error('s_lambda is empty. Please use save_all flag when estimating group prior.')
    end

    sub = final.Params.s_lambda(:, :, j);

    [~, labels] = max(sub');
    labels(sum(sub, 2) == 0) = 0; % Label medial wall as 0
    lh_labels = labels(1:vertex_num)';
    rh_labels = labels((vertex_num+1):end)';

    output_file = fullfile(output_dir, ['Ind_parcellation_MSHBM_' subs{j} '.mat']);
    save(output_file, 'lh_labels', 'rh_labels', 'num_clusters');
end

disp('Step 7 DONE!');

%% Step 8 可视化

vis_dir = fullfile(project_dir, 'visual');
if ~exist(vis_dir, 'dir')
    mkdir(vis_dir);
end

data_dir = fullfile(project_dir, 'output', 'ind_parcellation');
for j = 1:nsub

    load(fullfile(data_dir, ['Ind_parcellation_MSHBM_', subs{j}, '.mat']));

    temp = cifti_read(fullfile(temp_dir, 'DU15NET_consensus_fsLR_32k.dlabel.nii'));
    temp.cdata = [lh_labels; rh_labels];

    fn = fullfile(vis_dir, [subs{j}, '_', num2str(num_clusters), 'Net_HCP.dlabel.nii']); 
    cifti_write(temp, fn);
end

disp("Step 8 DONE!")
