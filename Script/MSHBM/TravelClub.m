%% Pipeline of Xue-Ru's anaysis: Mapping individual networks
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 0 Set up parameters

clear;clc

%—————————————————————————— 下面这里需要手动定义

project_dir = '/Users/xuerufan/Downloads/3RB2_test';
code_dir = '/Users/xuerufan/TravelClub/Script/MSHBM'; 
vis_dir = '/Users/xuerufan/TravelClub/Script/Templet';

site = {'A', 'C', 'D'};
subid = [1 2 4 15 17];


%—————————————————————————— 下面不用管
mex -setup C
mex -setup C++
nsite = length(site);
subs = arrayfun(@(x) {sprintf('sub-%03d', x)}, subid);
subnums = arrayfun(@(x) {sprintf('%03d', x)}, subid);
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
        subfolderName = subs{j};
        subfolderPath = fullfile(siteFolderPath, subfolderName);
        
        if exist(subfolderPath, 'dir')
            niiFilePaths = {};

            runFolders = dir(fullfile(subfolderPath, 'MNINonLinear', 'Results', '*PA_run*'));

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
                txtFileName = fullfile(outputFolder, sprintf('%s_site-%s.txt', subfolderName, ...
                    siteLetter));
                fileID = fopen(txtFileName, 'w');
                for n = 1:length(niiFilePaths)
                    fprintf(fileID, '%s\n', niiFilePaths{n});
                end
                fclose(fileID);
            end
        end
    end
end

disp('Step 1 DONE!');

%% Step 2 set input
% see Kong2019_MSHBM README.md for input suggestions

seed_mesh = 'fs_LR_900';
targ_mesh = 'fs_LR_32k';
split_flag = '0';
out_dir = fullfile(project_dir, 'output/generate_profiles_and_ini_params');

disp('Step 2 DONE!');

%% Step 3 Generating profiles and initialization parameters
% make sure you have added the CBIG folder into your path
% 由于多回波ICA融合时会进行降噪，这里没有再设置10%的阈值，>0=1, <0=0

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
        TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, targ_mesh, profile_file, 'NONE', '0.1', fMRI_list, ...
            'NONE', 'NONE', split_flag);
        disp([subs{j} ' site-' site{i} ' DONE!' ]);
    else
        disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ]);
    end
 end
end

disp('Step 3 DONE!');

%% Step 4 To obtain the group averaged profiles
% modified from CBIG_MSHBM_avg_profiles.m

if(~exist(fullfile(out_dir,'profiles','avg_profile')))
    mkdir(fullfile(out_dir,'profiles','avg_profile'));
end
num_data = 0;
for j = 1:nsub
    out_profile_dir = fullfile(out_dir,'profiles',['sub' subnums{j}],['sess' site{i}]);       
    avg_profile_file = fullfile(out_dir,'profiles','avg_profile',[targ_mesh '_roi' seed_mesh '_avg_profile.mat']);
    profile_file = fullfile(out_profile_dir,['sub' subnums{j} '_sess' site{i} '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.mat']);      
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
num_initialization = 1000; % !!!!!!!!!!!!!use 1000 for real analysis
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
disp('DONE!');
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
work_dir = fullfile(project_dir, 'output', 'estimate_group_priors');
% Load mat file
final = load(fullfile(work_dir, 'priors', 'Params_Final.mat'));
vertex_num = size(final.Params.s_lambda, 1) / 2;
num_clusters = size(final.Params.s_lambda, 2);
sub_num = size(final.Params.s_lambda, 3);

output_dir = fullfile(work_dir, 'ind_parcellation');
if(~exist(output_dir))
    mkdir(output_dir);
end

% Extract individual parcellation
for i = 1:nsub
    if(isempty(final.Params.s_lambda))
        error('s_lambda is empty. Please use save_all flag when estimating group prior.')
    end
    sub = final.Params.s_lambda(:, :, i);
    [~, labels] = max(sub');
    labels(sum(sub, 2) == 0) = 0; % Label medial wall as 0
    lh_labels = labels(1:vertex_num)';
    rh_labels = labels((vertex_num+1):end)';
    output_file = fullfile(output_dir, ['Ind_parcellation_MSHBM_sub' subnums{i} '.mat']);
    group_file = fullfile(project_dir, 'group', 'group.mat');
    if(exist(group_file, 'file'))
        group = load(group_file);
        if(isfield(group, 'colors')) % Use group color table for individual parcellations
            colors = group.colors;
            save(output_file, 'lh_labels', 'rh_labels', 'colors', 'num_clusters');
        else
            save(output_file, 'lh_labels', 'rh_labels', 'num_clusters');
        end
    else
        save(output_file, 'lh_labels', 'rh_labels', 'num_clusters');
    end
end

%% Step 10 可视化
if ~exist(vis_dir, 'dir')
    mkdir(vis_dir);
end
data_dir = fullfile(project_dir, 'output', 'estimate_group_priors', 'ind_parcellation');
for s = 1:nsub
    % 导入个体分区文件
    indi_data = load(fullfile(data_dir, ['Ind_parcellation_MSHBM_sub', subnums{s}, '.mat']));
    % 导入皮层模版nii文件
    temp = cifti_read(fullfile(code_dir, 'DU15NET_consensus_fsLR_32k.dlabel.nii'));
    temp.cdata = [indi_data.lh_labels;indi_data.rh_labels];

    fn = fullfile(vis_dir, ['sub', subnums{s}, '_Ind_Net_', num2str(num_clusters), '.dlabel.nii']); 
    cifti_write(temp, fn);
end
disp("Finish!")