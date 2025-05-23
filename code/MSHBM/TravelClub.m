%% Pipeline of individual networks mapping
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 0 Set up parameters

clear;clc

% 启用并行处理
% delete(gcp('nocreate'));
% parpool('Processes', 56, 'IdleTimeout', Inf);

%—————————— MAC
setenv('CBIG_CODE_DIR', '/Users/xuerufan/matlab-toolbox/CBIG');
setenv('FREESURFER_HOME', '/Applications/freesurfer');
project_dir = '/Volumes/Zuolab_Xueru/TravelClub';
code_dir = '/Volumes/Zuolab_Xueru/TravelClub/code/MSHBM'; 
temp_dir = '/Volumes/Zuolab_Xueru/TravelClub/templet';

addpath('/Applications/freesurfer');
addpath('/Applications/workbench/bin_macosx64/');
addpath(genpath('/home/ubuntu/homes/FanXueru/toolbox-matlab/CBIG'));
addpath(genpath('/media/ubuntu/Zuolab_Xueru/TravelClub'));

%—————————— Linux
% setenv('CBIG_CODE_DIR', '/home/ubuntu/homes/FanXueru/toolbox-matlab/CBIG');
% setenv('FREESURFER_HOME', '/home/ubuntu/Softwares/FREESURFER/freesurfer');
% project_dir = '/media/ubuntu/Zuolab_Xueru/TravelClub';
% code_dir = '/media/ubuntu/Zuolab_Xueru/TravelClub/code/MSHBM'; 
% temp_dir = '/media/ubuntu/Zuolab_Xueru/TravelClub/templet';
% 
% addpath('/home/ubuntu/Softwares/FREESURFER/freesurfer');
% addpath('/home/ubuntu/Softwares/workbench/bin_linux64/');
% addpath(genpath('/home/ubuntu/homes/FanXueru/toolbox-matlab/CBIG'));
% addpath(genpath('/media/ubuntu/Zuolab_Xueru/TravelClub'));

site = {'A', 'B', 'C', 'D', 'E', 'G'};
subid = [2 4 10]; % 注意只有一个site扫描的人去掉
% subid = [1:16 18:24 26:28 30 32 34:39 42];
seed_mesh = 'fs_LR_900';
targ_mesh = 'fs_LR_32k';
networks = [15 10];%[7:14 16];
threshold = [10 20]/100; %0.1; %参考du的方案（也就是默认方案)
niter = 4; % 正式改成1000


mex -setup C
mex -setup C++
nsite = length(site);
subs = arrayfun(@(x) {sprintf('sub-%03d', x)}, subid);
% subs = arrayfun(@(x) {sprintf('%d', x)}, subid); % test
nsub = length(subid);
xnet = length(networks);
xthe = length(threshold);
cd(code_dir)


disp('Step 0 DONE!');

%% Step 1 prepare the fMRI list files

dataFolder = fullfile(project_dir, 'data');
outputFolder = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params', 'data_list', 'fMRI_list');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

siteFolders = dir(dataFolder);
siteFolders = siteFolders([siteFolders.isdir]);
siteFolders = siteFolders(~ismember({siteFolders.name}, {'.', '..'}));

% 测试阶段不需要下面这部分
% if length(siteFolders) ~= nsite
%     error('Number of site folders doesnt match the expected number of sites. Please check!');
% end



%---- 并行处理
% 由于parfor里不能进行内层并行，这里创建一个扁平化的索引，然后对这个扁平化的循环使用 parfor
flatIndices = [];
for i = 1:nsite
    for j = 1:nsub
        flatIndices = [flatIndices; i, j];
    end
end

parfor idx = 1:size(flatIndices, 1)

    i = flatIndices(idx, 1); % 站点索引
    j = flatIndices(idx, 2); % 被试索引

    try
        siteFolderName = siteFolders(i).name;
        siteLetter = siteFolderName(end);
        siteFolderPath = fullfile(dataFolder, siteFolderName);
        
        subFolderName = subs{j};
        subfolderPath = fullfile(siteFolderPath, subFolderName);
                
        if exist(subfolderPath, 'dir')
            niiFilePaths = {};
            runFolders = dir(fullfile(subfolderPath, '*PA_run*'));

            for k = 1:length(runFolders)
                runFolderName = runFolders(k).name;
                runFolderPath = fullfile(subfolderPath, runFolderName);
                niiFile = dir(fullfile(runFolderPath, '*.dtseries.nii'));

                if isempty(niiFile)
                    error('This sub doesnt have this data. Please check!');
                end

                niiFilePaths{end+1} = fullfile(runFolderPath, niiFile(1).name);
            end

            if ~isempty(niiFilePaths)
                txtFileName = fullfile(outputFolder, sprintf('%s_site-%s.txt', subFolderName, siteLetter));
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
    catch ME
        fprintf('Error processing subject %s in site %s: %s\n', siteFolderName, ME.message);
    end
end

%---- 单机处理
for i = 1:nsite
%     i = 1
    siteFolderName = siteFolders(i).name;
    siteLetter = siteFolderName(end);
    siteFolderPath = fullfile(dataFolder, siteFolderName);
    
    for j = 1:nsub
%         j = 1
        subFolderName = subs{j};
        subfolderPath = fullfile(siteFolderPath, subFolderName);
        
        if exist(subfolderPath, 'dir')
            niiFilePaths = {};

%             runFolders = dir(fullfile(subfolderPath, 'MNINonLinear', 'Results', '*PA_run*'));
            runFolders = dir(fullfile(subfolderPath, '*PA_run*'));

            for k = 1:length(runFolders)
%                 k = 1
                runFolderName = runFolders(k).name;
%                 runFolderPath = fullfile(subfolderPath, 'MNINonLinear', 'Results', ...
%                     runFolderName);
                runFolderPath = fullfile(subfolderPath, runFolderName);
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
% 由于多回波ICA融合时会进行降噪，在step0里需要小心设置阈值

out_dir = fullfile(project_dir, 'output', 'generate_profiles_and_ini_params');

flatIndices = [];
for j = 1:nsub
    for i = 1:nsite
        flatIndices = [flatIndices; j, i];
    end
end

parfor idx = 1:size(flatIndices, 1)

    j = flatIndices(idx, 1); % 被试索引
    i = flatIndices(idx, 2); % 站点索引

    fMRI_list = fullfile(out_dir,'data_list','fMRI_list',[subs{j} '_site-' site{i} '.txt']);
    
    if exist(fMRI_list, 'file')

        for t = 1:xthe 
            try
                out_profile_dir = fullfile(out_dir,'profiles', subs{j}, ['site-' site{i}], ['t-' num2str(threshold(t))]);
        
                if ~exist(out_profile_dir, 'dir')
                    mkdir(out_profile_dir);
                end
        
%                 profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.mat']);
                profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_t-' num2str(threshold(t)) '.mat']);

                TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, targ_mesh, profile_file, 'NONE', num2str(threshold(t)), ...
                   fMRI_list, 'NONE', 'NONE', '0');
                %disp(['Threshold ', num2str(threshold(t)), ' ', subs{j} ' site-' site{i} ' DONE!' ]);
                
            catch ME
                fprintf('Error processing subject %s in site %s: %s\n', subs{j}, site{i}, ME.message);
            end
        end
    else
        disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ]);
    end
end


% for j = 1:nsub
%  for i = 1:nsite
% %      i = 1
% %      j = 1
% 
%     out_profile_dir = fullfile(out_dir,'profiles', subs{j}, ['site-' site{i}]);
% 
%     if ~exist(out_profile_dir)
%         mkdir(out_profile_dir);
%     end
% 
%     profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_' targ_mesh ...
%         '_roi' seed_mesh '.surf2surf_profile.mat']);
%     fMRI_list = fullfile(out_dir,'data_list','fMRI_list',[subs{j} '_site-' site{i} '.txt']);
% 
%     if exist(fMRI_list)
%         TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, targ_mesh, profile_file, 'NONE', num2str(threshold), ...
%             fMRI_list, 'NONE', 'NONE', '0');
%         disp([subs{j} ' site-' site{i} ' DONE!' ]);
%     else
%         disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ]);
%     end
%  end
% end


disp('Step 2 DONE!');

%% Step 3 To obtain the group averaged profiles
% modified from CBIG_MSHBM_avg_profiles.m

out_dir = fullfile(project_dir, 'output/generate_profiles_and_ini_params');

avg_profile_dir = fullfile(out_dir, 'profiles', 'avg_profile');

if ~exist(avg_profile_dir)
    mkdir(avg_profile_dir);
end

for t = 1:xthe

    t_str = num2str(threshold(t));
    
    num_data = 0;

    avg_profile_file = fullfile(avg_profile_dir, ['t-' t_str '_avg_profile.mat']);

    for j = 1:nsub
        for i = 1:nsite
    
            out_profile_dir = fullfile(out_dir, 'profiles', subs{j}, ['site-', site{i}], ['t-', t_str]);  

%             profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_' targ_mesh ...
%                 '_roi' seed_mesh '.surf2surf_profile.mat']);
            profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_t-' t_str '.mat']);
        
            if exist(profile_file)
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
end

disp('Step 3 DONE!');

%% Step 4 To run Yeo2011 clustering algorithm for generate our own group prior
% modified from CBIG_MSHBM_generate_ini_params.m

% tic; % 开始计时

out_dir = fullfile(project_dir, 'output/generate_profiles_and_ini_params');

gro_dir = fullfile(project_dir, 'output', 'estimate_group_priors', 'group');

if(~exist(fullfile(gro_dir)))
    mkdir(fullfile(gro_dir));
end

avg_profile_dir = fullfile(out_dir, 'profiles', 'avg_profile');

% 并行工作 还不好使
% parfor n_idx = 1:xnet   
%     n_str = num2str(networks(n));
%     for t_idx = 1:xthe
%         t_str = num2str(threshold(t));
%         
%         TravelClub_CBIG_MSHBM_generate_ini_params(targ_mesh, current_network, ...
%             fullfile(gro_dir, ['t-' t_str '_n-' n_str '_group.mat']), ...
%             fullfile(avg_profile_dir, ['t-' t_str '_avg_profile.mat']), ...
%             niter);
%     end
% end


% 线性工作
avg_profile_dir = fullfile(out_dir, 'profiles', 'avg_profile');

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        output_file = fullfile(gro_dir, ['t-' t_str '_n-' n_str '_group.mat']);

        % avg_profile_file = fullfile(out_dir,'profiles','avg_profile',[targ_mesh '_roi' seed_mesh '_avg_profile.mat']);
        avg_profile_file = fullfile(avg_profile_dir, ['t-' t_str '_avg_profile.mat']);
        
        CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, '', networks(n), output_file, ...
            avg_profile_file, 'NONE', 0, niter, 0, 1000, 1);
        
        load(output_file);
        clustered.mtc = mtc;
        clustered.lowerbound = lowerbound;
        clustered.lambda = lambda;
        save(output_file,'lh_labels','rh_labels','clustered');

    end
end

disp('Step 4 DONE!');

% toc;
% elapsedTime = toc/60;
% disp(['耗时', num2str(elapsedTime), '分钟'])

%% Step 5 把每个被试的功能连接结果路径按照要求整理成txt文件
% 这部分具体的整理要求参阅CBIG_MSHBM_estimate_group_priors.m里的描述

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

        for t = 1:xthe%%%%%%%%%%%%%%%%%%%

            theFolderName = num2str(threshold(t));
            thefolderPath = fullfile(subfolderPath, ['t-' theFolderName]);
        
            matFiles = dir(fullfile(thefolderPath, '*.mat'));
            if isempty(matFiles)
                matFilePaths{end+1}  = 'NONE';
            else
                matFilePaths{end+1}  = fullfile(thefolderPath, matFiles(1).name);
            end
    
            fileID = fopen(txtFileName, 'w');
            for n = 1:length(matFilePaths)
                fprintf(fileID, '%s\n', matFilePaths{n});
            end
            fclose(fileID);

        end
    end  
end

disp('Step 5 DONE!');

%% Step 6 估计模型的参数同时得出个体分区

work_dir = fullfile(project_dir, 'output', 'estimate_group_priors');

% estimate group priors 这里要用雪如修改过的代码
Params = TravelClub_CBIG_MSHBM_estimate_group_priors(work_dir, targ_mesh, num2str(nsub), num2str(nsite), ...
    site, num2str(networks)) % Du 设置了 max_iter 5

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

if(isempty(final.Params.s_lambda))
        error('s_lambda is empty. Please use save_all flag when estimating group prior.')
end

parfor j = 1:nsub
    try
        sub = [];
        labels = [];
        lh_labels = [];
        rh_labels = [];

        sub = final.Params.s_lambda(:, :, j);
    
        [~, labels] = max(sub');
        labels(sum(sub, 2) == 0) = 0; % Label medial wall as 0
        lh_labels = labels(1:vertex_num)';
        rh_labels = labels((vertex_num + 1):end)';
    
        output_file = fullfile(output_dir, ['Ind_parcellation_MSHBM_' subs{j} '.mat']);
        save(output_file, 'lh_labels', 'rh_labels', 'networks');

        catch ME
        fprintf('Error processing subject %s: %s\n', subs{j}, ME.message);
    end
end

disp('Step 7 DONE!');

%% Step 8 可视化

vis_dir = fullfile(project_dir, 'visual');
if ~exist(vis_dir, 'dir')
    mkdir(vis_dir);
end

data_dir = fullfile(project_dir, 'output', 'ind_parcellation');

parfor j = 1:nsub
    try
        % 显式初始化临时变量
        lh_labels = [];
        rh_labels = [];
        temp = [];
        fn = '';

        load(fullfile(data_dir, ['Ind_parcellation_MSHBM_', subs{j}, '.mat']));
    
        temp = cifti_read(fullfile(temp_dir, 'DU15NET_consensus_fsLR_32k.dlabel.nii'));
        temp.cdata = [lh_labels; rh_labels];
    
        fn = fullfile(vis_dir, [subs{j}, '_', num2str(networks), 'Net_HCP.dlabel.nii']); 
        cifti_write(temp, fn);

        catch ME
        fprintf('Error processing subject %s: %s\n', subs{j}, ME.message);
    end
end

disp("Step 8 DONE!")
