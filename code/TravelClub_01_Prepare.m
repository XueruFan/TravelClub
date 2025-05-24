%% Pipeline of individual networks mapping
% Step 0~4
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 0 Set up parameters

clear;clc

% 启用并行处理
% delete(gcp('nocreate'));
% parpool('Processes', 56, 'IdleTimeout', Inf); % 第二个参数设置worker数量

%—————————— For Xueru's MAC
setenv('CBIG_CODE_DIR', '/Users/xuerufan/matlab-toolbox/CBIG');
setenv('FREESURFER_HOME', '/Applications/freesurfer');
data_dir = '/Volumes/FanXueru/TravelClub';
code_dir = '/Users/xuerufan/TravelClub/code'; 
temp_dir = '/Users/xuerufan/TravelClub/templet';

addpath('/Applications/freesurfer');
addpath('/Applications/workbench/bin_macosx64/');
addpath(genpath('/Users/xuerufan/matlab-toolbox/CBIG'));
addpath(genpath(data_dir));

%—————————— For Linux @205
% setenv('CBIG_CODE_DIR', '/home/ubuntu/homes/FanXueru/toolbox-matlab/CBIG');
% setenv('FREESURFER_HOME', '/home/ubuntu/Softwares/FREESURFER/freesurfer');
% data_dir = '/media/ubuntu/Zuolab_Xueru/TravelClub';
% code_dir = '/media/ubuntu/Zuolab_Xueru/TravelClub/code'; 
% temp_dir = '/media/ubuntu/Zuolab_Xueru/TravelClub/templet';
% 
% addpath('/home/ubuntu/Softwares/FREESURFER/freesurfer');
% addpath('/home/ubuntu/Softwares/workbench/bin_linux64/');
% addpath(genpath('/home/ubuntu/homes/FanXueru/toolbox-matlab/CBIG'));
% addpath(genpath(data_dir));

site = {'A', 'B', 'C', 'D', 'E', 'G'}; 
subid = [2 4 10]; % 注意只有一个site扫描的人去掉，不用于脑图谱绘制的分析
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
% 生成站点被试的fMRI数据地址list文件

dataFolder = fullfile(data_dir, 'data');
outputFolder = fullfile(data_dir, 'output', 'generate_profiles_and_ini_params', 'data_list', 'fMRI_list');

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


% %------------------------------ 并行处理使用下面这部分代码 ------------------------------%
% flatIndices = [];
% for i = 1:nsite
%     for j = 1:nsub
%         flatIndices = [flatIndices; i, j];
%     end
% end
% 
% parfor idx = 1:size(flatIndices, 1)
% 
%     i = flatIndices(idx, 1); % 站点索引
%     j = flatIndices(idx, 2); % 被试索引
% 
%     try
%         siteFolderName = siteFolders(i).name;
%         siteLetter = siteFolderName(end);
%         siteFolderPath = fullfile(dataFolder, siteFolderName);
%         
%         subFolderName = subs{j};
%         subfolderPath = fullfile(siteFolderPath, subFolderName);
%                 
%         if exist(subfolderPath, 'dir')
%             niiFilePaths = {};
%             runFolders = dir(fullfile(subfolderPath, '*PA_run*'));
% 
%             for k = 1:length(runFolders)
%                 runFolderName = runFolders(k).name;
%                 runFolderPath = fullfile(subfolderPath, runFolderName);
%                 niiFile = dir(fullfile(runFolderPath, '*.dtseries.nii'));
% 
%                 if isempty(niiFile)
%                     error('This sub doesnt have this data. Please check!');
%                 end
% 
%                 niiFilePaths{end+1} = fullfile(runFolderPath, niiFile(1).name);
%             end
% 
%             if ~isempty(niiFilePaths)
%                 txtFileName = fullfile(outputFolder, sprintf('%s_site-%s.txt', subFolderName, siteLetter));
%                 fileID = fopen(txtFileName, 'w');
%                 for n = 1:length(niiFilePaths)
%                     if n == length(niiFilePaths)
%                         fprintf(fileID, '%s', niiFilePaths{n});
%                     else
%                         fprintf(fileID, '%s ', niiFilePaths{n});
%                     end
%                 end
%                 fclose(fileID);
%             end
%         end
%     catch ME
%         fprintf('Error processing subject %s in site %s: %s\n', siteFolderName, ME.message);
%     end
% end

%------------------------------ 单机处理使用下面这部分代码 ------------------------------%
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
            runFolders = dir(fullfile(subfolderPath, '*PA_run*'));

            for k = 1:length(runFolders)
%                 k = 1
                runFolderName = runFolders(k).name;
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
% 根据不同的阈值生成功能连接矩阵，由于多回波ICA融合时会进行降噪，在step0里需要小心设置阈值
% modified from CBIG_ComputeCorrelationProfile.m

out_dir = fullfile(data_dir, 'output', 'generate_profiles_and_ini_params');

% %------------------------------ 并行处理使用下面这部分代码 ------------------------------%
% flatIndices = [];
% for j = 1:nsub
%     for i = 1:nsite
%         flatIndices = [flatIndices; j, i];
%     end
% end
% 
% parfor idx = 1:size(flatIndices, 1)
% 
%     j = flatIndices(idx, 1); % 被试索引
%     i = flatIndices(idx, 2); % 站点索引
% 
%     fMRI_list = fullfile(out_dir,'data_list','fMRI_list',[subs{j} '_site-' site{i} '.txt']);
%     
%     if exist(fMRI_list, 'file')
% 
%         for t = 1:xthe 
%             try
%                 out_profile_dir = fullfile(out_dir,'profiles', subs{j}, ['site-' site{i}], ['t-' num2str(threshold(t))]);
%         
%                 if ~exist(out_profile_dir, 'dir')
%                     mkdir(out_profile_dir);
%                 end
%         
% %                 profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.mat']);
%                 profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_t-' num2str(threshold(t)) '.mat']);
% 
%                 TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, targ_mesh, profile_file, 'NONE', num2str(threshold(t)), ...
%                    fMRI_list, 'NONE', 'NONE', '0');
%                 %disp(['Threshold ', num2str(threshold(t)), ' ', subs{j} ' site-' site{i} ' DONE!' ]);
%                 
%             catch ME
%                 fprintf('Error processing subject %s in site %s: %s\n', subs{j}, site{i}, ME.message);
%             end
%         end
%     else
%         disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ]);
%     end
% end

%------------------------------ 单机处理使用下面这部分代码 ------------------------------%
for j = 1:nsub
    for i = 1:nsite
%      i = 1;j = 1;

        fMRI_list = fullfile(out_dir,'data_list','fMRI_list',[subs{j} '_site-' site{i} '.txt']);
        
        if exist(fMRI_list, 'file')
   
            for t = 1:xthe 
%                 t = 1;
                out_profile_dir = fullfile(out_dir,'profiles', subs{j}, ['site-' site{i}], ['t-' num2str(threshold(t))]);
        
                if ~exist(out_profile_dir, 'dir')
                    mkdir(out_profile_dir);
                end
        
                profile_file = fullfile(out_profile_dir, [subs{j} '_site-' site{i} '_t-' num2str(threshold(t)) '.mat']);

                TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, targ_mesh, profile_file, 'NONE', num2str(threshold(t)), ...
                   fMRI_list, 'NONE', 'NONE', '0');
            end
            
        else
            disp([subs{j} ' site-' site{i} ' DOES NOT EXIT' ]);
        end
    end
end

disp('Step 2 DONE!');

%% Step 3 把每个被试的功能连接结果路径按照要求整理成txt文件 Prepare the individual FC list files
% 分阈值分别生成list文件，这部分处理简单，单机就可以
% 这部分具体的整理要求参阅CBIG_MSHBM_estimate_group_priors.m里的描述

profilesFolder = fullfile(data_dir, 'output', 'generate_profiles_and_ini_params', 'profiles');

for t = 1:xthe
    t_str = num2str(threshold(t));
    for n = 1:xnet
        n_str = num2str(networks(n));

        outputFolder = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'estimate_group_priors', ...
        'profile_list', 'training_set');
    
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
    
        for i = 1:nsite
            siteFolderName = ['site-' site{i}];
         
            matFilePaths = {};
        
            txtFileName = fullfile(outputFolder, [siteFolderName '.txt']);
        
            for j = 1:nsub               
                subFolderName = subs{j};
                subfolderPath = fullfile(profilesFolder, subFolderName, siteFolderName, ['t-' t_str]);
            
                matFiles = dir(fullfile(subfolderPath, '*.mat'));
                matFiles = matFiles(~startsWith({matFiles.name}, '._'));
    
                if isempty(matFiles)
                    matFilePaths{end+1}  = 'NONE';
                else
                    matFilePaths{end+1}  = fullfile(subfolderPath, matFiles(1).name);
                end
        
                fileID = fopen(txtFileName, 'w');
                for m = 1:length(matFilePaths)
                    fprintf(fileID, '%s\n', matFilePaths{m});
                end
                fclose(fileID);
        
            end
        end 
    end
end

disp('Step 3 DONE!');

%% Step 4 To obtain the group averaged profiles
% 分阈值计算组水平功能连接，这部分处理简单，单机就可以
% modified from CBIG_MSHBM_avg_profiles.m

out_dir = fullfile(data_dir, 'output', 'generate_profiles_and_ini_params');

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

disp('Step 4 DONE!');