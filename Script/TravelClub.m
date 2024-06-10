%% Pipeline of Xue-Ru's anaysis
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University
% Email:xueru@mail.bnu.edu.cn
%% Section 1 Mapping individual networks
    %% Step 0
    clear;clc
    setenv('CBIG_CODE_DIR', 'C:\CBIG');
    mex -setup C
    project_dir = 'E:\PhDproject\HCP\HCP_test4MSHBM';
    code_dir = 'C:\TravelClub\Script';
   
    % addpath
    addpath(genpath("C:\CBIG"));
    addpath(genpath("C:\matlab-toolbox"))
    addpath(genpath("C:\Program Files\MATLAB\R2022b\toolbox\"))
    savepath
    %% Step 1
    TravelClub_MSHBM_Prep(project_dir)
    %% Step 2
    % see Kong2019_MSHBM README.md for input suggestions
%     TravelClub_MSHBM_SetInput
    seed_mesh = 'fs_LR_900';
    targ_mesh = 'fs_LR_32k';
    nsub = '2';
    maxsess = '2';
    split_flag = '0';
    project_dir = 'E:\PhDproject\HCP\HCP_test4MSHBM\output\generate_profiles_and_ini_params';
    %% Step 3
    % Generating profiles and initialization parameters
    % make sure you have added the CBIG folder into your path
    for sub = 1:str2num(nsub)
     for sess = 1:str2num(maxsess)
        CBIG_MSHBM_generate_profiles(seed_mesh, targ_mesh, project_dir, ...
            num2str(sub), num2str(sess), '0');
     end
    end
    %% Step 4
    % To obtain the group averaged profiles
    CBIG_MSHBM_avg_profiles(seed_mesh, targ_mesh, project_dir, nsub, maxsess);
    %% Step 5
    % To run Yeo2011 clustering algorithm for generate or own group prior
    % 这里用的就是自己的数据
    num_clusters = 15;
    num_initialization = 2; % !!!!!!!!!!!!!use 1000 for real analysis
    CBIG_MSHBM_generate_ini_params(seed_mesh, targ_mesh, num_clusters, ...
        num_initialization, project_dir)
    % 这一步Du用HCP的40个人的数据，从2~21个网络都run了一遍，各做了一个组水平的
    % 先验，最终是结合后面基于ROI的连接图，确定了15个网络
    %% Step 6
    % 把每个被试的功能连接结果路径按照要求整理成txt文件
    % 这部分具体的真理要求参阅CBIG_MSHBM_estimate_group_priors.m里的描述
    TravelClub_MSHBM_ProfilePrep(project_dir);
    %% Step 7
    % 把生成的group.mat文件和文件夹拷贝到MSHBM_Step2的文件夹里
    TravelClub_MSHBM_CopyGroup(project_dir)
    %% Step 8 估计模型的参数同时得出个体分区
    % set up input
    work_dir = fullfile(project_dir, 'output', 'estimate_group_priors');
    mesh = 'fs_LR_32k';
    num_sub = '3';
    num_sessions = '2';
    num_clusters = '15';
    % estimate group priors
    Params = CBIG_MSHBM_estimate_group_priors(work_dir, mesh, num_sub, ...
        num_sessions, num_clusters);
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