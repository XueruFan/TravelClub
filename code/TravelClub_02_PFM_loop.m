%% Pipeline of individual networks mapping
% 如果使用循环依次处理 不同阈值、不同网络个数 使用本代码
% Step 5~6
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 5 To run Yeo2011 clustering algorithm for generate our own group prior
% modified from CBIG_MSHBM_generate_ini_params.m

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        out_dir = fullfile(data_dir, 'output', 'generate_profiles_and_ini_params');
        gro_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'estimate_group_priors', 'group');
        
        if(~exist(fullfile(gro_dir)))
            mkdir(fullfile(gro_dir));
        end
        
        avg_profile_dir = fullfile(out_dir, 'profiles', 'avg_profile');
        avg_profile_file = fullfile(avg_profile_dir, ['t-' t_str '_avg_profile.mat']);
        output_file = fullfile(gro_dir, 'group.mat');
        
        CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, '', n_str, output_file, ...
            avg_profile_file, 'NONE', 0, niter, 0, 1000, 1);
        
        load(output_file);
        clustered.mtc = mtc;
        clustered.lowerbound = lowerbound;
        clustered.lambda = lambda;
        save(output_file,'lh_labels','rh_labels','clustered');
    end
end

disp('Step 5 DONE!');

%% Step 6 估计模型的参数同时得出个体分区 estimate individual mapping

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        work_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'estimate_group_priors');
        
        % estimate group priors 这里要用雪如修改过的代码
        Params = TravelClub_CBIG_MSHBM_estimate_group_priors(work_dir, targ_mesh, num2str(nsub), num2str(nsite), ...
            site, n_str); % Du 设置了 max_iter 5
    end
end

% mex -v -largeArrayDims mtimesx.c -lmwblas -lmwlapack 如果需要编译
% Add the result directory and its subfolders to the MATLAB path
addpath(genpath(data_dir));
savepath

disp('Step 6 DONE!');
