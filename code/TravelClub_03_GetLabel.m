%% Pipeline of individual networks mapping
% Step 7~8
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 7 提取出个体分区的结果label Extract individual network labels
% 这部分处理简单，不需要并行，循环逐个提取出来即可
% modified from CBIG_IndCBM_extract_MSHBM_result.m

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));

        work_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'estimate_group_priors');
        final = load(fullfile(work_dir, 'priors', 'Params_Final.mat'));
        vertex_num = size(final.Params.s_lambda, 1) / 2;
        
        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');
        if(~exist(output_dir))
            mkdir(output_dir);
        end
        
        if(isempty(final.Params.s_lambda))
                error('s_lambda is empty. Please use save_all flag when estimating group prior.')
        end
        
        for j = 1:nsub
            sub = [];
            labels = [];
            lh_labels = [];
            rh_labels = [];
    
            sub = final.Params.s_lambda(:, :, j);
        
            [~, labels] = max(sub');
            labels(sum(sub, 2) == 0) = 0; % Label medial wall as 0
            lh_labels = labels(1:vertex_num)';
            rh_labels = labels((vertex_num + 1):end)';
        
            output_file = fullfile(output_dir, [subs{j} '.mat']);
            save(output_file, 'lh_labels', 'rh_labels');
        end
    end
end

disp('Step 7 DONE!');

%% Step 8 可视化 Make individual visualizations
% 同不需要并行

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));

        vis_dir = fullfile(data_dir, 'visual', ['t-' t_str '_n-' n_str]);

        if ~exist(vis_dir, 'dir')
            mkdir(vis_dir);
        end

        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');

        for j = 1:nsub

            lh_labels = [];
            rh_labels = [];
            temp = [];
            fn = '';
    
            load(fullfile(output_dir, [subs{j}, '.mat']));
        
            temp = cifti_read(fullfile(temp_dir, 'DU15NET_consensus_fsLR_32k.dlabel.nii'));
            temp.cdata = [lh_labels; rh_labels];
        
            fn = fullfile(vis_dir, [subs{j}, '_t-', t_str, '_n-', n_str, '.dlabel.nii']); 
            cifti_write(temp, fn);
        end
    end
end

disp("Step 8 DONE!")
disp("-----------------  ALL DONE !  -----------------")
