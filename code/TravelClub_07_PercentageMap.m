%% Pipeline of individual networks mapping
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2025-06-03 @ Beijing Normal University

%% Step 14 Percentage maps for each network

% n=1;t=1;
for n = 1:xnet
    n_str = num2str(networks(n));
    numNetworks = networks(n);  % 获取当前网络数量
    
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');
        matFiles = dir(fullfile(output_dir, 'sub-*.mat'));
        
        % 加载第一个被试的数据以获取维度信息
        sample_data = load(fullfile(output_dir, matFiles(1).name));
        numLHVoxels = length(sample_data.lh_labels);
        numRHVoxels = length(sample_data.rh_labels);

        lh_counts = zeros(numLHVoxels, numNetworks);
        rh_counts = zeros(numRHVoxels, numNetworks);
        
        % 遍历所有被试，统计每个网络在每个体素的出现次数
        for i = 1:nsub
            data = load(fullfile(output_dir, matFiles(i).name));

            for v = 1:numLHVoxels % 左半球统计
                network_idx = data.lh_labels(v);
                if network_idx >= 1 && network_idx <= numNetworks
                    lh_counts(v, network_idx) = lh_counts(v, network_idx) + 1;
                end
            end
           
            for v = 1:numRHVoxels % 右半球统计
                network_idx = data.rh_labels(v);
                if network_idx >= 1 && network_idx <= numNetworks
                    rh_counts(v, network_idx) = rh_counts(v, network_idx) + 1;
                end
            end
        end
        
        lh_percentage = (lh_counts / nsub) * 100;
        rh_percentage = (rh_counts / nsub) * 100;
        
        prob_file = fullfile(output_dir, 'Percentage.mat');
        save(prob_file, 'lh_percentage', 'rh_percentage');
        
    end
end

disp("Step 14 DONE!");

%% Step 15 Visualize probability maps for each network
for n = 1:xnet
    n_str = num2str(networks(n));
    numNetworks = networks(n);
    
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        vis_dir = fullfile(data_dir, 'visual', ['t-' t_str '_n-' n_str]);
        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');
        
        load(fullfile(output_dir, 'Percentage.mat'));

        temp = cifti_read(fullfile(temp_dir, "DU15NET_AgreeProb_fsLR_32k.dscalar.nii"));
        
        for net = 1:numNetworks

            lh_net_prob = lh_percentage(:, net);
            rh_net_prob = rh_percentage(:, net);
            
            net_temp = temp;
            net_temp.cdata = [lh_net_prob; rh_net_prob];
            
            output_file = fullfile(vis_dir, sprintf('Percentage_Network%d.dscalar.nii', net));
            cifti_write(net_temp, output_file);
        end       
    end
end

disp("Step 15 DONE!");