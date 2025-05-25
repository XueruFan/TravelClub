%% Pipeline of individual networks mapping
% Step 10
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University

%% Step 10 Consensus maps according to winner-take-all strategy

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');
        
        log_file = fullfile(output_dir, 'NetowrkAssignmentLog.txt');
        fid_log = fopen(log_file, 'w');
        fprintf(fid_log, '处理日志: 网络%s 阈值%s\n生成时间: %s\n\n', n_str, t_str, datestr(now));
        fclose(fid_log);

        logmessage = @(msg) fprintf(fopen(log_file, 'a'), '%s: %s\n', datestr(now), msg);
        
        matFiles = dir(fullfile(output_dir, 'sub-*.mat'));
        numSubjects = length(matFiles);
        
        sample_data = load(fullfile(output_dir, matFiles(1).name));
        numVoxels = length(sample_data.lh_labels);
        lh_all_labels = zeros(numVoxels, numSubjects);
        rh_all_labels = zeros(numVoxels, numSubjects);
              
        for i = 1:numSubjects
            data = load(fullfile(output_dir, matFiles(i).name));
            lh_all_labels(:, i) = data.lh_labels;
            rh_all_labels(:, i) = data.rh_labels;
        end
        
        tie_msg = sprintf('平局: lh: %d处, rh: %d处', lh_ties, rh_ties);
        log_message(log_file, tie_msg);

        [lh_consensus, lh_ties] = process_hemisphere(lh_all_labels, 'lh', log_file); 
        [rh_consensus, rh_ties] = process_hemisphere(rh_all_labels, 'rh', log_file);
        
        output_file = fullfile(output_dir, 'Consensus.mat');
        save(output_file, 'lh_consensus', 'rh_consensus');

    end
end

disp("Step 10 DONE!")

%% Step 11 Make consensus map visualizations

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
%         n_str = num2str(networks(1)); t_str = num2str(threshold(1));

        vis_dir = fullfile(data_dir, 'visual', ['t-' t_str '_n-' n_str]);
        load(fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation', 'Consensus.mat'));
    
        temp = cifti_read(fullfile(vis_dir, [subs{1}, '_t-', t_str, '_n-', n_str, '.dlabel.nii']));
        temp.cdata = [lh_consensus; rh_consensus];

        cifti_write(temp, fullfile(vis_dir, 'Consensus.dlabel.nii'));

    end
end

disp("Step 11 DONE!")


%% 平局处理函数（带日志记录）
function [consensus_labels, tie_count] = process_hemisphere(all_labels, hemisphere_name, log_file)
    [numVoxels, ~] = size(all_labels);
    consensus_labels = zeros(numVoxels, 1);
    tie_count = 0;
    
    all_possible_labels = unique(all_labels(:));
    
    for v = 1:numVoxels
        labels = all_labels(v, :);
        counts = histcounts(labels, [all_possible_labels; max(all_possible_labels)+1]);
        max_count = max(counts);
        candidates = all_possible_labels(counts == max_count);
        
        if length(candidates) == 1
            consensus_labels(v) = candidates;
        else
            selected_label = candidates(randi(length(candidates))); % 随机选择一个平局标签
            consensus_labels(v) = selected_label;
            
             log_message(log_file, sprintf(...
                 '%s %05d 候选标签: %s (次数=%d) 随机选择: %d', ...
                hemisphere_name, v, strjoin(arrayfun(@num2str, candidates(:), 'UniformOutput', false), ', '), ...
                max_count, selected_label));
            
            tie_count = tie_count + 1;
        end
    end
end

%% 日志记录辅助函数
function log_message(log_file, message)
    fid = fopen(log_file, 'a');
    if fid == -1
        error('无法打开日志文件: %s', log_file);
    end
    fprintf(fid, '%s\n', message);  % 移除了时间部分
    fclose(fid);
end