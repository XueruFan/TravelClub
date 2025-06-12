%% Pipeline of individual networks mapping
% 3R-BRAIN Project
% Modified with majority consensus condition
% Written by Xue-Ru Fan 2025-05-23 @ Beijing Normal University

%% Step 10 Consensus maps with majority condition

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');
        
        log_file = fullfile(output_dir, 'NetworkAssignmentLog.txt');
        fid_log = fopen(log_file, 'w');
        fprintf(fid_log, '处理日志: 网络%s 阈值%s\n生成时间: %s\n\n', n_str, t_str, datestr(now));
        fclose(fid_log);

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
        
        % 处理左右半球（新增numSubjects参数）
        [lh_consensus, lh_stats] = process_hemisphere(lh_all_labels, numSubjects, 'lh', log_file); 
        [rh_consensus, rh_stats] = process_hemisphere(rh_all_labels, numSubjects, 'rh', log_file);
        
        % 记录统计信息
        log_message(log_file, sprintf('左半球统计: 未分配体素=%d (%.1f%%), 平局=%d', ...
            lh_stats.unassigned, lh_stats.unassigned_pct, lh_stats.ties));
        
        log_message(log_file, sprintf('右半球统计: 未分配体素=%d (%.1f%%), 平局=%d', ...
            rh_stats.unassigned, rh_stats.unassigned_pct, rh_stats.ties));

        output_file = fullfile(output_dir, 'Consensus_withThreshold.mat');
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
        load(fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation', 'Consensus_withThreshold.mat'));
    
        temp = cifti_read(fullfile(vis_dir, [subs{1}, '_t-', t_str, '_n-', n_str, '.dlabel.nii']));
        temp.cdata = [lh_consensus; rh_consensus];

        cifti_write(temp, fullfile(vis_dir, 'Consensus_withThreshold.dlabel.nii'));

    end
end

disp("Step 11 DONE!")

%% 新版处理函数：带多数共识条件
function [consensus_labels, stats] = process_hemisphere(all_labels, numSubjects, hemisphere_name, log_file)
    [numVoxels, ~] = size(all_labels);
    consensus_labels = zeros(numVoxels, 1);
    
    % 初始化统计信息
    stats = struct();
    stats.unassigned = 0;    % 未分配网络的体素数
    stats.ties = 0;          % 满足多数条件但出现平局的体素数
    
    for v = 1:numVoxels
        labels = all_labels(v, :);
        [unique_labels, ~, ic] = unique(labels);
        counts = accumarray(ic, 1);
        
        [max_count, max_idx] = max(counts);
        majority_threshold = ceil(numSubjects / 2);  % 超过一半的最小整数
        
        % 检查是否满足多数条件
        if max_count >= majority_threshold
            candidates = unique_labels(counts == max_count);
            
            if numel(candidates) == 1
                consensus_labels(v) = candidates;
            else
                % 满足多数条件但出现平局：随机选择
                selected_label = candidates(randi(numel(candidates)));
                consensus_labels(v) = selected_label;
                stats.ties = stats.ties + 1;
                
                % 记录平局信息
                log_message(log_file, sprintf('%s 体素 %d: 平局 (%.0f%%), 候选网络: %s, 选择: %d', ...
                    hemisphere_name, v, 100*max_count/numSubjects, ...
                    strjoin(cellstr(num2str(candidates(:))), ', '), ...
                    selected_label));
            end
        else
            % 未满足多数条件：标记为0（未分配）
            consensus_labels(v) = 0;
            stats.unassigned = stats.unassigned + 1;
        end
    end
    
    % 计算未分配体素的百分比
    stats.unassigned_pct = 100 * stats.unassigned / numVoxels;
end

%% 日志记录辅助函数
function log_message(log_file, message)
    fid = fopen(log_file, 'a');
    if fid == -1
        error('无法打开日志文件: %s', log_file);
    end
    fprintf(fid, '%s: %s\n', datestr(now), message);
    fclose(fid);
end