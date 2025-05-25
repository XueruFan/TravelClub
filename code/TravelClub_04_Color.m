%% Pipeline of individual networks mapping
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2025-05-24 @ Home

% Step 9 修改网络名称和颜色 Modify name and color for networks
% 每个阈值、网络的都不一样，全部需要手动修改

%% 15 networks, Threshold 0.1
n_str = '15'; t_str = '0.1';

vis_dir = fullfile(data_dir, 'visual', ['t-' t_str '_n-' n_str]);

for j = 1:nsub

    fn = fullfile(vis_dir, [subs{j}, '_t-', t_str, '_n-', n_str, '.dlabel.nii']);
    temp = cifti_read(fn);

    % 修改元数据
    temp.metadata = [];

    % 自定义每个网络名称和颜色对应的标签（原文件里的颜色和网络名称是对应的）
    idx = find(strcmp({temp.diminfo{1, 2}.maps.table.name}, 'FPN-A'));
    temp.diminfo{1, 2}.maps.table(idx).key = 2;

    idx = find(strcmp({temp.diminfo{1, 2}.maps.table.name}, 'CG-OP'));
    temp.diminfo{1, 2}.maps.table(idx).key = 11;
    
    %%%%%%%%%%%%%%%%%%%% 这里依次改好所有网络

    % 把table里的key按照顺序排放
    table_data = temp.diminfo{1, 2}.maps.table;
    keys = [table_data.key];
    [sorted_keys, sort_idx] = sort(keys);
    sorted_table = table_data(sort_idx);
    temp.diminfo{1, 2}.maps.table = sorted_table;

    cifti_write(temp, fn);
end

disp("Step 9 DONE!")
