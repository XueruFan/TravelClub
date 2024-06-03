function TravelClub_MSHBM_VisNet(data_dir, vis_dir, num_sub, num_clusters)
% To visualize the individual networks results
% Input:
%   - data_dir: where is the individual network mat files
%   - vis_dir: where to save the results. remember you neet to copy the
%              templet file within this folder
%   - num_sub: how many subjects
%   - num_clusters: how many networks
% Written by Xue-Ru Fan 2024-06-03 @ Beijing Normal University
% Email:xueru@mail.bnu.edu.cn

for s = 1:str2double(num_sub)
    % 导入个体分区文件
    indi_data = load(fullfile(data_dir, ['Ind_parcellation_MSHBM_sub', num2str(s), '.mat']));
    % 导入皮层模版nii文件
    temp = cifti_read(fullfile(vis_dir, 'templet', 'DU15NET_consensus_fsLR_32k.dlabel.nii'));
    temp.cdata = [indi_data.lh_labels;indi_data.rh_labels];

    fn = fullfile(vis_dir, ['sub', num2str(s), '_Ind_Net_', num_clusters, '.dlabel.nii']); 
    cifti_write(temp, fn);
end
disp("Finish!")
end