
%% 1
clc;clear;close;

data = cifti_read("D:\TravelClub\data\3RB2_A\sub-002\EyeOpen_PA_run1\EyeOpen_PA_run1_Atlas_s0.dtseries.nii");
bold_ts = data.cdata(1:20,1:100);

figure('Position', [100 100 200 200]);
hold on;
for i = 1:20
    plot(bold_ts(i,1:100), 'LineWidth', 0.5);
end

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k'); % 黑色边框
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\BOLD1.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');


fc_matrix = corr(bold_ts');

figure('Position', [100 100 200 200]); % 设置图形大小
imagesc(fc_matrix);
colormap(jet); 
caxis([-1 1]);
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC1.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

% 保留前10%强度的连接并二值化 ---
% 步骤1：将矩阵的下三角部分（不含对角线）展开为向量
mask = tril(true(size(fc_matrix)), -1); % 下三角掩码（排除对角线）
fc_values = fc_matrix(mask);            % 提取下三角的所有连接强度值

% 步骤2：计算前10%强度的阈值
threshold = prctile(fc_values, 90);   % 找到90百分位数（即保留前10%）

% 步骤3：生成二值化矩阵（保留前10%的连接设为1，其余为0）
binary_fc = zeros(size(fc_matrix));
binary_fc(fc_matrix >= threshold & ~eye(size(fc_matrix))) = 1; % 排除对角线

% --- 绘制二值化热图 ---
figure('Position', [100 100 200 200]);
imagesc(binary_fc);
colormap(parula);
caxis([0 1]);
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC1_threshold.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

%% 2
clc;clear;close;

data = cifti_read("D:\TravelClub\data\3RB2_A\sub-002\EyeClose_PA_run2\EyeClose_PA_run2_Atlas_s0.dtseries.nii");
bold_ts = data.cdata(1:20,1:100);

figure('Position', [100 100 200 200]);
hold on;
for i = 1:20
    plot(bold_ts(i,1:100), 'LineWidth', 0.5);
end

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k'); % 黑色边框
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\BOLD2.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');


fc_matrix = corr(bold_ts');

figure('Position', [100 100 200 200]); % 设置图形大小
imagesc(fc_matrix);
colormap(jet); 
caxis([-1 1]);
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC2.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

mask = tril(true(size(fc_matrix)), -1);
fc_values = fc_matrix(mask); 
threshold = prctile(fc_values, 90); 
binary_fc = zeros(size(fc_matrix));
binary_fc(fc_matrix >= threshold & ~eye(size(fc_matrix))) = 1; % 排除对角线

figure('Position', [100 100 200 200]);
imagesc(binary_fc);
colormap(parula);
caxis([0 1]);        
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC2_threshold.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

%% 3
clc;clear;close;

data = cifti_read("D:\TravelClub\data\3RB2_A\sub-002\EyeClose_PA_run3\EyeClose_PA_run3_Atlas_s0.dtseries.nii");
bold_ts = data.cdata(1:20,1:100);

figure('Position', [100 100 200 200]);
hold on;
for i = 1:20
    plot(bold_ts(i,1:100), 'LineWidth', 0.5);
end

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k'); % 黑色边框
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\BOLD3.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');


fc_matrix = corr(bold_ts');

figure('Position', [100 100 200 200]); % 设置图形大小
imagesc(fc_matrix);
colormap(jet); 
caxis([-1 1]);
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC3.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

mask = tril(true(size(fc_matrix)), -1);
fc_values = fc_matrix(mask); 
threshold = prctile(fc_values, 90); 
binary_fc = zeros(size(fc_matrix));
binary_fc(fc_matrix >= threshold & ~eye(size(fc_matrix))) = 1; % 排除对角线

figure('Position', [100 100 200 200]);
imagesc(binary_fc);
colormap(parula);
caxis([0 1]);        
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC3_threshold.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

%% 4
clc;clear;close;

data = cifti_read("D:\TravelClub\data\3RB2_A\sub-002\EyeOpen_PA_run4\EyeOpen_PA_run4_Atlas_s0.dtseries.nii");
bold_ts = data.cdata(1:20,1:100);

figure('Position', [100 100 200 200]);
hold on;
for i = 1:20
    plot(bold_ts(i,1:100), 'LineWidth', 0.5);
end

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k'); % 黑色边框
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\BOLD4.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');


fc_matrix = corr(bold_ts');

figure('Position', [100 100 200 200]); % 设置图形大小
imagesc(fc_matrix);
colormap(jet); 
caxis([-1 1]);
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC4.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

mask = tril(true(size(fc_matrix)), -1);
fc_values = fc_matrix(mask); 
threshold = prctile(fc_values, 90); 
binary_fc = zeros(size(fc_matrix));
binary_fc(fc_matrix >= threshold & ~eye(size(fc_matrix))) = 1; % 排除对角线

figure('Position', [100 100 200 200]);
imagesc(binary_fc);
colormap(parula);
caxis([0 1]);        
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC4_threshold.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');


%% 5: 计算并绘制四个FC矩阵的平均热图
close all
clc; clearvars -except data; close all;  % 保留之前的data变量

fc_matrices = cell(1,4);
datasets = {
    "D:\TravelClub\data\3RB2_A\sub-002\EyeOpen_PA_run1\EyeOpen_PA_run1_Atlas_s0.dtseries.nii";
    "D:\TravelClub\data\3RB2_A\sub-002\EyeClose_PA_run2\EyeClose_PA_run2_Atlas_s0.dtseries.nii";
    "D:\TravelClub\data\3RB2_A\sub-002\EyeClose_PA_run3\EyeClose_PA_run3_Atlas_s0.dtseries.nii";
    "D:\TravelClub\data\3RB2_A\sub-002\EyeOpen_PA_run4\EyeOpen_PA_run4_Atlas_s0.dtseries.nii"
};

for i = 1:4
    data = cifti_read(datasets{i});
    bold_ts = data.cdata(1:20, 1:100);
    fc_matrix = corr(bold_ts');
    
    mask = tril(true(size(fc_matrix)), -1);
    fc_values = fc_matrix(mask);
    threshold = prctile(fc_values, 90);
    binary_fc = zeros(size(fc_matrix));
    binary_fc(fc_matrix >= threshold & ~eye(size(fc_matrix))) = 1;
    
    binary_fc_all{i} = binary_fc;
end

mean_binary_fc = mean(cat(3, binary_fc_all{:}), 3);

figure('Position', [100 100 200 200]);
imagesc(mean_binary_fc);
colormap(parula);
caxis([0 1]);
axis square;

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); 
set(gca, 'Box', 'on', 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k');
title('');
colorbar('off');
box on;

exportgraphics(gcf, 'E:\Documents\Work\文章投稿\3RB\Figure\elements\FC_Mean.png', ...
    'Resolution', 600, 'BackgroundColor', 'none');

close all;