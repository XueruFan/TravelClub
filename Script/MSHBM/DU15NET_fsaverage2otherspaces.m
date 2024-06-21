%% This is the script which projects consensus map from fsaverage6 to different spaces - written by Jingnan Du 20240117

addpath(genpath('/ncf/sba10/MS-HBM/CBIG_CODE'))
addpath(genpath('/ncf/sba10/MS-HBM/code'))
addpath(genpath('/ncf/tools/0.10.0/code/templates/surface/fsaverage6'))
addpath(genpath('/n/sw/ncf/apps/freesurfer/6.0.0'));
addpath('/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/')
addpath('/ncf/cnl03/25/users/430/HCP_WB_Tutorial_1.0/gifti-1.6/')

wkdir = '/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/DU15NET_Final_Atlas/';
%atlases = {'DU15NET_consensus';'DU15NET_AgreeProb';'DU15NET_Agree80per';'DU15NET_Agree53per';'DU15NET_consensus_raw';'DU15NET_Prior'};
atlases = {'DU15NET_Tripartite'};

for j = 1:length(atlases)
    

    atlas_name = atlases{j};
    map_file = [wkdir atlas_name '.dscalar.nii']; % the map file should be a atlas in fsaverage6 space
 
    g = ciftiopen(map_file,'/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command',1);
    img = g.cdata;
    lh_FS_data = img(1:40962,:)';
    rh_FS_data = img(40963:end,:)';
    folder_to_write = [wkdir 'tmp'];
    if exist(folder_to_write,'dir')
        rmdir(folder_to_write,'s')
    end
    
    if ~strcmp(atlas_name,'DU15NET_AgreeProb')
        [lh_label_fsLR_32k,rh_label_fsLR_32k,lh_label_fsLR_164k,rh_label_fsLR_164k]=CBIG_project_fsaverage2fsLR(lh_FS_data',rh_FS_data','fsaverage6','label',folder_to_write);
    elseif strcmp(atlas_name,'DU15NET_AgreeProb')
        [lh_label_fsLR_32k,rh_label_fsLR_32k,lh_label_fsLR_164k,rh_label_fsLR_164k]=CBIG_project_fsaverage2fsLR(lh_FS_data',rh_FS_data','fsaverage6','metric',folder_to_write);
        lh_label_fsLR_32k = lh_label_fsLR_32k./15;
        rh_label_fsLR_32k = rh_label_fsLR_32k./15;
        lh_label_fsLR_164k = lh_label_fsLR_164k./15;
        rh_label_fsLR_164k = rh_label_fsLR_164k./15;
    end
    
    % save fsLR_32k
    g = ciftiopen('/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/HCP_WB_Tutorial_1.5_Pr_kN3mg/Q1-Q6_R440.sulc.32k_fs_LR.dscalar.nii','/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command',1);
    g.cdata = [lh_label_fsLR_32k;rh_label_fsLR_32k];
    FILENAME = [wkdir 'HCP/fsLR_32k/' atlas_name '_fsLR_32k'];
    ciftisavereset(g,[FILENAME,'.dscalar.nii'],'/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command');
    if ~strcmp(atlas_name,'DU15NET_Tripartite')
        colorfile='/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/ColorMap_16_uncertain.txt ';
    else
        colorfile='/ncf/sba10/MS-HBM/ColorMap/Flechsig.txt ';  
    end
    system(['/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command -cifti-label-import ' fullfile([FILENAME '.dscalar.nii ']) colorfile fullfile([FILENAME '.dlabel.nii'])])
    dlabel = [FILENAME '.dlabel.nii ']; % create .gii
    system(['/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command -cifti-separate ' dlabel ' COLUMN -label CORTEX_LEFT ' dlabel(1:end-12) '_lh.label.gii -label CORTEX_RIGHT ' dlabel(1:end-12) '_rh.label.gii'])
    
    % save fsLR_164k
    g = ciftiopen('/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/HCP_WB_Tutorial_1.5_Pr_kN3mg/Q1-Q6_R440.sulc.164k_fs_LR.dscalar.nii','/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command',1);
    g.cdata = [lh_label_fsLR_164k;rh_label_fsLR_164k];
    FILENAME = [wkdir 'HCP/fsLR_164k/' atlas_name '_fsLR_164k'];
    ciftisavereset(g,[FILENAME,'.dscalar.nii'],'/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command');
      if ~strcmp(atlas_name,'DU15NET_Tripartite')
        colorfile='/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/ColorMap_16_uncertain.txt ';
    else
        colorfile='/ncf/sba10/MS-HBM/ColorMap/Flechsig.txt ';  
    end
    system(['/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command -cifti-label-import ' fullfile([FILENAME '.dscalar.nii ']) colorfile fullfile([FILENAME '.dlabel.nii'])])
    %OUTDIR = '/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/DU15NET_Final_Atlas';
    dlabel = [FILENAME '.dlabel.nii '];  % create .gii
    system(['/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command -cifti-separate ' dlabel ' COLUMN -label CORTEX_LEFT ' dlabel(1:end-12) '_lh.label.gii -label CORTEX_RIGHT ' dlabel(1:end-12) '_rh.label.gii'])
    
    %% fsaverage6 to fsaverage5 and also fsaverage space - Jingnan successfully project fsaverage6 to fsaverage and fsaverage5 !! -  Date: 240120
    wkdir = '/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/DU15NET_Final_Atlas/';
    g = ciftiopen(map_file,'/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command',1);
    img = g.cdata;
    lh_label = img(1:40962,:); rh_label = img(40963:end,:);
    %lh_labels6 = lh_label; rh_labels6 = rh_label;
    % ------------------------------------------------- upsample from fsaverage6 to fsaverage
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex'); rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'white', 'cortex');
    lh_mesh6 = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'white', 'cortex'); rh_mesh6 = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'white', 'cortex');
    lh_labels7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices,lh_mesh6,lh_label');
    rh_labels7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices,rh_mesh6,rh_label');
    % ------------------------------------------------- upsample from fsaverage6 to fsaverage5
    lh_mesh5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'white', 'cortex'); rh_mesh5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'white', 'cortex');
    lh_mesh6 = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'white', 'cortex'); rh_mesh6 = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'white', 'cortex');
    lh_labels5 = MARS_NNInterpolate_kdTree(lh_mesh5.vertices,lh_mesh6,lh_label');
    rh_labels5 = MARS_NNInterpolate_kdTree(rh_mesh5.vertices,rh_mesh6,rh_label');
    % ------------------------------------------------- write fsaverage annot using .mat file which contains labels in fsaverage5
    lh_labels = lh_labels5'; rh_labels  = rh_labels5';
    matfile = [wkdir atlas_name '_fsaverage5.mat'];
    save(matfile,'lh_labels','rh_labels');
    lh_output_fsaverage = [wkdir 'FreeSurfer/fsaverage/label/lh.' atlas_name '_fsaverage.annot'];
    rh_output_fsaverage = [wkdir 'FreeSurfer/fsaverage/label/rh.' atlas_name '_fsaverage.annot'];
    %c This function works!!But note that I change the color to be consistent
    % with DU15NET / Originally it used FreeSurfer color
    addpath('/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET') % cd here because I am using a modified function, not the original one
    
    if ~strcmp(atlas_name,'DU15NET_AgreeProb')
        CBIG_SaveParcellationToFreesurferAnnotation(matfile,lh_output_fsaverage,rh_output_fsaverage)
        %-------------------------------------------- write annot files for fsaverage5 and fsaverage6
        [~, lh, lh_s] = read_annotation(lh_output_fsaverage);
        [~, rh, rh_s] = read_annotation(rh_output_fsaverage);
        % fsaverage 5
        disp('Writing annot files for fsaverage5 and fsaverage6.');
        write_annotation([wkdir 'FreeSurfer/fsaverage5/label/lh.' atlas_name '_fsaverage5.annot'], 0:10241, lh(1:10242), lh_s);
        write_annotation([wkdir 'FreeSurfer/fsaverage5/label/rh.' atlas_name '_fsaverage5.annot'], 0:10241, rh(1:10242), rh_s);
        % fsaverage 6
        write_annotation([wkdir 'FreeSurfer/fsaverage6/label/lh.' atlas_name '_fsaverage6.annot'], 0:40961, lh(1:40962), lh_s);
        write_annotation([wkdir 'FreeSurfer/fsaverage6/label/rh.' atlas_name '_fsaverage6.annot'], 0:40961, rh(1:40962), rh_s);
        
    elseif ~strcmp(atlas_name,'DU15NET_AgreeProb')
    end
    %-------------------------------------------- done
    %% save text files
    g = ciftiopen(map_file,'/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command',1);
    img = g.cdata;
    lh_labels6 = img(1:40962,:); rh_labels6 = img(40963:end,:);
    fsdir = [wkdir 'FreeSurfer/'];
    type = 'fsaverage6';
    
    if ~strcmp(atlas_name,'DU15NET_AgreeProb')
        lh = lh_labels6 ; rh =  rh_labels6;
    elseif ~strcmp(atlas_name,'DU15NET_AgreeProb')
        lh = lh_labels6./15 ; rh =  rh_labels6./15;
    end
    
    dlmwrite([fsdir type '/lh.' atlas_name '_' type '.txt'], lh, 'delimiter','\t','newline','pc');
    dlmwrite([fsdir type '/rh.' atlas_name '_' type '.txt'], rh, 'delimiter','\t','newline','pc');
    
    type = 'fsaverage5';
    if ~strcmp(atlas_name,'DU15NET_AgreeProb')
        lh = lh_labels5' ; rh =  rh_labels5';
    elseif ~strcmp(atlas_name,'DU15NET_AgreeProb')
        lh = lh_labels5'./15 ; rh =  rh_labels5'./15;
    end
    dlmwrite([fsdir type '/lh.' atlas_name '_' type '.txt'], lh, 'delimiter','\t','newline','pc');
    dlmwrite([fsdir type '/rh.' atlas_name '_' type '.txt'], rh, 'delimiter','\t','newline','pc');
    
    type = 'fsaverage';
    if ~strcmp(atlas_name,'DU15NET_AgreeProb')
        lh = lh_labels7'; rh =  rh_labels7';
    elseif ~strcmp(atlas_name,'DU15NET_AgreeProb')
        lh = lh_labels7'./15; rh =  rh_labels7'./15;
    end
    dlmwrite([fsdir type '/lh.' atlas_name '_' type '.txt'], lh, 'delimiter','\t','newline','pc');
    dlmwrite([fsdir type '/rh.' atlas_name '_' type '.txt'], rh, 'delimiter','\t','newline','pc');
    %% surface (fsaverage6) to volume (MNI) - Jingnan successfully project fsaverage6 to MNI -  Date: 240213
    % I use the registration methods from Wu2017
    % /ncf/sba10/MS-HBM/CBIG_CODE/stable_projects/registration/Wu2017_RegistrationFusion/examples
    % T1: /ncf/sba10/MS-HBM/CBIG_CODE/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz
    
    g = ciftiopen(map_file,'/n/sw/fasrcsw/apps/Core/workbench/1.0-fasrc01/bin_rh_linux64/wb_command',1);
    img = g.cdata; lh_label = img(1:40962,:); rh_label = img(40963:end,:);
    % ------------------------------------------------- upsample from fsaverage6 to fsaverage
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');
    rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'white', 'cortex');
    lh_mesh6 = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'white', 'cortex');
    rh_mesh6 = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'white', 'cortex');
    lh_labels7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices,lh_mesh6,lh_label');
    rh_labels7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices,rh_mesh6,rh_label');
    
    if ~strcmp(atlas_name,'DU15NET_AgreeProb')
        addpath('/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET')
        Du_project_to_MNI(lh_label', rh_label', [wkdir 'MNI/' atlas_name])
        
    elseif strcmp(atlas_name,'DU15NET_AgreeProb')
        addpath('/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET')
        Du_project_to_MNI(lh_label', rh_label', [wkdir 'MNI/' atlas_name])
        for i = 1:3
            tmp = MRIread(['/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/DU15NET_Final_Atlas/MNI/DU15NET_AgreeProb_MNI152_' num2str(i) 'mm.nii.gz']);
            tmp_vol = tmp.vol;
            tmp.vol = tmp_vol./15;
            MRIwrite(tmp,['/ncf/sba09/MS-HBM/iCere_Full/Analysis/Cortical/DU15NET/DU15NET_Final_Atlas/MNI/DU15NET_AgreeProb_MNI152_' num2str(i) 'mm.nii.gz']);
        end
    end
   
end