# TravelClub PFM Analysis Pipeline

Written by: Xueru Fan

依次运行下列代码

1. [`copy_3RB2_data.sh`](code/copy_3RB2_data.sh)
   
首先，把要用的数据从预处理好的文件夹中拷贝出来，方便后续分析时读入

2. [`TravelClub_01_Prepare.m`](code/TravelClub_01_Prepare.m)
   
Step 0: Set up environment and parameters (设置环境变量，提供站点、被试编号，同时需要设置网络分区个数和卡阈值的大小)

Step 1: prepare the fMRI list files(生成站点被试的fMRI数据地址list文件)

Step 2: Generating profiles and initialization parameters (根据不同的阈值生成功能连接矩阵，由于多回波ICA融合时会进行降噪，在step0里需要小心设置卡阈值的范围)

Step 3: Prepare the individual FC list files (把每个被试的功能连接结果路径按照要求整理成txt文件)

3. [`TravelClub_02_PFM_loop.m`](code/TravelClub_02_PFM_loop.m) or [`TravelClub_02_PFM_uni.m`](c0de/TravelClub_02_PFM_uni.m)
   
如果使用循环依次处理 不同阈值、不同网络个数 使用loop代码；如果开多个窗口同时处理 不同阈值、不同网络个数 使用uni代码

Step 5: To run Yeo2011 clustering algorithm for generate our own group prior

Step 6: Estimate individual mapping (估计模型的参数同时得出个体分区)

4. [`TravelClub_03_GetLabel.m`](code/TravelClub_03_GetLabel.m)
   
Step 7: Label Extract individual network labels (提取出个体分区的结果)

Step 8: Make individual visualizations (可视化)

5. .m

自定义不同网络的名称和颜色，还没写这个代码
