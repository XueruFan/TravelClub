# 本代码用来统计绘图旅行者们每次MRI扫描的情况
# 范雪如 Xue-Ru Fan 2024/5/28 @ Beijing Normal University

############################## 加载必要的包和数据文件 ##############################################
rm(list=ls())
library(ggplot2)
library(dplyr)
library(lubridate)
library(patchwork)

filefolder <- "C:/TravelClub/Data"
setwd(filefolder)
TestDate <- read.csv("TestDate.csv")
BasicInfo <- read.csv("BasicInfo.csv")

############################## 整理原始数据信息 ####################################################

# 重新排序TestDate的扫描次序
TestDate_clean <- TestDate %>%
  filter(Site != "H") %>% # 删去H机器的信息
  filter(!(P == 22 & Site == "J" & Date == "20240108")) %>% # 删去P22的第一次J扫描信息
  arrange(P, T) %>%
  group_by(P) %>%
  mutate(MT = row_number()) %>%
  ungroup() %>%
  select(-T)

# 按出生日期排序并生成新的被试编号
BasicInfo_new <- BasicInfo %>%
  arrange(desc(BirthDate)) %>%
  mutate(NP = row_number())

# 合并两个数据框
MergedData <- TestDate_clean %>%
  left_join(BasicInfo_new, by = "P")

# 计算年龄并添加Age列
MergedData <- MergedData %>%
  mutate(BirthDate = ymd(BirthDate),
         Date = ymd(Date),
         Age = round(interval(BirthDate, Date) / years(1), 2))

# 添加Interval_days列，计算每次测试和第一次测试之间的年龄差距（天数）
MergedData <- MergedData %>%
  group_by(P) %>%
  mutate(Interval_days = time_length(interval(first(Date), Date), "days")) %>%
  ungroup()

# 添加Interval_gap列，计算每次测试和上一次测试之间的年龄差距（天数）
MergedData <- MergedData %>%
  group_by(P) %>%
  mutate(Interval_gap = time_length(interval(lag(Date), Date), "days")) %>%
  ungroup()

# 添加AgeBase列，表示每个被试在MT=1时的年龄
MergedData <- MergedData %>%
  group_by(NP) %>%
  mutate(AgeBase = first(Age)) %>%
  ungroup()

# 重新排序列，将NP放在第一列
MergedData <- MergedData %>%
  select(NP, everything())


############################## 画扫描时间间隔图 ####################################################

ggplot(MergedData, aes(x = Interval_days, y = AgeBase, color = Site, group = NP)) +
  coord_fixed(ratio = 12, xlim = c(0, 160), ylim = c(20, 35)) + 
  geom_line(color = "lightgray", size = 1.5, alpha = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#86b5a1", "#e47159", "#7976a2", "#b95a58", "#3d5c6f", "#f9ae78",
                                "#963f5e", "#4282c6", "#00a391")) +
  labs(x = "", y = "", color = "Scanner Index") +
  scale_x_continuous(breaks = seq(0, 160, 40)) + 
  scale_y_continuous(breaks = seq(20, 35, 3)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'lines'),  # 调整图例键的大小
        legend.text = element_text(size = 12),  # 调整图例文本的大小
        legend.title = element_text(size = 14), # 调整图例标题的大小
        axis.text = element_text(size = 14)) + # 调整坐标轴刻度文字的大小
  guides(color = guide_legend(override.aes = list(size = 4)))  # 调整图例中点的大小

filename <- 'C:/TravelClub/Figure/MRIscanInterval.png'
ggsave(filename, width = 8, height = 8, units = "in", dpi = 300)


############################## 画志愿者第一次扫描时的年龄概率密度图 ################################

# 过滤出 MT = 1 时的 AgeBase 数据
age_base_data <- MergedData %>%
  filter(MT == 1)  %>%
  select(AgeBase, Sex)

# 绘制概率密度图 Ver 1 位于主图左侧
ggplot(age_base_data, aes(x = AgeBase, fill = Sex)) +
  geom_density(alpha = 0.4) +
  scale_y_reverse() +  # 翻转y轴，使图面向左边
  coord_flip() +
  scale_x_continuous(breaks = seq(20, 35, 3), limits = c(20, 35)) +
  scale_fill_manual(values = c("Male" = "#0064b5", "Female" = "#ff6347")) +  # 自定义颜色
  theme_minimal() +
  labs(x = "", y = "") +
  theme(legend.position = c(0.3, 0.85),  # 调整图例位置
    legend.background = element_rect(fill = alpha('white', 0.5), color = NA),  # 设置图例背景透明度并去掉黑框
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, 'lines'),  # 调整图例键的大小
    legend.text = element_text(size = 12),  # 调整图例文本的大小
    legend.title = element_text(size = 14) # 调整图例标题的大小
  )

filename <- 'C:/TravelClub/Figure/MRIparticipantV1.png'
ggsave(filename, width = 2, height = 8, units = "in", dpi = 300)


# 绘制概率密度图 Ver 2 作为独立图
ggplot(age_base_data, aes(x = AgeBase, fill = Sex)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(breaks = seq(20, 35, 3), limits = c(20, 35)) +
  scale_fill_manual(values = c("Male" = "#0064b5", "Female" = "#ff6347")) +  # 自定义颜色
  theme_minimal() +
  labs(x = "", y = "") +
  theme(legend.position = c(0.88, 0.7),  # 调整图例位置
        legend.background = element_rect(fill = alpha(NA, 0.5), color = NA),  # 设置图例背景透明度并去掉黑框
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key.size = unit(1.5, 'lines'),  # 调整图例键的大小
        legend.text = element_text(size = 12),  # 调整图例文本的大小
        legend.title = element_text(size = 14) # 调整图例标题的大小
  )

filename <- 'C:/TravelClub/Figure/MRIparticipantV2.png'
ggsave(filename, width = 8, height = 2, units = "in", dpi = 300)

############################## 计算 AgeBase 的平均数和标准差 #######################################

age_base_summary <- age_base_data %>%
  summarize(
    mean_AgeBase = round(mean(AgeBase, na.rm = TRUE), 2),
    sd_AgeBase = round(sd(AgeBase, na.rm = TRUE), 2)
  )

############################## 画两次扫描间隔的概率密度图 ##########################################

interval_days <- MergedData %>%
  filter(!is.na(Interval_gap)) %>%
  select(Interval_gap)

# 绘制概率密度图
ggplot(interval_days, aes(x = Interval_gap)) +
  scale_x_continuous(breaks = seq(2, 86, 14), limits = c(2, 85)) +
  geom_density(fill = "lightgray") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14), # 调整坐标轴刻度文字的大小
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())

filename <- 'C:/TravelClub/Figure/MRIscanDay.png'
ggsave(filename, width = 8, height = 2, units = "in", dpi = 300)

############################## 保存处理好的表格 ####################################################
write.csv(MergedData, file.path(filefolder, "MRIdetailInfo.csv"), row.names = F)
