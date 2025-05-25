#!/bin/bash
# Copy the necessary files from 44T
# Written by Xue-Ru Fan 2025-02-27 @ Beijing Normal University

# 定义源路径和目标路径的基础部分
SOURCE_BASE="D:/3RB2_X"
DEST_BASE="D:/TravelClub/data"

# 遍历所有3RB2_X文件夹
for X in A; do # for X in A B C D E G; do
    # 定义当前3RB2_X的路径
    CURRENT_3RB2="$SOURCE_BASE/$X"
    
    # 检查当前3RB2_X文件夹是否存在
    if [ -d "$CURRENT_3RB2" ]; then
        # 遍历当前3RB2_X下的所有sub-XXX文件夹
        for SUB in "$CURRENT_3RB2"/sub-*; do
            # 检查是否是文件夹
            if [ -d "$SUB" ]; then
                SOURCE_PATH="$SUB/ciftify/$(basename "$SUB")/MNINonLinear/Results"
                
                # 检查源路径是否存在
                if [ -d "$SOURCE_PATH" ]; then
                    # 遍历源路径下的所有文件夹
                    for FOLDER in "$SOURCE_PATH"/*; do
                        # 检查是否是文件夹并且不包含“echo”
                        if [ -d "$FOLDER" ] && [[ "$(basename "$FOLDER")" != *"echo"* ]]; then
                            # 创建目标路径
                            DEST_PATH="$DEST_BASE/$X/$(basename "$SUB")/$(basename "$FOLDER")"
                            mkdir -p "$DEST_PATH"
                            # 复制文件夹
                            cp -r "$FOLDER" "$DEST_PATH"
                            echo "Copied $FOLDER to $DEST_PATH"
                        fi
                    done
                else
                    echo "Source path $SOURCE_PATH does not exist."
                fi
            fi
        done
    else
        echo "3RB2_X folder $CURRENT_3RB2 does not exist."
    fi
done