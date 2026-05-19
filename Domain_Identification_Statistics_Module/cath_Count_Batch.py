import pandas as pd
import sys
import os
import re

def process_cath_data(file_path):
    """
    处理Excel表格中的CATH标签数据并统计出现次数
    
    参数:
    file_path: Excel文件路径
    
    返回:
    dict: CATH标签统计字典
    """
    try:
        # 读取Excel文件
        df = pd.read_excel(file_path)
        
        # 检查数据列是否存在
        if df.shape[1] < 3:
            print(f"错误: Excel文件至少需要3列，但只有{df.shape[1]}列")
            return None
        
        # 获取第三列数据（索引为2），跳过第一行标题
        cath_data = df.iloc[1:, 2]  # 第一行是标题，从第二行开始
        
        # 初始化统计字典
        cath_count = {}
        
        # 处理每一行数据
        for item in cath_data:
            # 跳过空白行
            if pd.isna(item) or str(item).strip() == '':
                continue
                
            item_str = str(item).strip()
            
            # 按分号分割不同的CATH标签
            cath_entries = [entry.strip() for entry in item_str.split(';')]
            
            for entry in cath_entries:
                if entry:  # 确保不是空字符串
                    # 处理包含逗号的情况，只取逗号前的部分
                    if ',' in entry:
                        # 分割逗号，只取第一个部分
                        first_part = entry.split(',')[0].strip()
                        # 验证格式是否符合a.b.c或a.b.c.d
                        if is_valid_cath_format(first_part):
                            cath_count[first_part] = cath_count.get(first_part, 0) + 1
                    else:
                        # 直接验证格式
                        if is_valid_cath_format(entry):
                            cath_count[entry] = cath_count.get(entry, 0) + 1
        
        return cath_count
    
    except Exception as e:
        print(f"处理文件时出错: {e}")
        return None

def is_valid_cath_format(cath_string):
    """
    验证字符串是否符合CATH标签格式(a.b.c或a.b.c.d)
    
    参数:
    cath_string: 要验证的字符串
    
    返回:
    bool: 是否符合格式
    """
    # 使用正则表达式验证格式
    
    pattern = r'^\d+\.\d+\.\d+(\.\d+)?$'
    return bool(re.match(pattern, cath_string))

def save_results_to_excel(cath_dict, output_file):
    """
    将统计结果保存到Excel文件
    
    参数:
    cath_dict: CATH标签统计字典
    output_file: 输出文件路径
    """
    # 创建DataFrame
    result_df = pd.DataFrame(list(cath_dict.items()), 
                            columns=['CATH标签', '出现次数'])
    
    # 按出现次数降序排列
    result_df = result_df.sort_values('出现次数', ascending=False)
    
    # 保存到Excel
    result_df.to_excel(output_file, index=False)
    
    print(f"结果已保存到: {output_file}")
    print("\n统计结果:")
    print(result_df)

def print_usage():
    """打印使用说明"""
    print("用法: python cath_Count_Batch.py <输入Excel文件路径> [输出Excel文件路径]")
    print("示例:")
    print("  python cath_Count_Batch.py data.xlsx")
    print("  python cath_Count_Batch.py data.xlsx results.xlsx")

def main():

    # 检查命令行参数
    if len(sys.argv) < 2:
        print("错误: 请提供输入文件路径")
        print_usage()
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # 检查输入文件是否存在
    if not os.path.isfile(input_file):
        print(f"错误: 文件 '{input_file}' 不存在")
        sys.exit(1)
    
    # 设置输出文件路径
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # 如果没有指定输出文件，使用输入文件名加上_statistics后缀
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_statistics.xlsx"
    
    print(f"正在处理文件: {input_file}")
    
    # 处理数据
    cath_statistics = process_cath_data(input_file)
    
    if cath_statistics is None:
        print("处理失败，程序退出")
        sys.exit(1)
    
    # 保存结果
    save_results_to_excel(cath_statistics, output_file)
    
    # 打印详细统计信息
    print(f"\n统计摘要:")
    print(f"总共统计了 {len(cath_statistics)} 个不同的CATH标签")
    print(f"总出现次数: {sum(cath_statistics.values())}")

if __name__ == "__main__":
    main()