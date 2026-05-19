import sys
import pandas as pd
import numpy as np
from scipy import stats
import os

def calculate_confidence_interval(data, confidence):
    """
    计算数据的T分布置信区间
    
    参数:
    data: 数据集 (list或array)
    confidence: 置信水平，默认为0.95 (95%)
    
    返回:
    tuple: (下限, 上限)
    """
    if len(data) < 2:
        # 数据点太少，无法计算置信区间
        return (min(data) if data else 0, max(data) if data else 0)
    
    # 计算均值和标准差
    mean_val = np.mean(data)
    std_val = np.std(data, ddof=1)  # 使用样本标准差
    
    # 计算标准误差
    n = len(data)
    se = std_val / np.sqrt(n)
    
    # 计算T分布的临界值
    # 使用双侧检验，所以是 (1 - confidence)/2
    t_critical = stats.t.ppf((1 + confidence) / 2, n - 1)
    
    # 计算置信区间
    lower_bound = mean_val - t_critical * se
    upper_bound = mean_val + t_critical * se
    
    return (lower_bound, upper_bound)

def is_outside_confidence_interval(data, value, confidence):
    """
    检查一个值是否在数据的T分布置信区间之外
    
    参数:
    data: 原始数据集
    value: 要检查的值
    confidence: 置信水平
    
    返回:
    bool: 如果在置信区间外返回True，否则返回False
    """
    lower, upper = calculate_confidence_interval(data, confidence)
    return value < lower or value > upper

def main():
    # 检查命令行参数
    if len(sys.argv) != 2:
        print("使用方法: python script.py <excel文件路径>")
        sys.exit(1)
    
    excel_file = sys.argv[1]
    
    # 检查文件是否存在
    if not os.path.exists(excel_file):
        print(f"错误: 文件 '{excel_file}' 不存在")
        sys.exit(1)
    
    try:
        # 读取Excel文件的第一工作表
        df = pd.read_excel(excel_file, sheet_name=0)
        print("成功读取Excel文件")
    except Exception as e:
        print(f"读取Excel文件时出错: {e}")
        sys.exit(1)
    
    # 获取列标题 (从第二列开始)
    headers = df.columns[1:].tolist()
    
    # 打印列标题
    print("\n列标题列表:")
    for i, header in enumerate(headers, 1):
        print(f"{i}. {header}")
    
    # 获取用户输入的分界线
    try:
        split_point = int(input(f"\n请输入分界线数字 (2-{len(headers)}): "))
        if split_point < 2 or split_point > len(headers):
            print(f"错误: 请输入2到{len(headers)}之间的数字")
            sys.exit(1)
    except ValueError:
        print("错误: 请输入有效的数字")
        sys.exit(1)
    
    # 获取用户输入的置信度
    try:
        confidence_level = float(input("\n请输入置信度 (例如: 0.95 表示95%, 0.99 表示99%): "))
        if confidence_level <= 0 or confidence_level >= 1:
            print("错误: 置信度必须在0和1之间 (例如0.95)")
            sys.exit(1)
    except ValueError:
        print("错误: 请输入有效的置信度数值")
        sys.exit(1)
    
    # 确定前组和后组的列索引
    # range是左闭右开的，所以前组是1到split_point-1，后组是split_point到最后
    front_group_indices = list(range(1, split_point))  # 列索引从1到split_point-1
    back_group_indices = list(range(split_point, len(headers) + 1))  # 列索引从split_point到最后
    
    # 检查前后组是否都有数据
    if not front_group_indices:
        print("错误: 前组没有列，请重新选择分界线")
        sys.exit(1)
        
    if not back_group_indices:
        print("错误: 后组没有列，请重新选择分界线")
        sys.exit(1)
    
    print(f"\n前组列数: {len(front_group_indices)} (列 {min(front_group_indices)} 到 {max(front_group_indices)})")
    print(f"后组列数: {len(back_group_indices)} (列 {min(back_group_indices)} 到 {max(back_group_indices)})")
    print(f"置信度: {confidence_level*100}%")
    
    # 存储结果的字典 - 现在会包含所有行的结果
    result_dict = {}
    
    # 统计各类结果的数量
    count_larger = 0
    count_smaller = 0
    count_no_difference = 0
    count_skipped = 0
    
    # 处理每一行数据
    for index, row in df.iterrows():
        cath_id = str(row[0])  # 第一列是CATH标签
        
        # 获取前组数据
        front_data = [row[i] for i in front_group_indices if pd.notna(row[i])]
        
        # 如果前组数据不足2个，无法计算置信区间，标记为跳过
        if len(front_data) < 2:
            result_dict[cath_id] = "数据不足"
            count_skipped += 1
            continue
        
        # 检查后组中是否有数据
        back_data_exists = False
        back_data = []
        for i in back_group_indices:
            if i < len(row) and pd.notna(row[i]):
                back_data_exists = True
                back_data.append(row[i])
        
        # 如果后组没有数据，标记为无后组数据
        if not back_data_exists:
            result_dict[cath_id] = "无后组数据"
            count_skipped += 1
            continue
        
        # 检查后组中的每个数字是否都在前组置信区间外
        all_outliers = True
        for value in back_data:
            # 检查该值是否在前组的置信区间内，使用用户指定的置信度
            if not is_outside_confidence_interval(front_data, value, confidence_level):
                all_outliers = False
                break
        
        # 根据条件设置结果
        if all_outliers:
            front_mean = np.mean(front_data)
            back_mean = np.mean(back_data)
            
            if back_mean > front_mean:
                result_dict[cath_id] = "偏大"
                count_larger += 1
            else:
                result_dict[cath_id] = "偏小"
                count_smaller += 1
        else:
            result_dict[cath_id] = "无显著差异"
            count_no_difference += 1
    
    # 输出结果统计
    print(f"\n分析完成 (置信度: {confidence_level*100}%):")
    print(f"- 偏大: {count_larger} 个")
    print(f"- 偏小: {count_smaller} 个")
    print(f"- 无显著差异: {count_no_difference} 个")
    print(f"- 跳过/无效: {count_skipped} 个")
    print(f"- 总计: {len(result_dict)} 个CATH标签")
    
    # 将结果保存到Excel文件
    if result_dict:
        output_df = pd.DataFrame(list(result_dict.items()), columns=['CATH标签', '比较结果'])
        output_file = f"t_distribution_analysis_result_{int(confidence_level*100)}percent.xlsx"
        output_df.to_excel(output_file, index=False)
        print(f"\n结果已保存到: {output_file}")
        

if __name__ == "__main__":
    main()