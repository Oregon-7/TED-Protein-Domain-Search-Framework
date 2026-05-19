import requests
import pandas as pd
import sys
import time
from typing import List, Dict, Any
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def search_ted_database(uniprot_id: str) -> Dict[str, Any]:
    """从TED数据库数据中提取CATH标签
    
    Args:
        uniprot_id (str): UniProt ID

    Returns:
        Dict[str, Any]: TED数据库返回的JSON数据
    """ 
    try:
        url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_id}"
        response = requests.get(url, verify=False)
        
        if response.status_code == 200:
            return response.json()
        else:
            print(f"error:TED请求失败，状态码：{response.status_code}")
            return None
            
    except Exception as e:
        print(f"error:TED搜索错误: {e}")
        return None

def extract_cath_labels(ted_data: Dict[str, Any]) -> List[str]:
    cath_labels = []
    
    if 'data' in ted_data and isinstance(ted_data['data'], list):
        for item in ted_data['data']:
            if 'cath_label' in item and item['cath_label'] != "-":
                cath_labels.append(item['cath_label'])
    
    return cath_labels

def process_uniprot_list(excel_file_path: str) -> Dict[str, List[str]]:
    """
    处理Excel文件中的UniProt编号列表，返回UniProt编号到CATH标签的映射字典
    
    Args:
        excel_file_path: Excel文件路径
        
    Returns:
        字典，键为UniProt编号，值为CATH标签列表
    """
    # 读取Excel文件
    try:
        df = pd.read_excel(excel_file_path)
    except Exception as e:
        print(f"读取Excel文件时出错: {e}")
        return {}
    
    uniprot_column = df.iloc[:, 0]
    uniprot_ids = uniprot_column.dropna().astype(str).tolist()
    
    print(f"从Excel文件中读取了 {len(uniprot_ids)} 个UniProt编号")
    print(f"UniProt编号示例: {uniprot_ids[:5]}")  # 显示前5个UniProt编号
    
    # 初始化结果字典
    result = {}
    
    # 处理每个UniProt编号
    for i, uniprot_id in enumerate(uniprot_ids, 1):
        print(f"\n处理进度: {i}/{len(uniprot_ids)} - 当前UniProt编号: {uniprot_id}")
        
        # 直接获取TED数据库数据
        ted_data = search_ted_database(uniprot_id)
        if not ted_data:
            print(f"  警告: 无法获取UniProt编号 {uniprot_id} 的TED数据")
            result[uniprot_id] = []
            continue
        
        # 提取CATH标签
        cath_labels = extract_cath_labels(ted_data)
        result[uniprot_id] = cath_labels
        
        print(f"  提取到 {len(cath_labels)} 个CATH标签: {cath_labels}")

    
    return result

def save_results_to_excel(results: Dict[str, List[str]], output_file: str = "results.xlsx"):
    """
    将结果字典保存到Excel文件
    
    Args:
        results: UniProt编号到CATH标签列表的映射字典
        output_file: 输出Excel文件名
    """
    # 准备数据用于创建DataFrame
    data = []
    for uniprot_id, cath_labels in results.items():
        # 将CATH标签列表转换为字符串，用分号分隔
        # 如果没有CATH标签，则留空
        cath_labels_str = "; ".join(cath_labels) if cath_labels else ""
        data.append({
            "UniProt编号": uniprot_id,
            "CATH标签": cath_labels_str,
            "CATH标签数量": len(cath_labels)
        })
    
    # 创建DataFrame
    df = pd.DataFrame(data)
    
    # 保存到Excel文件
    try:
        df.to_excel(output_file, index=False)
        print(f"\n结果已保存到 {output_file}")
        print(f"总计处理了 {len(results)} 个UniProt编号")
    except Exception as e:
        print(f"保存Excel文件时出错: {e}")

def main():
    """
    主函数：处理命令行输入
    """
    if len(sys.argv) < 2:
        print("使用方法: python cath_Batch_byUniprotID.py <Excel文件路径>")
        print("示例: python cath_Batch_byUniprotID.py uniprot_ids.xlsx")
        sys.exit(1)
    
    excel_file_path = sys.argv[1]
    
    print("=" * 60)
    print("UniProt编号到CATH标签批量处理工具")
    print("=" * 60)
    
    # 处理UniProt编号列表
    result = process_uniprot_list(excel_file_path)
    
    # 保存结果到Excel
    save_results_to_excel(result)
    
    # 显示结果摘要
    print("\n处理摘要:")
    total_uniprots = len(result)
    uniprots_with_cath = sum(1 for cath_labels in result.values() if cath_labels)
    print(f"总UniProt编号数: {total_uniprots}")
    print(f"找到CATH标签的UniProt编号数: {uniprots_with_cath}")
    print(f"未找到CATH标签的UniProt编号数: {total_uniprots - uniprots_with_cath}")

if __name__ == "__main__":
    main()