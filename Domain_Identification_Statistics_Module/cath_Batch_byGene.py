import requests
import pandas as pd
import sys
import time
from typing import List, Dict, Any

def get_uniprot_id(gene_id: str) -> str:
    """从基因编号映射对应的Uniprot编号

    Args:
        gene_id (str): 基因编号

    Returns:
        str: Uniprot编号，失败则返回None
    """
    try:
        query = f'gene:{gene_id}'
        response = requests.get(
            'https://rest.uniprot.org/uniprotkb/search',
            params={
                'query': query,
                'fields': 'accession',
                'format': 'json',
                'size': 1
            },
            timeout=30
        )
        response.raise_for_status()
        
        data = response.json()
        if data.get('results') and len(data['results']) > 0:
            return data['results'][0]['primaryAccession']
        else:
            return None
            
    except Exception as e:
        print(f"error:{e}")
        return None


def search_ted_database(uniprot_id: str) -> Dict[str, Any]:
    """从Uniprot编号映射对应的TED数据库数据

    Args:
        uniprot_id (str): Uniprot编号

    Returns:
        Dict[str, Any]: TED数据库返回的json文件，失败则返回None
    """
    try:
        url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            return response.json()
        else:
            print(f"error:TED请求失败，状态码：{response.status_code}")
            return None
            
    except Exception as e:
        print(f"error:TED搜索错误: {e}")
        return None

def extract_cath_labels(ted_data: Dict[str, Any]) -> List[str]:
    """从TED数据库数据中提取CATH标签
    Args:
        ted_data (Dict[str, Any]): TED数据库返回的json文件

    Returns:
        List[str]: CATH标签列表，忽略默认空值与"-"空值
    """
    cath_labels = []
    
    if 'data' in ted_data and isinstance(ted_data['data'], list):
        for item in ted_data['data']:
            if 'cath_label' in item and item['cath_label'] != "-":
                cath_labels.append(item['cath_label'])
    
    return cath_labels

def process_gene_list(excel_file_path: str) -> Dict[str, List[str]]:
    """处理Excel文件中的基因列表，返回基因编号到CATH标签的映射字典
    
    Args:
        excel_file_path (str): Excel文件路径
        
    Returns:
         (dict)，键为基因编号，值为CATH标签列表
    """
    # 读取Excel文件
    try:
        df = pd.read_excel(excel_file_path)
    except Exception as e:
        print(f"读取Excel文件时出错: {e}")
        return {}
    
    # 假设基因编号在第一列
    gene_column = df.iloc[:, 0]
    gene_ids = gene_column.dropna().astype(str).tolist()
    
    print(f"从Excel文件中读取了 {len(gene_ids)} 个基因编号")
    print(f"基因编号示例: {gene_ids[:5]}")  # 显示前5个基因编号
    
    # 初始化结果字典
    result = {}
    
    # 处理每个基因编号
    for i, gene_id in enumerate(gene_ids, 1):
        print(f"\n处理进度: {i}/{len(gene_ids)} - 当前基因: {gene_id}")
        
        # 获取UniProt编号
        uniprot_id = get_uniprot_id(gene_id)
        if not uniprot_id:
            print(f"  警告: 无法获取基因 {gene_id} 的UniProt编号")
            result[gene_id] = []
            continue
        
        print(f"  获取到UniProt编号: {uniprot_id}")
        
        # 获取TED数据库数据
        ted_data = search_ted_database(uniprot_id)
        if not ted_data:
            print(f"  警告: 无法获取UniProt编号 {uniprot_id} 的TED数据")
            result[gene_id] = []
            continue
        
        # 提取CATH标签
        cath_labels = extract_cath_labels(ted_data)
        result[gene_id] = cath_labels
        
        print(f"  提取到 {len(cath_labels)} 个CATH标签: {cath_labels}")
        
        # 添加延迟以避免请求过于频繁
        time.sleep(0.2)
    
    return result

def save_results_to_excel(results: Dict[str, List[str]], output_file: str = "results.xlsx"):
    """
    将结果字典保存到Excel文件
    
    Args:
        results: 基因编号到CATH标签列表的映射字典
        output_file: 输出Excel文件名
    """
    # 准备数据用于创建DataFrame
    data = []
    for gene_id, cath_labels in results.items():
        # 将CATH标签列表转换为字符串，用分号分隔
        # 如果没有CATH标签，则留空
        cath_labels_str = "; ".join(cath_labels) if cath_labels else ""
        data.append({
            "基因编号": gene_id,
            "CATH标签": cath_labels_str,
            "CATH标签数量": len(cath_labels)
        })
    
    # 创建DataFrame
    df = pd.DataFrame(data)
    
    # 保存到Excel文件
    try:
        df.to_excel(output_file, index=False)
        print(f"\n结果已保存到 {output_file}")
        print(f"总计处理了 {len(results)} 个基因编号")
    except Exception as e:
        print(f"保存Excel文件时出错: {e}")

def main():
    """
    主函数：处理命令行输入
    """
    if len(sys.argv) < 2:
        print("使用方法: python cath_Batch_byGene.py <Excel文件路径>")
        print("示例: python cath_Batch_byGene.py gene_List.xlsx")
        sys.exit(1)
    
    excel_file_path = sys.argv[1]
    
    print("=" * 60)
    print("基因编号到CATH标签批量处理工具")
    print("=" * 60)
    
    # 处理基因列表
    result = process_gene_list(excel_file_path)
    
    # 保存结果到Excel
    save_results_to_excel(result)
    
    # 显示结果摘要
    print("\n处理摘要:")
    total_genes = len(result)
    genes_with_cath = sum(1 for cath_labels in result.values() if cath_labels)
    print(f"总基因数: {total_genes}")
    print(f"找到CATH标签的基因数: {genes_with_cath}")
    print(f"未找到CATH标签的基因数: {total_genes - genes_with_cath}")

if __name__ == "__main__":
    main()