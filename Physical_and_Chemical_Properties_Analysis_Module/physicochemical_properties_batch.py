import requests
import pandas as pd
import time
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pathlib import Path
from requests.exceptions import RequestException


REQUEST_TIMEOUT = 30
REQUEST_RETRIES = 3
RETRY_BACKOFF_SECONDS = 1.0


def parse_chopping(chopping: str) -> list:
    """
    解析chopping字符串，返回残基区间列表。
    
    格式示例：
        "9-85" → [(9, 85)]
        "9-85_100-133" → [(9, 85), (100, 133)]
    
    Args:
        chopping: CATH结构域位置字符串
    
    Returns:
        列表，每个元素为(start, end)元组，均为基于1的索引
    """
    if not isinstance(chopping, str) or not chopping.strip():
        raise ValueError(f"Invalid chopping value: {chopping}")

    segments = chopping.split('_')
    ranges = []
    for seg in segments:
        seg = seg.strip()
        parts = seg.split('-')
        if len(parts) != 2:
            raise ValueError(f"Invalid chopping segment format: {seg}")
        start, end = map(int, parts)
        if start <= 0 or end <= 0 or start > end:
            raise ValueError(f"Invalid chopping segment range: {seg}")
        ranges.append((start, end))

    if not ranges:
        raise ValueError(f"No valid chopping ranges found: {chopping}")
    return ranges


def get_with_retries(url: str, timeout: int = REQUEST_TIMEOUT, retries: int = REQUEST_RETRIES):
    """对GET请求进行重试，处理瞬时网络错误和服务端暂时性错误。"""
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=timeout)
            # 429和5xx通常可重试
            if response.status_code == 429 or 500 <= response.status_code < 600:
                last_error = Exception(f"HTTP {response.status_code}")
            else:
                return response
        except RequestException as e:
            last_error = e

        if attempt < retries:
            time.sleep(RETRY_BACKOFF_SECONDS * attempt)

    raise Exception(f"Request failed after {retries} attempts: {url}; reason={last_error}")


def extract_cath_domains(protein_data: dict, target_cath: str) -> list:
    """
    从CATH API返回的数据中提取与目标CATH编号匹配的结构域信息。
    
    Args:
        protein_data: CATH API返回的JSON数据
        target_cath: 用户输入的CATH编号（如"1.10.10"或"1.10.150.10"）
    
    Returns:
        列表，每个元素为(dict)包含'ted_id'和'chopping'等信息
    """
    matching_domains = []
    if 'data' not in protein_data:
        return matching_domains
    
    for domain in protein_data['data']:
        if domain.get('cath_label') == target_cath:
            matching_domains.append({
                'ted_id': domain.get('ted_id'),
                'chopping': domain.get('chopping'),
                'nres_domain': domain.get('nres_domain'),
                'plddt': domain.get('plddt')
            })
    return matching_domains


def fetch_protein_sequence(uniprot_id: str) -> str:
    """
    从UniProt API获取蛋白质的完整一级序列。
    
    Args:
        uniprot_id: UniProt ID
    
    Returns:
        单字母表示的蛋白质序列字符串
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = get_with_retries(url)
    
    if response.status_code != 200:
        raise Exception(f"Failed to fetch sequence for {uniprot_id}: HTTP {response.status_code}")
    
    # 解析FASTA格式：第一行是注释，后续行是序列
    lines = response.text.strip().split('\n')
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    if not sequence:
        raise Exception(f"Empty sequence returned for {uniprot_id}")
    
    return sequence


def extract_domain_sequence(full_sequence: str, chopping: str) -> str:
    """
    根据chopping位置从完整序列中截取结构域序列。
    
    支持单段（如"9-85"）和双段（如"9-85_100-133"）格式。
    蛋白序列位置从1开始编号。
    
    Args:
        full_sequence: 完整的蛋白质序列
        chopping: 结构域位置字符串
    
    Returns:
        截取后的结构域序列
    """
    ranges = parse_chopping(chopping)
    if not ranges:
        raise ValueError(f"No ranges parsed from chopping: {chopping}")
    domain_parts = []
    
    for start, end in ranges:
        # 转换为Python的0-based索引
        start_idx = start - 1
        end_idx = end  # Python切片是左闭右开，所以直接用end作为结束索引
        
        if start_idx < 0 or end_idx > len(full_sequence):
            raise Exception(f"Chopping range {start}-{end} out of sequence length {len(full_sequence)}")
        
        domain_parts.append(full_sequence[start_idx:end_idx])
    
    domain_sequence = ''.join(domain_parts)
    if not domain_sequence:
        raise ValueError(f"Empty domain sequence extracted from chopping: {chopping}")
    return domain_sequence


def calculate_protein_properties(sequence: str) -> dict:
    """
    使用Biopython ProtParam计算蛋白质的理化性质。
    
    Args:
        sequence: 单字母表示的蛋白质序列
    
    Returns:
        包含各项理化性质的字典
    """
    analysis = ProteinAnalysis(sequence)
    
    # 基础性质
    properties = {
        'length': len(sequence),
        'molecular_weight': round(analysis.molecular_weight(), 2),
        'isoelectric_point': round(analysis.isoelectric_point(), 2),
    }
    
    # 氨基酸计数和百分比
    aa_count = analysis.count_amino_acids()
    if hasattr(analysis, 'get_amino_acids_percent'):
        # 兼容旧版Biopython
        aa_percent_raw = analysis.get_amino_acids_percent()
    elif hasattr(analysis, 'amino_acids_percent'):
        # 兼容新版Biopython
        aa_percent_raw = analysis.amino_acids_percent
    else:
        raise AttributeError("ProteinAnalysis missing amino acid percent API")

    # 某些版本返回0-1比例，某些版本返回0-100百分比，这里统一到0-100。
    if aa_percent_raw and max(aa_percent_raw.values()) <= 1.0:
        aa_percent = {aa: round(pct * 100, 2) for aa, pct in aa_percent_raw.items()}
    else:
        aa_percent = {aa: round(float(pct), 2) for aa, pct in aa_percent_raw.items()}
    properties['aa_count'] = aa_count
    properties['aa_percent'] = aa_percent
    
    # 带电残基统计
    charged_residues = analysis.charge_at_pH(7.0) if hasattr(analysis, 'charge_at_pH') else None
    
    # 芳香性指数
    properties['aromaticity'] = round(analysis.aromaticity(), 4)
    
    # 不稳定性指数
    properties['instability_index'] = round(analysis.instability_index(), 2)
    
    # GRAVY（亲水性指数）
    properties['gravy'] = round(analysis.gravy(), 4)
    
    # 消光系数（还原态和氧化态）
    ec_reduced, ec_oxidized = analysis.molar_extinction_coefficient()
    properties['extinction_coefficient_reduced'] = ec_reduced
    properties['extinction_coefficient_oxidized'] = ec_oxidized
    
    # 二级结构倾向分数（螺旋、转角、折叠）
    helix, turn, sheet = analysis.secondary_structure_fraction()
    properties['helix_fraction'] = round(helix, 4)
    properties['turn_fraction'] = round(turn, 4)
    properties['sheet_fraction'] = round(sheet, 4)
    
    # 带电残基计数
    neg_charged = aa_count.get('D', 0) + aa_count.get('E', 0)
    pos_charged = aa_count.get('K', 0) + aa_count.get('R', 0) + aa_count.get('H', 0)
    properties['negatively_charged_count'] = neg_charged
    properties['positively_charged_count'] = pos_charged
    
    return properties


def properties_to_flat_dict(seq_key: str, sequence: str, properties: dict) -> dict:
    """
    将嵌套的性质字典展平为单行记录。
    
    Args:
        seq_key: 序列的键名
        sequence: 蛋白质序列
        properties: calculate_protein_properties返回的嵌套字典
    
    Returns:
        展平后的字典，适合写入DataFrame
    """
    flat = {
        'key': seq_key,
        'sequence': sequence,
        'length': properties['length'],
        'molecular_weight': properties['molecular_weight'],
        'isoelectric_point': properties['isoelectric_point'],
        'aromaticity': properties['aromaticity'],
        'instability_index': properties['instability_index'],
        'gravy': properties['gravy'],
        'extinction_coefficient_reduced': properties['extinction_coefficient_reduced'],
        'extinction_coefficient_oxidized': properties['extinction_coefficient_oxidized'],
        'helix_fraction': properties['helix_fraction'],
        'turn_fraction': properties['turn_fraction'],
        'sheet_fraction': properties['sheet_fraction'],
        'negatively_charged_count': properties['negatively_charged_count'],
        'positively_charged_count': properties['positively_charged_count'],
    }
    
    # 添加20种标准氨基酸的百分比
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        flat[f'aa_percent_{aa}'] = properties['aa_percent'].get(aa, 0.0)
    
    return flat


def main():
    """主函数：执行完整的批处理流程"""
    
    # 0. 用户输入
    excel_path = input("请输入Excel文件路径: ").strip()
    target_cath = input("请输入三到四位的CATH编号（如 1.10.10 或 1.10.150.10）: ").strip()

    # 验证CATH编号格式：3级或4级数字层级
    parts = target_cath.split('.')
    if len(parts) not in (3, 4) or any((not p.isdigit()) for p in parts):
        print(f"错误：CATH编号格式不合法 '{target_cath}'，示例：1.10.10 或 1.10.150.10")
        return
    
    # 验证文件存在
    if not Path(excel_path).exists():
        print(f"错误：文件 '{excel_path}' 不存在")
        return
    
    # 1. 创建空字典存储内容
    domain_sequences = {}
    failed_records = []
    
    # 2. 读取Excel文件的第一列（跳过第一行表头）
    print(f"\n正在读取Excel文件: {excel_path}")
    try:
        df_input = pd.read_excel(excel_path, header=0)  # header=0表示第一行为列名
        uniprot_ids = [str(x).strip() for x in df_input.iloc[:, 0].dropna().tolist()]
        uniprot_ids = [x for x in uniprot_ids if x]
        print(f"读取到 {len(uniprot_ids)} 个UniProt ID")
    except Exception as e:
        print(f"读取Excel文件失败: {e}")
        return
    
    # 3-5. 遍历每个UniProt ID进行处理
    print(f"\n目标CATH编号: {target_cath}")
    print("-" * 50)
    
    for idx, uniprot_id in enumerate(uniprot_ids, 1):
        print(f"[{idx}/{len(uniprot_ids)}] 处理: {uniprot_id}")
        
        # 3. 访问CATH API
        cath_url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_id}"
        
        try:
            response = get_with_retries(cath_url)
            if response.status_code != 200:
                print(f"  ⚠️  CATH API失败: HTTP {response.status_code}，跳过")
                failed_records.append({
                    'uniprot_id': uniprot_id,
                    'stage': 'fetch_cath',
                    'reason': f'HTTP {response.status_code}'
                })
                continue
            
            cath_data = response.json()
            
            # 查找匹配目标CATH的结构域
            matching_domains = extract_cath_domains(cath_data, target_cath)
            
            if not matching_domains:
                print(f"  ⚠️  未找到匹配的CATH编号 {target_cath}，跳过")
                failed_records.append({
                    'uniprot_id': uniprot_id,
                    'stage': 'match_cath',
                    'reason': f'No matching cath_label={target_cath}'
                })
                continue
            
            # 4. 获取完整蛋白序列
            full_sequence = fetch_protein_sequence(uniprot_id)
            print(f"  ✅ 获取完整序列: {len(full_sequence)} 个残基")
            
            # 5. 为每个匹配的结构域创建键值对
            for dom_idx, domain in enumerate(matching_domains, 1):
                # 生成键名：UniProt ID + 序号
                key = f"{uniprot_id}-{dom_idx}"
                chopping = domain['chopping']

                if not chopping:
                    failed_records.append({
                        'uniprot_id': uniprot_id,
                        'stage': 'extract_domain',
                        'reason': f'Missing chopping for domain {dom_idx}'
                    })
                    print(f"    ⚠️  结构域 #{dom_idx} 缺少chopping，跳过")
                    continue
                
                # 截取结构域序列
                domain_seq = extract_domain_sequence(full_sequence, chopping)
                domain_sequences[key] = domain_seq
                
                print(f"    结构域 #{dom_idx}: chopping={chopping}, 序列长度={len(domain_seq)}")
            
            # 避免请求过快
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  ❌ 处理失败: {e}")
            failed_records.append({
                'uniprot_id': uniprot_id,
                'stage': 'pipeline',
                'reason': str(e)
            })
            continue
    
    # 6. 输出第一步结果（序列Excel）
    if not domain_sequences:
        print("\n错误：没有成功提取到任何结构域序列")
        return
    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    input_file_stem = Path(excel_path).stem
    output_seq_path = f"domain_sequences_{input_file_stem}_{timestamp}.xlsx"
    df_seq = pd.DataFrame(list(domain_sequences.items()), columns=['key', 'sequence'])
    df_seq.to_excel(output_seq_path, index=False)
    print(f"\n✅ 第一步完成！已输出序列文件: {output_seq_path}")
    print(f"   共提取 {len(domain_sequences)} 个结构域序列")
    
    # 7. 调用ProCaliper/Biopython计算理化性质
    print("\n开始计算理化性质...")
    print("-" * 50)
    
    all_properties = []
    
    for key, sequence in domain_sequences.items():
        print(f"计算: {key} (序列长度: {len(sequence)})")
        try:
            # 计算性质
            props = calculate_protein_properties(sequence)
            # 展平为记录
            flat_record = properties_to_flat_dict(key, sequence, props)
            all_properties.append(flat_record)
            print(f"  ✅ pI={props['isoelectric_point']}, MW={props['molecular_weight']:.1f}")
        except Exception as e:
            print(f"  ❌ 计算失败: {e}")
            failed_records.append({
                'uniprot_id': key,
                'stage': 'calculate_properties',
                'reason': str(e)
            })
            # 添加空记录以保持一致性
            all_properties.append({
                'key': key,
                'sequence': sequence,
                **{col: None for col in ['length', 'molecular_weight', 'isoelectric_point', 
                                         'aromaticity', 'instability_index', 'gravy',
                                         'extinction_coefficient_reduced', 'extinction_coefficient_oxidized',
                                         'helix_fraction', 'turn_fraction', 'sheet_fraction',
                                         'negatively_charged_count', 'positively_charged_count'] +
                        [f'aa_percent_{aa}' for aa in 'ACDEFGHIKLMNPQRSTVWY']}
            })
    
    # 写入最终结果Excel
    df_results = pd.DataFrame(all_properties)
    
    # 重新排列列顺序：key, sequence, 然后其他性质
    cols = ['key', 'sequence'] + [c for c in df_results.columns if c not in ['key', 'sequence']]
    df_results = df_results[cols]
    
    output_final_path = f"domain_properties_{input_file_stem}_{timestamp}.xlsx"
    df_results.to_excel(output_final_path, index=False)

    if failed_records:
        failed_path = f"domain_failed_records_{input_file_stem}_{timestamp}.xlsx"
        pd.DataFrame(failed_records).to_excel(failed_path, index=False)
    else:
        failed_path = None
    
    print(f"\n{'='*50}")
    print(f"🎉 批处理完成！")
    print(f"   最终结果文件: {output_final_path}")
    print(f"   共处理 {len(all_properties)} 个结构域序列")
    success_ratio = len(all_properties) / max(len(uniprot_ids), 1)
    print(f"   成功率（按至少产出1个结构域计）: {success_ratio:.2%}")
    if failed_path:
        print(f"   失败明细文件: {failed_path}（共 {len(failed_records)} 条）")
    print(f"\n输出文件包含以下列：")
    print("   - key: 结构域标识符（UniProt ID-序号）")
    print("   - sequence: 结构域序列（单字母表示）")
    print("   - length: 序列长度")
    print("   - molecular_weight: 分子量 (Da)")
    print("   - isoelectric_point: 等电点 (pI)")
    print("   - aromaticity: 芳香性指数")
    print("   - instability_index: 不稳定性指数")
    print("   - gravy: 亲水性指数 (GRAVY)")
    print("   - extinction_coefficient_*: 摩尔消光系数")
    print("   - helix/turn/sheet_fraction: 二级结构倾向分数")
    print("   - aa_percent_*: 20种标准氨基酸的百分比 (%)")


if __name__ == "__main__":
    main()