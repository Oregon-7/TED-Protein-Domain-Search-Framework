import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import SelectFromModel, RFECV
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import warnings
warnings.filterwarnings('ignore')

class ProteinExpressionAnalyzer:
    def __init__(self, file_path):
        """
        初始化分析器
        Args:
            file_path: Excel文件路径
        """
        self.df = pd.read_excel(file_path, index_col=0)
        self.protein_ids = self.df.index.tolist()
        self.species = self.df.columns.tolist()
        self.scaled_data = None
        self.feature_importance = None
        self.selected_proteins = None
        
    def normalize_data(self, method='robust'):
        """
        数据标准化
        Args:
            method: 标准化方法，可选 'standard', 'robust', 'minmax'
        """
        if method == 'standard':
            scaler = StandardScaler()
        elif method == 'robust':
            scaler = RobustScaler()  # 对异常值更鲁棒
        else:
            scaler = StandardScaler()
            
        self.scaled_data = pd.DataFrame(
            scaler.fit_transform(self.df.T).T,
            index=self.df.index,
            columns=self.df.columns
        )
        print(f"数据已使用 {method} 方法标准化")
        return self.scaled_data
    
    def get_user_labels(self):
        """
        获取用户输入的正类和负类标签
        """
        print("\n" + "="*50)
        print("可用的物种列表:")
        for i, species in enumerate(self.species, 1):
            print(f"{i}. {species}")
        
        print("\n请选择正类物种（高表达组）:")
        pos_indices = input("输入物种编号（多个用逗号分隔）: ").strip().split(',')
        pos_indices = [int(i.strip()) - 1 for i in pos_indices if i.strip().isdigit()]
        positive_species = [self.species[i] for i in pos_indices if i < len(self.species)]
        
        print("\n请选择负类物种（低表达/对照组）:")
        neg_indices = input("输入物种编号（多个用逗号分隔）: ").strip().split(',')
        neg_indices = [int(i.strip()) - 1 for i in neg_indices if i.strip().isdigit()]
        negative_species = [self.species[i] for i in neg_indices if i < len(self.species)]
        
        print(f"\n正类物种: {positive_species}")
        print(f"负类物种: {negative_species}")
        
        return positive_species, negative_species
    
    def create_labels(self, positive_species, negative_species):
        """
        创建二分类标签
        """
        y = []
        X_data = []
        
        # 合并所有物种数据
        all_species = positive_species + negative_species
        
        for species in all_species:
            if species in self.scaled_data.columns:
                X_data.append(self.scaled_data[species].values)
                if species in positive_species:
                    y.append(1)  # 正类
                else:
                    y.append(0)  # 负类
        
        X = np.array(X_data).T  # 转置为 (蛋白质数 × 物种数)
        y = np.array(y)
        
        return X, y, all_species
    
    def statistical_screening(self, X, y, alpha=0.05):
        """
        使用统计检验进行初步筛选
        """
        pos_indices = np.where(y == 1)[0]
        neg_indices = np.where(y == 0)[0]
        
        p_values = []
        fold_changes = []
        
        for i in range(X.shape[0]):
            pos_data = X[i, pos_indices]
            neg_data = X[i, neg_indices]
            
            # 计算fold change
            fc = np.mean(pos_data) - np.mean(neg_data)
            fold_changes.append(fc)
            
            # t检验
            if len(pos_data) > 1 and len(neg_data) > 1:
                t_stat, p_val = stats.ttest_ind(pos_data, neg_data, equal_var=False)
                p_values.append(p_val)
            else:
                p_values.append(1.0)  # 样本不足，设为不显著
        
        # FDR校正
        reject, pvals_corrected = fdrcorrection(p_values, alpha=alpha)
        
        # 筛选条件：显著且fold change为正
        significant_idx = np.where((reject) & (np.array(fold_changes) > 0))[0]
        
        return significant_idx, pvals_corrected, fold_changes
    
    def machine_learning_screening(self, X, y, method='random_forest', 
                                   feature_selection='importance'):
        """
        使用机器学习方法筛选重要蛋白质
        """
        if method == 'random_forest':
            model = RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=42
            )
        elif method == 'svm':
            model = SVC(kernel='linear', probability=True, random_state=42)
        else:
            model = RandomForestClassifier(n_estimators=100, random_state=42)
        
        # 划分训练测试集
        X_train, X_test, y_train, y_test = train_test_split(
            X.T, y, test_size=0.3, random_state=42, stratify=y
        )
        
        # 训练模型
        model.fit(X_train, y_train)
        
        # 评估模型
        y_pred = model.predict(X_test)
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        
        print("\n模型性能评估:")
        print(classification_report(y_test, y_pred))
        print(f"AUC-ROC: {roc_auc_score(y_test, y_pred_proba):.3f}")
        
        # 获取特征重要性
        if method == 'random_forest':
            importance = model.feature_importances_
        elif method == 'svm':
            importance = np.abs(model.coef_[0])  # SVM系数的绝对值
        else:
            importance = model.feature_importances_
        
        self.feature_importance = pd.DataFrame({
            'Protein_ID': self.protein_ids,
            'Importance': importance
        }).sort_values('Importance', ascending=False)
        
        # 特征选择
        if feature_selection == 'importance':
            # 选择重要性前N的特征
            threshold = np.percentile(importance, 75)  # 前25%
            selected_idx = np.where(importance > threshold)[0]
            
        elif feature_selection == 'recursive':
            # 递归特征消除
            selector = RFECV(model, cv=5)
            selector.fit(X_train, y_train)
            selected_idx = np.where(selector.support_)[0]
            
        elif feature_selection == 'from_model':
            # 基于模型的特征选择
            selector = SelectFromModel(model, threshold='median')
            selector.fit(X_train, y_train)
            selected_idx = np.where(selector.get_support())[0]
        
        return selected_idx, model
    
    def integrated_selection(self, positive_species, negative_species, 
                            method='random_forest', top_n=50):
        """
        综合统计和机器学习的集成选择方法
        """
        # 数据标准化
        if self.scaled_data is None:
            self.normalize_data()
        
        # 创建标签和数据
        X, y, all_species = self.create_labels(positive_species, negative_species)
        
        print(f"\n数据维度: {X.shape}")
        print(f"样本分布 - 正类: {sum(y==1)}, 负类: {sum(y==0)}")
        
        # 1. 统计筛选
        print("\n[阶段1] 进行统计检验筛选...")
        stat_idx, pvals, fold_changes = self.statistical_screening(X, y)
        print(f"统计筛选得到 {len(stat_idx)} 个显著差异蛋白质")
        
        # 2. 机器学习筛选
        print("\n[阶段2] 进行机器学习筛选...")
        ml_idx, model = self.machine_learning_screening(
            X, y, method=method, feature_selection='importance'
        )
        print(f"机器学习筛选得到 {len(ml_idx)} 个重要蛋白质")
        
        # 3. 集成选择（取交集或并集）
        # 这里我们取并集，但优先选择在两种方法中都显著的
        all_idx = list(set(stat_idx) | set(ml_idx))
        
        # 计算综合评分
        scores = []
        for idx in all_idx:
            # 结合p值、fold change和机器学习重要性
            if idx < len(pvals):
                p_score = -np.log10(pvals[idx] + 1e-10)  # p值负对数
            else:
                p_score = 0
                
            if idx < len(fold_changes):
                fc_score = fold_changes[idx]
            else:
                fc_score = 0
                
            # 查找机器学习重要性
            protein_id = self.protein_ids[idx]
            ml_row = self.feature_importance[self.feature_importance['Protein_ID'] == protein_id]
            if not ml_row.empty:
                imp_score = ml_row['Importance'].values[0]
            else:
                imp_score = 0
                
            # 综合评分（可根据需要调整权重）
            combined_score = p_score * 0.3 + fc_score * 0.3 + imp_score * 0.4
            scores.append(combined_score)
        
        # 创建结果DataFrame
        results = []
        for idx, score in zip(all_idx, scores):
            protein_id = self.protein_ids[idx]
            
            # 计算正类负类的平均表达量
            pos_data = X[idx, np.where(y == 1)[0]]
            neg_data = X[idx, np.where(y == 0)[0]]
            
            pos_mean = np.mean(pos_data)
            neg_mean = np.mean(neg_data)
            fold_change = pos_mean - neg_mean
            
            # 检查是否在统计筛选中显著
            is_stat_sig = idx in stat_idx
            
            results.append({
                'Protein_ID': protein_id,
                'Combined_Score': score,
                'Fold_Change': fold_change,
                'Positive_Mean': pos_mean,
                'Negative_Mean': neg_mean,
                'Is_Statistically_Significant': is_stat_sig
            })
        
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('Combined_Score', ascending=False)
        
        # 选择Top N
        self.selected_proteins = results_df.head(top_n)
        
        return self.selected_proteins, results_df
    
    def save_results(self, output_file='selected_proteins.xlsx'):
        """
        保存结果到Excel文件
        """
        if self.selected_proteins is not None:
            # 保存筛选出的蛋白质
            self.selected_proteins.to_excel(output_file, index=False)
            
            # 保存所有蛋白质的重要性排名
            importance_file = output_file.replace('.xlsx', '_importance.xlsx')
            self.feature_importance.to_excel(importance_file, index=False)
            
            print(f"\n结果已保存:")
            print(f"- 筛选出的蛋白质: {output_file}")
            print(f"- 所有蛋白质重要性: {importance_file}")
        else:
            print("没有筛选出的蛋白质可保存")
    
    def visualize_top_proteins(self, top_n=20):
        """
        可视化Top蛋白质的表达模式（可选）
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            if self.selected_proteins is not None:
                top_ids = self.selected_proteins.head(top_n)['Protein_ID'].tolist()
                
                # 获取原始数据
                top_data = self.df.loc[top_ids]
                
                # 绘制热图
                plt.figure(figsize=(12, 8))
                sns.heatmap(top_data, cmap='RdBu_r', center=0,
                           cbar_kws={'label': 'Expression Level'})
                plt.title(f'Top {top_n} Proteins Expression Pattern')
                plt.xlabel('Species')
                plt.ylabel('Protein ID')
                plt.tight_layout()
                plt.savefig('protein_expression_heatmap.png', dpi=300)
                plt.show()
                
                print(f"热图已保存为 'protein_expression_heatmap.png'")
                
        except ImportError:
            print("需要安装matplotlib和seaborn进行可视化")
        except Exception as e:
            print(f"可视化时出错: {e}")

def main():
    """
    主函数：交互式分析流程
    """
    print("="*60)
    print("蛋白质表达差异分析工具")
    print("="*60)
    
    # 1. 读取数据
    file_path = input("请输入Excel文件路径: ").strip()
    
    try:
        analyzer = ProteinExpressionAnalyzer(file_path)
        print(f"成功读取数据: {analyzer.df.shape[0]} 个蛋白质, {analyzer.df.shape[1]} 个物种")
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return
    
    # 2. 数据标准化
    print("\n选择标准化方法:")
    print("1. StandardScaler (标准正态分布)")
    print("2. RobustScaler (鲁棒标准化，推荐)")
    
    norm_choice = input("请选择 (默认2): ").strip()
    if norm_choice == '1':
        analyzer.normalize_data(method='standard')
    else:
        analyzer.normalize_data(method='robust')
    
    # 3. 用户选择正类负类
    positive_species, negative_species = analyzer.get_user_labels()
    
    if not positive_species or not negative_species:
        print("错误: 必须至少选择一个正类和一个负类物种")
        return
    
    # 4. 选择分析方法
    print("\n选择机器学习方法:")
    print("1. Random Forest (随机森林)")
    print("2. SVM (支持向量机)")
    
    ml_choice = input("请选择 (默认1): ").strip()
    method = 'random_forest' if ml_choice != '2' else 'svm'
    
    # 5. 选择要输出的蛋白质数量
    top_n = input(f"\n要输出的蛋白质数量 (默认50): ").strip()
    try:
        top_n = int(top_n) if top_n else 50
    except:
        top_n = 50
    
    # 6. 执行分析
    print("\n开始分析...")
    selected_proteins, all_results = analyzer.integrated_selection(
        positive_species, negative_species,
        method=method,
        top_n=top_n
    )
    
    # 7. 显示结果
    print("\n" + "="*60)
    print("筛选结果:")
    print("="*60)
    print(f"共筛选出 {len(selected_proteins)} 个蛋白质")
    print("\nTop 10 蛋白质:")
    print(selected_proteins.head(10).to_string())
    
    # 8. 保存结果
    save_option = input("\n是否保存结果? (y/n, 默认y): ").strip().lower()
    if save_option != 'n':
        output_file = input("输出文件名 (默认: selected_proteins.xlsx): ").strip()
        if not output_file:
            output_file = 'selected_proteins.xlsx'
        analyzer.save_results(output_file)
    
    # 9. 可视化选项
    viz_option = input("\n是否生成可视化图表? (y/n, 默认n): ").strip().lower()
    if viz_option == 'y':
        analyzer.visualize_top_proteins(min(20, top_n))
    
    print("\n分析完成!")

if __name__ == "__main__":
    main()