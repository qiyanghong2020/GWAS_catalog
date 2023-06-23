import pandas as pd
from pandasgwas.get_studies import get_studies_by_efo_id
from pandasgwas import summary_statistics
import re
import os
import numpy as np
import gzip
import operator
from scipy import stats
import statsmodels.api as sm


respiratory_trait_keywords = [
            'respiratory',
            'lung',
            'trachea',
            'airway',
            'alveol',
            'breathing',
            'inhalation',
            'exhalation',
            'oxygen',
            'carbon',
            'asthma',
            'pulmonary',
            'bronch',
            'emphysema',
            'pneumonia',
            'tuberculosis',
            'pulmonary',
            'rhinitis',
            'sinusitis',
            'apnea',
            'pleural',
            'spirometry',
            'bronchoscopy',
            'ventilation',
            'intubation',
            'tracheostomy',
            'coughing',
            'wheezing',
            'breath',
            'sputum',
            'smok',
            'cigarettes',
            'snoring',
            'tobacco',
            ' air',
            'COVID',
            'Silicosis',
            'FEV',
            'FVC'
        ]


def read_file_to_df(filename):
    _, file_extension = os.path.splitext(filename)
    if file_extension == '.csv':
        return pd.read_csv(filename)
    elif file_extension == '.tsv':
        return pd.read_csv(filename, sep='\t')
    elif file_extension == '.r2':
        return pd.read_csv(filename, sep='\t')
    elif file_extension == '.gz':
        return pd.read_csv(filename, compression='gzip')
    elif file_extension == '.xlsx':
        return pd.read_excel(filename)
    else:
        raise ValueError(f'Unsupported file type {file_extension}')


def get_efo_by_keywords(keywords, save_name=None):
    # match all respiratory system traits in EFO, based on the related keywards
    efo_df = read_file_to_df('all_EFO_traits.csv')
    all_match_efo_df = pd.DataFrame()
    all_match_efo_df = efo_df[efo_df['trait'].str.contains('|'.join(keywords), regex=True)]

    # print(all_match_efo_df)
    if not save_name:
        save_name = '00.all_efo_items_related_with_respiratory_traits.csv'
    all_match_efo_df.to_csv(save_name, index=False)

    # After we get the respiratory traits in EFO, then we try to find the Study ID
    efo_ids = all_match_efo_df['shortForm']
    efo_list = list(efo_ids)
    return efo_list


def search_summary_statistics_and_download(accession_id, download=False):
    search_df = summary_statistics.search(study_accession_id=accession_id)
    print('\n-----\naccession_id=', accession_id)
    if search_df['file_name'].str.contains('tsv').any():
        print(search_df)
        if download:
            print('Wait for summary downloading...')
            summary_statistics.download(search_df)
    else:
        print('There are not summary statistics to download!')


def clean_string_for_linux_dir(s):
    # This regular expression pattern matches any character that is not a letter, a number, or a hyphen.
    pattern = r'[^a-zA-Z0-9\-_/]'
    # The sub() function replaces the matched characters with '_'
    return re.sub(pattern, '_', s)

def create_dirs_and_del_files(out_pre_path):
    # Check if directory exists, if not, create it
    if not os.path.exists(out_pre_path):
        os.makedirs(out_pre_path)

    # 遍历文件夹中的所有文件并删除
    for file_name in os.listdir(out_pre_path):
        file_path = os.path.join(out_pre_path, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

    print("文件夹已清空")


def comparing_pairwise(data_path, out_pre_path, trait1_col, trait2_col, traits, all=False):
    create_dirs_and_del_files(out_pre_path)

    traits_pairs_df = pd.DataFrame()
    # columns: p2	p1	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se	r2p	description_p1	description_p2
    df = read_file_to_df(data_path)
    for index, trait1 in enumerate(traits[:-1]):
        trait1_series = df[trait1_col].str.contains(trait1, regex=True, case=False)
        trait1_df = df[trait1_series].copy()
        for trait2 in traits[index + 1:]:
            trait2_series = trait1_df[trait2_col].str.contains(trait2, regex=True, case=False)
            trait2_df = trait1_df[trait2_series].copy()
            trait1_path = clean_string_for_linux_dir(trait1.split('|')[0])
            trait2_path = clean_string_for_linux_dir(trait2.split('|')[0])
            empty_flag = ''
            if trait2_df.empty:
                print("trait2_df is empty")
                empty_flag = 'Null.'
            else:
                print("trait2_df is not empty")
            trait2_df.to_excel(f'{out_pre_path}\{empty_flag}Association......{trait1_path}....{trait2_path}......xlsx',
                               index=False)
            if all:
                traits_pairs_df = pd.concat([traits_pairs_df, trait2_df])
            else:
                traits_pairs_df = pd.concat([traits_pairs_df, trait2_df[[trait1_col, trait2_col]]])
    traits_pairs_df.drop_duplicates().to_excel(f'{out_pre_path}/Two_trait_pairts_matched.xlsx', index=False)


def perm_test(all_traits_df, breath_traits_df, num_perm, replace=True):
    # 计算呼吸trait的显著性比例
    breath_ratio = np.sum(breath_traits_df['is_significant']) / len(breath_traits_df)

    # 进行排列试验，计算每次抽取的显著性比例
    ratios = []
    for _ in range(num_perm):
        sample_df = all_traits_df.sample(len(breath_traits_df), replace=replace)
        ratio = np.sum(sample_df['is_significant']) / len(sample_df)
        ratios.append(ratio)

    # 打印ratios
    print('Ratios: ', ratios)

    # 打印breath_ratio
    print('Breath ratio: ', breath_ratio)

    # 计算呼吸trait显著性在所有排列试验中的位置
    p_value = np.sum(np.array(ratios) >= breath_ratio) / num_perm

    return p_value


def calculate_association_traits(all_traits_file, breath_traits_file, thresholds, num_perms,
                                 replace, output_file, threshold_col='p', comparator='<', absolute=False,
                                 confounding_cols=('rg', 'se', 'z', 'h2_obs', 'h2_obs_se', 'h2_int',
                                                   'h2_int_se', 'gcov_int', 'gcov_int_se')):
    # 根据comparator参数决定使用的比较运算符
    ops = {'<': operator.lt, '<=': operator.le, '>': operator.gt, '>=': operator.ge, '==': operator.eq}
    if comparator not in ops:
        raise ValueError('Invalid comparator. Expected one of: {}'.format(ops.keys()))

    # 读取文件
    all_traits_df = read_file_to_df(all_traits_file)
    breath_traits_df = read_file_to_df(breath_traits_file)

    # 设置输出格式
    df_out = pd.DataFrame(columns=[f'{threshold_col}_threshold',
                                   'num_perm',
                                   'all_traits_rows',
                                   'all_traits_sig_rows',
                                   'all_traits_sig_ratio',
                                   'breath_traits_rows',
                                   'breath_traits_sig_rows',
                                   'breath_traits_ratio',
                                   'p_value_perm',
                                   'p_chi2',
                                   'p_value_binom'])

    for threshold in thresholds:
        for num_perm in num_perms:
            # 打印阈值、检验次数和抽样方式
            print(f"设定的显著性阈值为: {threshold}")
            print(f"设定的排列检验次数为: {num_perm}")
            print(f"设定的抽样方法为: {'有放回抽样' if replace else '不放回抽样'}")

            # 创建显著性列，p值小于阈值认为是显著的
            if absolute:
                all_traits_df['is_significant'] = ops[comparator](abs(all_traits_df[threshold_col]), threshold)
                breath_traits_df['is_significant'] = ops[comparator](abs(breath_traits_df[threshold_col]), threshold)
            else:
                all_traits_df['is_significant'] = ops[comparator](all_traits_df[threshold_col], threshold)
                breath_traits_df['is_significant'] = ops[comparator](breath_traits_df[threshold_col], threshold)


            # 获取行数并打印
            all_traits_rows = all_traits_df.shape[0]
            breath_traits_rows = breath_traits_df.shape[0]
            print(f"all_traits_df的总行数为: {all_traits_rows}")
            print(f"breath_traits_df的总行数为: {breath_traits_rows}")

            # 获取显著性的行数并打印
            all_traits_sig_rows = all_traits_df[all_traits_df['is_significant']].shape[0]
            breath_traits_sig_rows = breath_traits_df[breath_traits_df['is_significant']].shape[0]
            print(f"all_traits_df中显著性的行数为: {all_traits_sig_rows}")
            print(f"breath_traits_df中显著性的行数为: {breath_traits_sig_rows}")

            # 执行排列测试
            p_value_perm = perm_test(all_traits_df, breath_traits_df, num_perm=num_perm,
                                                                  replace=replace)

            # 进行卡方检验

            # 计算所有trait的显著性比例
            all_ratio = np.sum(all_traits_df['is_significant']) / len(all_traits_df)

            observed = [np.sum(breath_traits_df['is_significant']),
                        len(breath_traits_df) - np.sum(breath_traits_df['is_significant'])]
            expected = [all_ratio * len(breath_traits_df), (1 - all_ratio) * len(breath_traits_df)]
            chi2, p_chi2 = stats.chisquare(f_obs=observed, f_exp=expected)

            # 进行二项检验
            p_value_binom = stats.binom_test(np.sum(breath_traits_df['is_significant']), len(breath_traits_df),
                                             all_ratio)


            # 打印排列测试和二项测试的结果
            print("排列检验：", "拒绝零假设，认为呼吸相关的trait pair中显著的比例确实高于一般情况" if p_value_perm < 0.05 else "接受零假设，认为呼吸相关的trait pair中显著的比例并未高于一般情况")
            print("二项检验：", "拒绝零假设，认为呼吸相关的trait pair中显著的比例确实高于一般情况" if p_value_binom < 0.05 else "接受零假设，认为呼吸相关的trait pair中显著的比例并未高于一般情况")

            # 打印p值
            print("Permutation Test p-value: ", p_value_perm)
            print("Chi-square value: ", chi2)
            print("Chi-square Test p-value: ", p_chi2)
            print("Binomial Test p-value: ", p_value_binom)

            # 计算并保存显著性行数的比例
            all_traits_sig_ratio = (all_traits_sig_rows / all_traits_rows) * 100
            breath_traits_ratio = (breath_traits_sig_rows / breath_traits_rows) * 100

            # 保存计算结果到输出文件中
            row_dict = {
                f'{threshold_col}_threshold': threshold,
                'num_perm': num_perm,
                'all_traits_rows': all_traits_rows,
                'all_traits_sig_rows': all_traits_sig_rows,
                'all_traits_sig_ratio': "{:.2f}%".format(all_traits_sig_ratio),
                'breath_traits_rows': breath_traits_rows,
                'breath_traits_sig_rows': breath_traits_sig_rows,
                'breath_traits_ratio': "{:.2f}%".format(breath_traits_ratio),
                'p_value_perm': p_value_perm,
                'p_chi2': p_chi2,
                'p_value_binom': p_value_binom
            }


            # 查看混淆因素的影响
            # 假设 h2_obs 是一个在 all_traits_df 中的列
            # 使用逻辑回归模型来检查 h2_obs 对 'is_significant' 的影响
            for col in confounding_cols:
                all_traits_df[col] = all_traits_df[col].astype(float)
                X = sm.add_constant(all_traits_df[col])  # 添加常数项
                y = all_traits_df['is_significant']
                print('col=', col)
                # 拟合模型并查看结果
                model = sm.Logit(y, X)
                results = model.fit()
                print(results.summary())

                # 提取 p 值并添加到 row_dict
                p_value = results.pvalues[col]
                row_dict[f'confounding_{col}_pvalue'] = p_value

            # 使用concat()方法将新行添加到DataFrame中
            df_out = pd.concat([df_out, pd.DataFrame([row_dict])], ignore_index=True)

    df_out.to_excel(output_file, index=False)


