import numpy as np
from scipy import stats
import pandas as pd
from basefunctions import perm_test
import pandas as pd
from basefunctions import read_file_to_df
from basefunctions import calculate_association_traits


if __name__ == '__main__':

    # nealelab
    # columns: p2	p1	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se	r2p	description_p1	description_p2
    all_traits_file = r'C:\Users\hqy\Documents\Pumc\correlation_database\nealelab_UKBB\UKBB_ldsc_r2-master\r2_results\geno_correlation_sig.r2'
    breath_traits_file = r'C:\Users\hqy\Documents\Pumc\GWAS_catalog\Results\Associations\nealelab\respiratory_all_keyword_match\Two_trait_pairts_matched.xlsx'

    calculate_association_traits(all_traits_file=all_traits_file,
                                 breath_traits_file=breath_traits_file,
                                 thresholds=[0.2,  0.3, 0.4, 0.5, 0.6, 0.7],
                                 num_perms=[1, 5, 10, 50, 500],
                                 replace=True,
                                 output_file='Association_traits_significant_valuations_range_test.nealelab.xlsx',
                                 comparator='>',  # 筛选大于阈值的数据
                                 threshold_col='rg', absolute=True)

    # geneAtlas
    # 30280-0.0	Immature reticulocyte fraction	0.165565	21002-0.0	Weight	0.278915	0.27028756	0.255242467	0.277736415	0.222118419	0.025220884	0.196897535
    # Phenotype_1_ID	Phenotype_1_desc	Phenotype_1_h2	Phenotype_2_ID	Phenotype_2_desc	Phenotype_2_h2	r_y	r_g	r_e	cov_y	cov_g	cov_e
    # 23129-0.0	Trunk fat-free mass	0.319109	23116-0.0	Leg fat mass (left)	0.222917	0.518128174	0.532876419	0.516916135	2.492311834	0.420937187	2.071374647

    breath_traits_file = r'C:\Users\hqy\Documents\Pumc\GWAS_catalog\Results\Associations\geneATLAS\respiratory_all_keyword_match\Two_trait_pairts_matched.xlsx'
    all_traits_file = r'C:\Users\hqy\Documents\Pumc\correlation_database\geneAtlas\geneATLAS.corr.21_08_17.csv.gz'
    calculate_association_traits(all_traits_file=all_traits_file,
                                 breath_traits_file=breath_traits_file,
                                 thresholds=[0.2,  0.3, 0.4, 0.5, 0.6, 0.7],
                                 num_perms=[1, 5, 50, 500],
                                 replace=True,
                                 output_file='Association_traits_significant_valuations_range_test.geneAtlas.xlsx',
                                 comparator='>',  # 筛选大于阈值的数据
                                 threshold_col='r_g', absolute=True, confounding_cols=['cov_g'])
