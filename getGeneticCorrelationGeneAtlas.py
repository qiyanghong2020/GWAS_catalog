from basefunctions import respiratory_trait_keywords
from basefunctions import comparing_pairwise

data_path = r'C:\Users\hqy\Documents\Pumc\correlation_database\geneATLAS.corr.21_08_17.csv.gz'
out_pre_path = r'C:\Users\hqy\Documents\Pumc\GWAS_catalog\Results\Associations\geneATLAS'
trait1_col = 'Phenotype_1_desc'
trait2_col = 'Phenotype_2_desc'
traits = [
    # '[a-z]',
    "asthma|bronchial asthma",
    "chronic obstructive pulmonary disease|COPD|chronic obstructive lung disease",
    "Severe COVID-19 infection|COVID|Coronavirus disease 2019|SARS-CoV-2 infection",
    "lung adenocarcinoma|squamous cell lung carcinoma|small cell lung carcinoma|lung carcinoma|pulmonary carcinoma|lung cancer",
    "idiopathic pulmonary fibrosis|IPF|cryptogenic fibrosing alveolitis",
    "Eczema|atopic dermatitis|atopic eczema",
    "allergic rhinitis|hay fever",
    "interstitial lung disease|diffuse parenchymal lung disease",
    "apnea|apnoea|sleep-disordered breathing",
    "tuberculosis|phthisis|consumption",
]
traits = respiratory_trait_keywords

# 30280-0.0	Immature reticulocyte fraction	0.165565	21002-0.0	Weight	0.278915	0.27028756	0.255242467	0.277736415	0.222118419	0.025220884	0.196897535
# Phenotype_1_ID	Phenotype_1_desc	Phenotype_1_h2	Phenotype_2_ID	Phenotype_2_desc	Phenotype_2_h2	r_y	r_g	r_e	cov_y	cov_g	cov_e
# 23129-0.0	Trunk fat-free mass	0.319109	23116-0.0	Leg fat mass (left)	0.222917	0.518128174	0.532876419	0.516916135	2.492311834	0.420937187	2.071374647


# 根据关键词两两比对搜索匹配项目
def main():
    comparing_pairwise(data_path, out_pre_path, trait1_col, trait2_col, traits, all=True)


if __name__ == '__main__':
    main()
