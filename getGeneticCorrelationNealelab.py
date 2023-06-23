import pandas as pd
import os
from basefunctions import clean_string_for_linux_dir
from basefunctions import respiratory_trait_keywords
from basefunctions import create_dirs_and_del_files
from basefunctions import comparing_pairwise

data_path = r'C:\Users\hqy\Documents\Pumc\correlation_database\nealelab_UKBB\UKBB_ldsc_r2-master\r2_results\geno_correlation_sig.r2'
out_pre_path = r'C:\Users\hqy\Documents\Pumc\GWAS_catalog\Results\Associations\nealelab'
trait1_col = 'description_p1'
trait2_col = 'description_p2'
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


def main():
    comparing_pairwise(data_path, out_pre_path, trait1_col, trait2_col, traits, all=True)


if __name__ == '__main__':
    main()



