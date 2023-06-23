import pandas as pd
import re
from impact_factor.core import Factor
from basefunctions import get_efo_by_keywords

# #######  modify ##########
out_pre_name = 'sleep_apnoea.'
# out_pre_name = 'Associations_filter_by_EFO_we_want.'
trait_mode = 'efo_trait'  # 'efo'  'efo_trait' 'article_trait', if None, need to specify trait_column
trait_column = None   # 'LINK' None
trait_keywords = [
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
trait_keywords = [
            'apnoea'
        ]
efo_list = [
            "EFO_0000341",
            "EFO_0006953",
            "EFO_0000768",
            "EFO_0007614",
            "EFO_1002011",
            "EFO_0009759",
            "EFO_0010638",
            "EFO_0010049",
            "EFO_0008590",
            "EFO_0009448",
            "EFO_0004244",
            "EFO_0006505",
            "MONDO_0100096",
            "EFO_0600020",
            "EFO_0600019",
            "EFO_0803362",
            "EFO_0001361",
            "EFO_0000770",
            "EFO_0009680",
            "EFO_1000485",
            "EFO_0009637",
            "EFO_1000362",
            "MONDO_0001437",
            "EFO_0005220",
            "EFO_0005853",
            "EFO_0006505",
            "EFO_0009661",
            "EFO_1002018",
            "EFO_0007184",
            "EFO_0003918",
            "EFO_0007817",
            "EFO_0008455",
            "EFO_0008456",
            "EFO_0003877",
            "EFO_0006527",
            "EFO_0004318",
            "EFO_0004319",
            "EFO_0005670",
            "EFO_0005671",
            "EFO_0008361",
            "EFO_0009115",
            "EFO_0021784",
            "EFO_0010339",
            "EFO_0004713",
            "MONDO_0008903"
        ]
# Specify the threshold value
p_value_threshold = 5E-8
impact_fact_threshold = 0
population_threshold = 10
association_file = 'gwas_catalog_v1.0.2-associations_e109_r2023-05-20.tsv'

# #################  end modify ##########

search_list = trait_keywords
if trait_mode == 'article_trait':
    trait_column = 'DISEASE/TRAIT'
elif trait_mode == 'efo_trait':
    trait_column = 'MAPPED_TRAIT'
elif trait_mode == 'efo':
    trait_column = 'MAPPED_TRAIT_URI'
    search_list = efo_list



associations_df = pd.read_csv(association_file, sep='\t', low_memory=False)
# Drop rows with missing values in the 'MAPPED_TRAIT_URI' column
associations_df = associations_df.dropna(subset=['MAPPED_TRAIT_URI'])

# Filter #
# Traits filter
search_pattern = '|'.join(search_list)
pattern = re.compile(search_pattern, flags=re.IGNORECASE)
associations_df = associations_df[associations_df[trait_column].str.contains(pattern, regex=True)]


def find_matched_items(x):
    items = []
    for item in search_list:
        if item in x:
            items.append(item)
    return ', '.join(items)

associations_df['matched_items'] = associations_df['MAPPED_TRAIT_URI'].apply(find_matched_items)

# P value filter
associations_df = associations_df[associations_df['P-VALUE'] <= p_value_threshold]

# population filter
# population filter Extract the maximum numbers from each value in the "INITIAL SAMPLE SIZE" column
numbers_pattern = r'\d{1,3}(?:,\d{3})*(?:\.\d+)?'


def extract_max_number(value):
    numbers = [int(num.replace(',', '')) for num in re.findall(numbers_pattern, value)]
    return max(numbers) if numbers else 0


max_numbers = associations_df["INITIAL SAMPLE SIZE"].apply(extract_max_number)

associations_df['Max_population_numbers'] = max_numbers
associations_df = associations_df[max_numbers >= population_threshold]

# impact factor filter
fa = Factor()


def search_factor(journal):
    try:
        factor = fa.search(journal)[0]['factor']
        return factor
    except IndexError:
        print('->The journal Error is :', journal)
        return 999  # Return a default value if the impact factor is not found


associations_df['JOURNAL_IMPACT_FACTOR'] = associations_df['JOURNAL'].apply(search_factor)

associations_df = associations_df[associations_df['JOURNAL_IMPACT_FACTOR'] >= impact_fact_threshold]


# Remove duplicates based on the specified column
associations_df = associations_df.drop_duplicates(subset='SNPS')

# SNP counts filter
columns = ['LINK', 'DISEASE/TRAIT']

#   Calculate the count of each combination of items across the specified columns
count_series = associations_df[columns].apply(lambda x: tuple(x), axis=1).value_counts()
associations_df['SNP counts'] = associations_df[columns].apply(lambda x: tuple(x), axis=1).map(count_series)

# save
associations_df.to_csv(f'{out_pre_name}{trait_mode}_mode.csv', index=False)




