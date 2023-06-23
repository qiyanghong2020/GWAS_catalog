from pandasgwas.get_traits import get_traits_all
import pandas as pd
from pandasgwas.get_studies import get_studies_by_efo_id
from pandasgwas import summary_statistics
import sys


#traits = get_traits_all()
#traits.efo_traits.to_csv('all_EFO_traits.csv', index=False)

respiratory_terms = [
    'respiratory',
    'lung',
    'trachea',
    'airways',
    'alveoli',
    'breathing',
    'inhalation',
    'exhalation',
    'oxygen',
    'carbon',
    'asthma',
    'pulmonary',
    'bronchitis',
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
    'smoking',
    'cigarettes',
    'snoring',
    'tobacco',
    'smoke',
    ' air'
]



# match all respiratory system traits in EFO, based on the related keywards
efo_df = pd.read_csv('all_EFO_traits.csv')
all_match_efo_df = pd.DataFrame()
for respiratory_term in respiratory_terms:
    match_efo_df = efo_df[efo_df['trait'].str.contains(respiratory_term, case=False)]
    all_match_efo_df = pd.concat([all_match_efo_df, match_efo_df], ignore_index=True)
all_match_efo_df.drop_duplicates()
#print(all_match_efo_df)
all_match_efo_df.to_csv('00.all_efo_items_related_with_respiratory_traits.csv', index=False)

# After we get the respiratory traits in EFO, then we try to find the Study ID

efo_ids = all_match_efo_df['shortForm']
all_studies_df = pd.DataFrame()
c=0
for efo_id in efo_ids:
    c += 1
    print(c)
    print(efo_id)
    studies = get_studies_by_efo_id(efo_id)
    all_studies_df = pd.concat([all_studies_df, studies.studies], ignore_index=True)
    all_studies_df.to_csv('01.all_studies_related_with_respiratory_traits.csv', index=False)
all_studies_df.drop_duplicates(subset='accessionId')
all_studies_df.to_csv('01.all_studies_related_with_respiratory_traits.csv', index=False)

# get summary statistics based on study_accession_id.
all_studies_with_summary_statistics_df = pd.DataFrame()
for index, row in all_studies_df.iterrows():
    search_DF = summary_statistics.search(study_accession_id=row['accessionId'])
    if search_DF['file_name'].str.contains('tsv').any():
        row_df = pd.DataFrame([row])  # 将行转换为DataFrame对象
        all_studies_with_summary_statistics_df = pd.concat([all_studies_with_summary_statistics_df, row_df])
        #summary_statistics.download(search_DF)

all_studies_with_summary_statistics_df.to_csv('02.all_summary_statistics_related_with_respiratory_traits.csv', index=False)






