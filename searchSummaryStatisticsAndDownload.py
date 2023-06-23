from basefunctions import search_summary_statistics_and_download

accession_ids = [
    'GCST90255667',
    'GCST90134662',
    'GCST90255647',
    'GCST90018792',
    'GCST006409',
    'GCST90018863',
    'GCST011922',
    'GCST002810'
]

for accession_id in accession_ids:
    search_summary_statistics_and_download(accession_id=accession_id)