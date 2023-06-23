import pandas as pd
from basefunctions import clean_string_for_linux_dir

df = pd.read_excel(r'C:\Users\hqy\Documents\Pumc\Table1_specific_disorder\Download\manifest_GBMI_summary_statistics.slimfordownload.20230608.xlsx')

df_columns = ["phenotype", "phenotype_short", "sex", "ancestry", "Note", "biobank",
     "Genome-wide summary statistics (.gz)", "index file (.gz.tbi)",
     "QQ plot (.png)", "Manhattan plot (.png)"]


def format_bash_download_command(series):
    phenotype = series['phenotype']
    ancestry = series['ancestry']
    sex = series['sex']
    dir_path = f'/data/hongqy/database/GBMI/{phenotype}/{ancestry}/{sex}'
    dir_path = clean_string_for_linux_dir(dir_path)
    bash_command = f"""mkdir -p {dir_path} || true && cd {dir_path}\n
{series["Genome-wide summary statistics (.gz)"]}\n
{series["index file (.gz.tbi)"]}\n
{series["QQ plot (.png)"]}\n
{series["Manhattan plot (.png)"]}\n
"""
    return bash_command


out_serises = df.apply(format_bash_download_command, axis=1)
out_str = '\n'.join(out_serises)
print(out_str)
with open('GBMI_database_download.sh', 'w') as f:
    f.write(out_str)