import os
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

import pandas as pd


def download_ebi_drugs():
    pass


def load_ebi_drugs(file_path=f"{ABS_PATH}/resources/ebi_drugs.tsv"):
    df = pd.read_csv(file_path, sep='\t')
    return df.name.tolist(), df.chembl_id.tolist()


def load_charite_drug_goldstandard(file_path=f"{ABS_PATH}/resources/charite_drug_goldstandard.csv"):
    df = pd.read_csv(file_path)
    df = df[df['drug_class'] == 'no']

    drug_names = []
    chembl_ids = []

    for drug_name, df_group in df.groupby('treatment'):
        drug_names.append(drug_name)
        chembl_ids.append(df_group['chembl_id'].to_numpy().tolist())

    return drug_names, chembl_ids


def load_database_drug_goldstandard(file_path=f"{ABS_PATH}/resources/database_drug_goldstandard.csv"):
    df = pd.read_csv(file_path)
    df = df[df['Drug Class'] == 'no']

    drug_names = df['Drug Name'].to_numpy().tolist()
    chembl_ids = [[chembl_ids] for chembl_ids in df['Chembl ID'].apply(lambda x: None if x != x else x)]

    return drug_names, chembl_ids


def load_ctg_drug_goldstandard(file_path=f"{ABS_PATH}/resources/ctg_drug_goldstandard.tsv"):
    df = pd.read_csv(file_path, sep='\t')

    drug_names = df['intervention'].to_numpy().tolist()
    chembl_ids = [[chembl_ids] for chembl_ids in df['chembl_id'].apply(lambda x: None if x != x else x)]

    return drug_names, chembl_ids