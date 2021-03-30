import os
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

import pandas as pd


def load_ebi_drugs(file_path=f"{ABS_PATH}/resources/ebi_drugs.csv"):
    df = pd.read_csv(file_path, delimiter=';', low_memory=False)

    # reduce size
    df = df[["Name", "Synonyms", "ChEMBL ID"]]

    # filter df
    df = df[df.Name.notna() & df["ChEMBL ID"].notna()]

    drug_names, chembl_ids = [], []

    for idx, (name, synonyms, chembl_id) in df.iterrows():
        names = []

        # check for empty, add name
        if len(name) > 0:
            names.append(name)

        # check for nan/empty, add synonyms
        if isinstance(synonyms, str) and len(synonyms) > 0:
            for s in synonyms.split("|"):
                names.append(s)

        # check for empty, add chembl id
        if len(chembl_id) > 0:
            drug_names.extend(names)
            chembl_ids.extend([chembl_id] * len(names))

    return drug_names, chembl_ids


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
    df = df[df['drug_class'] == 'no']

    drug_names = df['treatment'].to_numpy().tolist()
    chembl_ids = [[chembl_ids] for chembl_ids in df['chembl_id'].apply(lambda x: None if x != x else x)]

    return drug_names, chembl_ids


def load_ctg_drug_goldstandard(file_path=f"{ABS_PATH}/resources/ctg_drug_goldstandard.csv"):
    df = pd.read_csv(file_path)

    drug_names = df['treatment'].to_numpy().tolist()
    chembl_ids = [[chembl_ids] for chembl_ids in df['chembl_id'].apply(lambda x: None if x != x else x)]

    return drug_names, chembl_ids