import os
import shutil

ABS_PATH = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import pandas as pd


def _store_resource_file(file_path, resource_file_name):
    '''
    Internal helper function to store a file in the resource directory.

    Parameters
    -----------
    :param file_path: the file that should be stored.
    :param resource_file_name: The file name in the resource folder.

    Examples
    -----------
    >>> _store_resource_file("/Users/username/Downloads/compounds.csv", "ebi_drugs.csv")
    '''
    resource_path = f"{ABS_PATH}/resources/"

    if not os.path.exists(resource_path):
        os.makedirs(resource_path)

    resource_file = os.path.join(resource_path, resource_file_name)
    shutil.copy(file_path, resource_file)


def store_ebi_drugs(file_path):
    '''
    Stores EBI drug names downloaded from https://www.ebi.ac.uk/chembl/g/#search_results/compounds
    as a CSV file in the local resources folder.

    Parameters
    -----------
    :param file_path: the file path at which the EBI compund file is located

    Examples
    -----------
    >>> store_ebi_drugs(file_path="/Users/Username/Downloads/compounds.csv")
    '''
    columns = ["Name", "Synonyms", "ChEMBL ID"]
    header = pd.read_csv(file_path, delimiter=';', nrows=0)

    if not all(column in header.columns for column in columns):
        raise ValueError(f"EBI drug names file must include: {columns}")

    _store_resource_file(file_path, "ebi_drugs.csv")


def load_ebi_drugs(file_path=f"{ABS_PATH}/resources/ebi_drugs.csv"):
    '''
    Loads and parses EBI drug names downloaded from https://www.ebi.ac.uk/chembl/g/#search_results/compounds.

    Parameters
    -----------
    :param file_path: the file path at which the EBI compund file is located
    :return: a tupel of drug names with associated chembl ids

    Examples
    -----------
    >>> drug_names, chembl_ids = load_ebi_drugs()
    '''
    df = pd.read_csv(file_path, delimiter=';', low_memory=False, usecols=["Name", "Synonyms", "ChEMBL ID"])

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


def store_drugbank_drugs(file_path):
    '''
    Stores DrugBank drug names downloaded from https://go.drugbank.com/releases/latest#open-data
    as a CSV file in the local resources folder.

    Parameters
    -----------
    :param file_path: the file path at which the DB compund file is located

    Examples
    -----------
    >>> store_drugbank_drugs(file_path="/Users/Username/Downloads/compounds.csv")
    '''
    columns = ["Common name", "Synonyms", "DrugBank ID"]
    header = pd.read_csv(file_path, delimiter=',', nrows=0)

    if not all(column in header.columns for column in columns):
        raise ValueError(f"DrugBank drug names file must include: {columns}")

    _store_resource_file(file_path, "drugbank_drugs.csv")


def load_drugbank_drugs(file_path=f"{ABS_PATH}/resources/drugbank_drugs.csv"):
    '''
    Loads and parses DrugBank drug names downloaded from https://go.drugbank.com/releases/latest#open-data.

    Parameters
    -----------
    :param file_path: the file path at which the DB compund file is located
    :return: a tupel of drug names with associated db ids

    Examples
    -----------
    >>> drug_names, db_ids = load_drugbank_drugs()
    '''
    df = pd.read_csv(file_path, delimiter=',', low_memory=False)

    drug_names, db_ids = [], []

    for _, row in df.iterrows():
        drug_names.append(row["Common name"])
        db_ids.append(row["DrugBank ID"])

        if row["Synonyms"] is np.nan:
            continue

        for synonym in row["Synonyms"].split(" | "):
            drug_names.append(synonym)
            db_ids.append(row["DrugBank ID"])

    return drug_names, db_ids


def load_charite_drug_goldstandard(file_path=f"{ABS_PATH}/resources/charite_drug_goldstandard.csv"):
    '''
    Loads our provided charite gold standard.

    Parameters
    -----------
    :param file_path: the file path at which the gold standard is located.
    :return: a tupel of drug names with associated chembl ids and drugbank ids

    Examples
    -----------
    >>> drug_names, chembl_ids = load_charite_drug_goldstandard()
    '''
    df = pd.read_csv(file_path, delimiter=';')
    df = df[df['drug_class'] == 'no']

    drug_names = []
    chembl_ids = []
    drugbank_ids = []

    for drug_name, df_group in df.groupby('treatment'):
        drug_names.append(drug_name)
        chembl_ids.append(df_group['chembl_id'].to_numpy().tolist())

        drugbank_id = df_group['drugbank_id'].to_numpy().tolist()
        drugbank_ids.append([None] if drugbank_id[0] is np.nan else drugbank_id)

    return drug_names, chembl_ids, drugbank_ids


def load_database_drug_goldstandard(file_path=f"{ABS_PATH}/resources/database_drug_goldstandard.csv"):
    '''
    Loads our provided database gold standard.

    Parameters
    -----------
    :param file_path: the file path at which the gold standard is located.
    :return: a tupel of drug names with associated chembl ids and drugbank ids

    Examples
    -----------
    >>> drug_names, chembl_ids = load_database_drug_goldstandard()
    '''
    df = pd.read_csv(file_path, delimiter=";")
    df = df[df['drug_class'] == 'no']

    drug_names = df['treatment'].to_numpy().tolist()
    chembl_ids = [[chembl_ids] for chembl_ids in df['chembl_id'].apply(lambda x: None if x != x else x)]
    drugbank_ids = [[drugbank_ids] for drugbank_ids in df['drugbank_id'].apply(lambda x: None if x != x else x)]

    return drug_names, chembl_ids, drugbank_ids


def load_ctg_drug_goldstandard(file_path=f"{ABS_PATH}/resources/ctg_drug_goldstandard.csv"):
    '''
    Loads our provided ctg gold standard.

    Parameters
    -----------
    :param file_path: the file path at which the gold standard is located.
    :return: a tupel of drug names with associated chembl ids and drugbank ids

    Examples
    -----------
    >>> drug_names, chembl_ids = load_ctg_drug_goldstandard()
    '''
    df = pd.read_csv(file_path, delimiter=";")

    drug_names = df['treatment'].to_numpy().tolist()
    chembl_ids = [[chembl_ids] for chembl_ids in df['chembl_id'].apply(lambda x: None if x != x else x)]
    drugbank_ids = [[drugbank_ids] for drugbank_ids in df['drugbank_id'].apply(lambda x: None if x != x else x)]

    return drug_names, chembl_ids, drugbank_ids