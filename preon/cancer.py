import os, urllib.request, pronto
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

import pandas as pd
import daproli as dp


def download_do_cancers(file_path=f"{ABS_PATH}/resources/do_cancers.obo"):
    url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"

    with urllib.request.urlopen(url) as file:
        ot = file.read().decode('utf-8')

    with open(file_path, "w") as file:
        file.write(ot)


def load_do_cancers(file_path=f"{ABS_PATH}/resources/do_cancers.obo", expand_doids=False):
    ot = pronto.Ontology(file_path)

    cancer_types = []
    doids = []

    cancer_ids = [term.id for term in ot["DOID:162"].subclasses()]

    for cancer_id in cancer_ids:
        if cancer_id == "DOID:162": continue

        term = ot[cancer_id]
        names = [term.name]

        # add doid synonyms
        for synonym in term.synonyms:
            if synonym.scope == "EXACT":
                names.append(synonym.description)

        names = set(names)

        for name in names:
            cancer_types.append(name)
            doids.append(term.id)

        if expand_doids is False:
            continue

        # add subclasses
        for subclass in term.subclasses():
            for name in names:
                cancer_types.append(name)
                doids.append(subclass.id)

        # add superclasses
        for superclass in term.superclasses():
            for name in names:
                cancer_types.append(name)
                doids.append(superclass.id)

    return cancer_types, doids


def download_or_load_do_cancers(file_path=f"{ABS_PATH}/resources/do_cancers.obo", expand_doids=False):
    if not os.path.exists(file_path):
        download_do_cancers(file_path)

    return load_do_cancers(file_path=file_path, expand_doids=expand_doids)


def load_do_flat_mapping(file_path=f"{ABS_PATH}/resources/do_cancers.obo"):
    ot = pronto.Ontology(file_path)
    mapping = dict()

    start_ids = ["DOID:0050687", "DOID:0050686"]

    for start_id in start_ids:
        for term in ot[start_id].subclasses(distance=1):
            for sub_term in term.subclasses():
                if sub_term.id not in mapping: mapping[sub_term.id] = list()
                mapping[sub_term.id].append(term.id)

    return mapping


def apply_do_flat_mapping_to_ontology(cancer_types, doids, do_flat_mapping):
    doids = dp.map(lambda doid: do_flat_mapping[doid], doids)
    expand_cancer_types, expanded_doids = [], []

    for cancer_type, entries in zip(cancer_types, doids):
        for doid in entries:
            expand_cancer_types.append(cancer_type)
            expanded_doids.append(doid)

    return expand_cancer_types, expanded_doids


def apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping):
    new_cancer_types, new_doids = [], []

    for cancer_type, entries in zip(cancer_types, doids):
        if entries != [None]:
            entries = dp.filter(lambda doid: doid in do_flat_mapping, entries, ret_type=list)
            entries = dp.map(lambda doid: do_flat_mapping[doid], entries, ret_type=list)
            entries = dp.flatten(entries, ret_type=list)

        if len(entries) > 0:
            new_cancer_types.append(cancer_type), new_doids.append(entries)

    return new_cancer_types, new_doids


def load_database_cancer_goldstandard(file_path=f"{ABS_PATH}/resources/database_cancer_goldstandard.csv"):
    df = pd.read_csv(file_path, sep=';')
    sources, cancer_types, doids = [], [], []

    for _, (source, cancer_type, doid) in df.iterrows():
        if doid != doid:
            doid = [None]
        else:
            doid = doid.split(',')

        sources.append(source)
        cancer_types.append(cancer_type)
        doids.append(doid)

    return cancer_types, doids


def load_ncbi_cancer_goldstandard(file_path=f"{ABS_PATH}/resources/ncbi_cancer_goldstandard.csv"):
    df = pd.read_csv(file_path)

    cancer_types = []
    dids = []

    for ncbi_name, df_group in df.groupby("ncbi_name"):
        cancer_types.append(ncbi_name)
        dids.append(df_group["disease_id"].unique().tolist())

    return cancer_types, dids