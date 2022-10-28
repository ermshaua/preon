import os, urllib.request, pronto
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

import pandas as pd
import daproli as dp

from lxml import etree


def download_do(file_path=f"{ABS_PATH}/resources/do.obo"):
    '''
    Downloads and stores the disease ontology.

    Parameters
    -----------
    :param file_path: the file path at which the disease ontology should be stored.

    Examples
    -----------
    >>> download_do()
    '''
    url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"

    with urllib.request.urlopen(url) as file:
        ot = file.read().decode('utf-8')

    with open(file_path, "w") as file:
        file.write(ot)


def download_mesh(file_path=f"{ABS_PATH}/resources/mesh.xml"):
    '''
    Downloads and stores the MeSH data.

    Parameters
    -----------
    :param file_path: the file path at which the mesh data should be stored.

    Examples
    -----------
    >>> download_mesh()
    '''
    url = "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2022.xml"

    with urllib.request.urlopen(url) as file:
        ot = file.read().decode('utf-8')

    with open(file_path, "w") as file:
        file.write(ot)


def load_do_cancers(file_path=f"{ABS_PATH}/resources/do.obo", expand_doids=False):
    '''
    Loads the disease ontology.

    Parameters
    -----------
    :param file_path: the file path at which the disease ontology is stored.
    :param expand_doids: a flag to decide whether cancer sub or superclasses should be considered.
    :return: a tupel of cancer types with associated doids.

    Examples
    -----------
    >>> cancer_types, doids = load_do_cancers()
    '''
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


def load_mesh_cancers(file_path=f"{ABS_PATH}/resources/mesh.xml"):
    '''
    Loads the mesh cancer descriptors.

    Parameters
    -----------
    :param file_path: the file path at which the mesh data is stored.
    :return: a tupel of mesh descriptors with associated mesh ids.

    Examples
    -----------
    >>> disease_types, mesh_ids = load_mesh_cancers()
    '''
    with open(file_path, 'r') as f:
        tree = etree.parse(f, etree.HTMLParser())

    disease_types, mesh_ids = [], []

    for desc_record in tree.getroot().xpath("//descriptorrecord"):
        # check that record is a neoplasm
        tree_entries = [entry.text for entry in desc_record.xpath("treenumberlist/treenumber")]
        if not any(entry.startswith("C04") for entry in tree_entries): continue

        disease_types.append(desc_record.xpath("descriptorname/string")[0].text)
        mesh_ids.append("MESH:" + desc_record.xpath("descriptorui")[0].text)

    return disease_types, mesh_ids


def download_or_load_do_cancers(file_path=f"{ABS_PATH}/resources/do_cancers.obo", expand_doids=False):
    '''
    Downloads or loads the disease ontology (depending whether it is stored or not).

    Parameters
    -----------
    :param file_path: the file path at which the disease ontology is/should be stored.
    :param expand_doids: a flag to decide whether cancer sub or superclasses should be considered.
    :return: a tupel of cancer types with associated doids.

    Examples
    -----------
    >>> cancer_types, doids = download_or_load_do_cancers()
    '''
    if not os.path.exists(file_path):
        download_do(file_path)

    return load_do_cancers(file_path=file_path, expand_doids=expand_doids)


def download_or_load_mesh_cancers(file_path=f"{ABS_PATH}/resources/mesh.xml"):
    '''
    Downloads or loads the mesh cancer data (depending whether it is stored or not).

    Parameters
    -----------
    :param file_path: the file path at which the mesh data is/should be stored.
    :return: a tupel of mesh descriptors with associated mesh ids.

    Examples
    -----------
    >>> disease_types, mesh_ids = download_or_load_mesh_cancers()
    '''
    if not os.path.exists(file_path):
        download_mesh(file_path)

    return load_mesh_cancers(file_path=file_path)


def load_do_flat_mapping(file_path=f"{ABS_PATH}/resources/do_cancers.obo"):
    '''
    Loads a dictionary that maps the first two layers from the cancer disease ontology
    to all of their sub classes.

    Parameters
    -----------
    :param file_path: the file path at which the disease ontology is stored.
    :return: a dictionary that maps the first two layers from the disease ontology
    to all of their sub classes.

    Examples
    -----------
    >>> flat_mapping = load_do_flat_mapping()
    '''
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
    '''
    Reduces the doids from the disease ontology to its first two cancer layers doids.

    Parameter
    -----------
    :param cancer_types: cancer types from the disease ontology
    :param doids: associated doids
    :param do_flat_mapping: the flat mapping to reduce the doids
    :return: a tupel of cancer types with associated doids.

    Examples
    -----------
    >>> cancer_types, doids = download_or_load_do_cancers()
    >>> do_flat_mapping = load_do_flat_mapping()
    >>> cancer_types, doids = apply_do_flat_mapping_to_ontology(cancer_types, doids, do_flat_mapping)
    '''
    doids = dp.map(lambda doid: do_flat_mapping[doid], doids)
    expand_cancer_types, expanded_doids = [], []

    for cancer_type, entries in zip(cancer_types, doids):
        for doid in entries:
            expand_cancer_types.append(cancer_type)
            expanded_doids.append(doid)

    return expand_cancer_types, expanded_doids


def apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping):
    '''
    Reduces the doids from the provided gold standards to the first two cancer layers
    from the disease ontology.

    Parameter
    -----------
    :param cancer_types: cancer types from a gold standard
    :param doids: associated doids
    :param do_flat_mapping: the flat mapping to reduce the doids
    :return: a tupel of cancer types with associated doids.

    Examples
    -----------
    >>> cancer_types, doids = load_database_cancer_goldstandard()
    >>> do_flat_mapping = load_do_flat_mapping()
    >>> cancer_types, doids = apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping)
    '''
    new_cancer_types, new_doids = [], []

    for cancer_type, entries in zip(cancer_types, doids):
        if entries != [None]:
            entries = dp.filter(lambda doid: doid in do_flat_mapping, entries, ret_type=list)
            entries = dp.map(lambda doid: do_flat_mapping[doid], entries, ret_type=list)
            entries = dp.flatten(entries, ret_type=list)

        if len(entries) > 0:
            new_cancer_types.append(cancer_type), new_doids.append(entries)

    return new_cancer_types, new_doids


def load_database_cancer_goldstandard(file_path=f"{ABS_PATH}/resources/database_cancer_goldstandard.csv", return_source=False):
    '''
    Loads our provided database gold standard.

    Parameters
    -----------
    :param file_path: the file path at which the gold standard is located.
    :param return_source: return the data base from which an entry originates.
    :return: a tupel of cancer types (and their sources) with associated doids and mesh ids.

    Examples
    -----------
    >>> cancer_types, doids, mesh_ids = load_database_cancer_goldstandard()
    '''
    df = pd.read_csv(file_path, sep=";")
    sources, cancer_types, doids, mesh_ids = [], [], [], []

    for _, (cancer_type, doid, source, mesh_id) in df.iterrows():
        if doid != doid:
            doid = [None]
        else:
            doid = doid.split(',')

        if mesh_id != mesh_id:
            mesh_id = [None]
        else:
            mesh_id = mesh_id.split(',')

        sources.append(source)
        cancer_types.append(cancer_type)
        doids.append(doid)
        mesh_ids.append(mesh_id)

    if return_source is True:
        return cancer_types, sources, doids, mesh_ids

    return cancer_types, doids, mesh_ids


def load_ncbi_cancer_goldstandard(file_path=f"{ABS_PATH}/resources/ncbi_cancer_goldstandard.csv"):
    '''
    Loads our provided ncbi gold standard.

    Parameters
    -----------
    :param file_path: the file path at which the gold standard is located.
    :return: a tupel of cancer types with associated doids and mesh ids.

    Examples
    -----------
    >>> cancer_types, doids, mesh_ids = load_ncbi_cancer_goldstandard()
    '''
    df = pd.read_csv(file_path, sep=";")

    cancer_types, dids, mids = [], [], []

    for ncbi_name, df_group in df.groupby("cancer"):
        cancer_types.append(ncbi_name)
        doids = df_group["doid"].unique().tolist()
        mesh_ids = df_group["mesh"].unique().tolist()

        if len(doids) == 1 and doids[0] != doids[0]:
            dids.append([None])
        else:
            dids.append(doids)

        if len(mesh_ids) == 1 and mesh_ids[0] != mesh_ids[0]:
            mids.append([None])
        else:
            mids.append(dp.flatten([_.split(',') for _ in mesh_ids if _ == _]))

    return cancer_types, dids, mids