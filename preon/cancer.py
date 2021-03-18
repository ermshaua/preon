import os
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

import preon as _pronto
import numpy as np
import pandas as pd


def download_do_cancers():
    pass


def load_do_cancers(file_path=f"{ABS_PATH}/resources/do_cancers.csv", expand_doids=True):
    ot = _pronto.Ontology(file_path)

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