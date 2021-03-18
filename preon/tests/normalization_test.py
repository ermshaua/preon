import unittest

from preon.normalization import PrecisionOncologyNormalizer
from preon.drug import load_ebi_drugs, load_charite_drug_goldstandard, load_database_drug_goldstandard, load_ctg_drug_goldstandard
from preon.cancer import load_do_cancers, load_do_flat_mapping, apply_do_flat_mapping_to_ontology, apply_do_flat_mapping_to_goldstandard, \
    load_database_cancer_goldstandard, load_ncbi_cancer_goldstandard
from preon.tests.utils import f1_score


class DrugNormalizationTest(unittest.TestCase):

    def test_charite_drug_goldstandard(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer().fit(drug_names, chembl_ids)

        drug_names, chembl_ids = load_charite_drug_goldstandard()
        df_eval = normalizer.evaluate(drug_names, chembl_ids)

        assert f1_score(df_eval) > .9

    def test_databse_drug_goldstandard(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer().fit(drug_names, chembl_ids)

        drug_names, chembl_ids = load_database_drug_goldstandard()
        df_eval = normalizer.evaluate(drug_names, chembl_ids)

        assert f1_score(df_eval) > .9

    def test_ctg_drug_goldstandard(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer().fit(drug_names, chembl_ids)

        drug_names, chembl_ids = load_ctg_drug_goldstandard()
        df_eval = normalizer.evaluate(drug_names, chembl_ids)

        assert f1_score(df_eval) > .9


class CancerNormalizationTest(unittest.TestCase):

    def test_database_cancer_goldstandard(self):
        cancer_types, doids = load_do_cancers()
        do_flat_mapping = load_do_flat_mapping()

        cancer_types, doids = apply_do_flat_mapping_to_ontology(cancer_types, doids, do_flat_mapping)
        normalizer = PrecisionOncologyNormalizer().fit(cancer_types, doids)

        cancer_types, doids = load_database_cancer_goldstandard()
        cancer_types, doids = apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping)
        df_eval = normalizer.evaluate(cancer_types, doids, n_grams=3)

        assert f1_score(df_eval) > .9

    def test_ncbi_cancer_goldstandard(self):
        cancer_types, doids = load_do_cancers()
        do_flat_mapping = load_do_flat_mapping()

        cancer_types, doids = apply_do_flat_mapping_to_ontology(cancer_types, doids, do_flat_mapping)
        normalizer = PrecisionOncologyNormalizer().fit(cancer_types, doids)

        cancer_types, doids = load_ncbi_cancer_goldstandard()
        cancer_types, doids = apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping)
        df_eval = normalizer.evaluate(cancer_types, doids, n_grams=3)

        assert f1_score(df_eval) > .8