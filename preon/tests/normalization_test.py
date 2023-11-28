import os
import shutil
import tempfile
import unittest

from preon.cancer import download_or_load_do_cancers, load_do_flat_mapping, apply_do_flat_mapping_to_ontology, \
    apply_do_flat_mapping_to_goldstandard, \
    load_database_cancer_goldstandard, load_ncbi_cancer_goldstandard
from preon.drug import load_ebi_drugs, load_charite_drug_goldstandard, load_database_drug_goldstandard, \
    load_ctg_drug_goldstandard, store_ebi_drugs, store_drugbank_drugs, load_drugbank_drugs
from preon.normalization import PrecisionOncologyNormalizer
from preon.tests.utils import f1_score


class DrugNormalizationTest(unittest.TestCase):

    def test_store_load_resources(self):
        abs_path = os.path.dirname(os.path.abspath(__file__))

        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = f"{abs_path}/../resources/ebi_drugs.csv"
            tmp_file = os.path.join(tmp_dir, "compounds.csv")
            shutil.copy(file_path, tmp_file)

            store_ebi_drugs(tmp_file)
            load_ebi_drugs()

            file_path = f"{abs_path}/../resources/drugbank_drugs.csv"
            tmp_file = os.path.join(tmp_dir, "compounds.csv")
            shutil.copy(file_path, tmp_file)

            store_drugbank_drugs(tmp_file)
            load_drugbank_drugs()

    def test_charite_drug_goldstandard(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer(enable_warnings=False).fit(drug_names, chembl_ids)

        drug_names, chembl_ids, _ = load_charite_drug_goldstandard()
        df_eval = normalizer.evaluate(drug_names, chembl_ids)

        assert f1_score(df_eval) > .9

    def test_databse_drug_goldstandard(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer(enable_warnings=False).fit(drug_names, chembl_ids)

        drug_names, chembl_ids, _ = load_database_drug_goldstandard()
        df_eval = normalizer.evaluate(drug_names, chembl_ids)

        assert f1_score(df_eval) > .9

    def test_ctg_drug_goldstandard(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer(enable_warnings=False).fit(drug_names, chembl_ids)

        drug_names, chembl_ids, _ = load_ctg_drug_goldstandard()
        df_eval = normalizer.evaluate(drug_names, chembl_ids)

        assert f1_score(df_eval) > .9

    def test_readme(self):
        drug_names, chembl_ids = load_ebi_drugs()
        normalizer = PrecisionOncologyNormalizer().fit(drug_names, chembl_ids)

        found_names, name_ids, meta_info = normalizer.query("Avastin")
        assert found_names == ['avastin']
        assert name_ids == [['CHEMBL1201583']]
        assert meta_info["match_type"] == "exact"

        found_names, name_ids, meta_info = normalizer.query("Ixabepilone Epothilone B analog")
        assert found_names == ['ixabepilone']
        assert name_ids == [['CHEMBL1201752']]
        assert meta_info["match_type"] == "substring"

        found_names, name_ids, meta_info = normalizer.query("Isavuconaconium")
        assert found_names == ['isavuconazonium']
        assert name_ids == [['CHEMBL1183349']]
        assert meta_info["match_type"] == "partial"
        assert meta_info["edit_distance"] == 0.067


class CancerNormalizationTest(unittest.TestCase):

    def test_database_cancer_goldstandard(self):
        cancer_types, doids = download_or_load_do_cancers()
        do_flat_mapping = load_do_flat_mapping()

        cancer_types, doids = apply_do_flat_mapping_to_ontology(cancer_types, doids, do_flat_mapping)
        normalizer = PrecisionOncologyNormalizer(enable_warnings=False).fit(cancer_types, doids)

        cancer_types, doids, _ = load_database_cancer_goldstandard()
        cancer_types, doids = apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping)
        df_eval = normalizer.evaluate(cancer_types, doids, n_grams=3)

        assert f1_score(df_eval) > .9

    def test_ncbi_cancer_goldstandard(self):
        cancer_types, doids = download_or_load_do_cancers()
        do_flat_mapping = load_do_flat_mapping()

        cancer_types, doids = apply_do_flat_mapping_to_ontology(cancer_types, doids, do_flat_mapping)
        normalizer = PrecisionOncologyNormalizer(enable_warnings=False).fit(cancer_types, doids)

        cancer_types, doids, _ = load_ncbi_cancer_goldstandard()
        cancer_types, doids = apply_do_flat_mapping_to_goldstandard(cancer_types, doids, do_flat_mapping)
        df_eval = normalizer.evaluate(cancer_types, doids, n_grams=3)

        assert f1_score(df_eval) > .7
