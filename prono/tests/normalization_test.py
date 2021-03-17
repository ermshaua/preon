import unittest

from prono.normalization import PrecisionOncologyNormalizer
from prono.drug import load_ebi_drugs, load_charite_drug_goldstandard, load_database_drug_goldstandard, load_ctg_drug_goldstandard
from prono.tests.utils import f1_score


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