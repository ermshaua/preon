# preon (PREcision Oncology Normalization)
preon is a fuzzy search tool for medical entities.

## Installation

You can install preon with PyPi:
`python -m pip install preon`

## Examples

Let's first import the normalizer and EBI drug names with CHEMBL ids.

```python3
>>> from preon.normalization import PrecisionOncologyNormalizer
>>> from preon.drug import store_ebi_drugs, load_ebi_drugs
```

Please download the <a href="https://www.ebi.ac.uk/chembl/g/#search_results/compounds">EBI compound CSV file</a> and store it as a local resource. This step only has to be performed when the resource file is created or updated. 

```python3
>>> store_ebi_drugs("/Users/Username/Downloads/compounds.csv")
```

Next, we can fit the normalizer with the drug names and ids as its reference data.

```python3
>>> drug_names, chembl_ids = load_ebi_drugs()
>>> normalizer = PrecisionOncologyNormalizer().fit(drug_names, chembl_ids)
```

We can now search for drug names and retrieve their CHEMBL ids. Let's search for the cancer drug "Avastin".

```python3
>>> normalizer.query("Avastin")
(['avastin'], [['CHEMBL1201583']], {'match_type': 'exact'})
```

As a result for our query, we get list of matching normalized drug names (in this case `['avastin']`), a list of associated CHEMBL ids for every returned drug name `[['CHEMBL1201583']]` and some meta information about the matching `{'match_type': 'exact'}`. We can also search for multi-token drug names like "Ixabepilone Epothilone B analog" and find CHEMBL ids for the relevant tokens.

```python3
>>> normalizer.query("Ixabepilone Epothilone B analog")
(['ixabepilone'], [['CHEMBL1201752']], {'match_type': 'substring'})
```

We find the relevant drug name `['ixabepilone']` and preon provides the meta information that the matching is based on a substring. On default, preon only looks for 1 matching token. It can also look for n-grams by setting the `n_grams` parameter in the query method. Let's take a harder example, say "Isavuconazonium", but misspell it as "Isavuconaconium".

```python3
>>> normalizer.query("Isavuconaconium")
(['isavuconazonium'], [['CHEMBL1183349']], {'match_type': 'partial', 'edit_distance': 0.067})
```

preon finds the correct drug "Isavuconazonium" and provides the meta information that it is a partial match with 7% distance. It returns drug names with a distance smaller than 20% on default. In order to change this parameter, set the `threshold` argument in the query method.

In a similar fashion you can normalize cancer types or genes. we provide gold standards for preon with which we test it. For more detail, see the example <a href="https://github.com/ermshaua/preon/tree/main/preon/examples">notebooks</a>. We also use preon in practice to normalize and integrate medical data in the <a href="https://predict.informatik.hu-berlin.de/">PREDICT</a> project.

## Citation

The preon package is actively maintained, updated and intended for application. If you use it in your scientific publication, we would appreciate the following <a href="https://doi.org/10.1101/2023.05.22.540912" target="_blank">citation</a>:

```
@article {preon2023,
	author = {Arik Ermshaus and Michael Piechotta and Gina R{\"u}ter and Ulrich Keilholz and Ulf Leser and Manuela Benary},
	title = {preon: Fast and accurate entity normalization for drug names and cancer types in precision oncology},
	year = {2023},
	doi = {10.1101/2023.05.22.540912},
	journal = {bioRxiv}
}
```
