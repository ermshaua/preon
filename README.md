# pronto (PRecision Oncology NormalizaTiOn)
pronto is a fuzzy search tool for medical entities.

## Installation

You can install pronto with PyPi:
`python -m pip install {todo}`

## Examples

Let's first import the normalizer and EBI drug names with CHEMBL ids.

```python3
>>> from pronto.normalization import PrecisionOncologyNormalizer
>>> from pronto.drug import load_ebi_drugs
```

Let's fit the normalizer with the drug names and ids as its reference data.

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

We find the relevant drug name `['ixabepilone']` and pronto provides the meta information that the matching is based on a substring. On default, pronto only looks for 1 matching token. It can also look for n-grams by setting the `n_grams` parameter in the query method. Let's take a harder example, say "Isavuconazonium", but misspell it as "Isavuconaconium".

```python3
>>> normalizer.query("Isavuconaconium")
(['isavuconazonium'], [['CHEMBL1183349']], {'match_type': 'partial', 'edit_distance': 0.0667})
```

pronto finds the correct drug "Isavuconazonium" and provides the meta information that it is a partial match with 7% distance. It returns drug names with a distance smaller than 20% on default. In order to change this parameter, set the `threshold` argument in the query method.
