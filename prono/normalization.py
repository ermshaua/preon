import numpy as np
import pandas as pd

import time
import jellyfish
import re

from tqdm import tqdm
from nltk import ngrams


class PrecisionOncologyNormalizer:

    def _transform_name(self, name):
        name = name.lower()
        name = re.sub("[^a-zA-Z0-9]", '', name)
        return name

    def fit(self, names, name_ids):
        self.names = dict()

        for name, name_id in zip(names, name_ids):
            name = self._transform_name(name)

            if name not in self.names: self.names[name] = list()
            self.names[name].append(name_id)

        return self

    def _get_query_result(self, found_names, meta_info):
        found_names = np.unique(found_names).tolist()
        name_ids = [np.unique(self.names[found_name]).tolist() for found_name in found_names]
        return found_names, name_ids, meta_info

    def query(self, query_name, match_type="all", threshold=.2, n_grams=1):
        _query_name = self._transform_name(query_name)
        meta_info = dict()

        if match_type in ("exact", "all"):
            # try to find the trivial match
            meta_info["match_type"] = "exact"

            if _query_name in self.names:
                return self._get_query_result([_query_name], meta_info)

        if match_type in ("substring", "all"):
            # try to find trivial substring match
            meta_info["match_type"] = "substring"

            substrings = query_name.split(" ")

            for i_gram in range(2, n_grams+1):
                grams = ngrams(substrings, i_gram)
                grams = [" ".join(gram) for gram in grams]
                substrings.extend(grams)

            matches = []

            for token in substrings:
                _token = self._transform_name(token)

                if len(_token) == 0: continue
                if _token in self.names: matches.append(_token)

            if len(matches) > 0:
                return self._get_query_result(matches, meta_info)

        if match_type in ("partial", "all"):
            # otherwise, try partial matching
            meta_info["match_type"] = "partial"
            names = list(self.names)

            distances = []

            # calculate distances for query string
            for name in names:
                dist = jellyfish.levenshtein_distance(_query_name, name) / max(len(_query_name), len(name))
                distances.append(dist)

            distances = np.array(distances)

            if np.min(distances) > threshold:
                return None

            meta_info["edit_distance"] = np.min(distances)
            names_idx = np.isin(distances, meta_info["edit_distance"])
            names = np.array(names)

            return self._get_query_result(names[names_idx].tolist(), meta_info)

        return None

    def transform(self, names, verbose=0, **query_args):
        df = []

        for name in tqdm(names, disable=verbose < 1):
            query_time = time.process_time()
            res = self.query(name, **query_args)
            query_time = time.process_time() - query_time

            if res is None:
                df.append((name, [], [[None]], "none", None, query_time))
                continue

            found_names, found_ids, meta_info = res
            df.append((name, found_names, found_ids, meta_info["match_type"], meta_info.get("edit_distance", None), query_time))

        df = pd.DataFrame.from_records(df, columns=["Name", "Found Names", "Found Name IDs", "Match Type", "Edit Distance", "Query Time"])
        return df

    def evaluate(self, X, y, verbose=0, **query_args):
        df = self.transform(X, verbose=verbose, **query_args)
        df["Name IDs"] = y

        # helper procedure
        flatten = lambda data: [item for sub in data for item in sub]

        correct_matches = []
        for name_ids, found_ids in zip(y, df["Found Name IDs"].tolist()):
            matches = 0

            for found_id in flatten(found_ids):
                if found_id in name_ids:
                    matches += 1

            correct_matches.append(matches)

        df["Correct Match"] = correct_matches

        df = df[["Name", "Found Names", "Name IDs", "Found Name IDs", "Match Type", "Correct Match", "Edit Distance",
                 "Query Time"]]
        return df