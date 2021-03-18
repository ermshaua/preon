

def true_positives(df_eval):
    return df_eval[(df_eval['Match Type'] != "none") & (df_eval['Correct Match'] >= 1)]


def false_positives(df_eval):
    return df_eval[(df_eval['Match Type'] != "none") & (df_eval['Correct Match'] == 0)]


def true_negatives(df_eval):
    return df_eval[(df_eval['Match Type'] == "none") & (df_eval['Correct Match'] >= 1)]


def false_negatives(df_eval):
    return df_eval[(df_eval['Match Type'] == "none") & (df_eval['Correct Match'] == 0)]


def precision_score(df_eval):
    tp, fp = true_positives(df_eval).shape[0], false_positives(df_eval).shape[0]
    return tp / (tp + fp)


def recall_score(df_eval):
    tp, fn = true_positives(df_eval).shape[0], false_negatives(df_eval).shape[0]
    return tp / (tp + fn)


def f1_score(df_eval):
    precision, recall = precision_score(df_eval), recall_score(df_eval)
    return 2 * (precision * recall) / (precision + recall)