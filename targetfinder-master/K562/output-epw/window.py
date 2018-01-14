# coding: utf-8
f = open("windows.bed")
f.close()
import pandas as pd
a = pd.read_csv("windows.bed")
a
a = pd.read_csv("windows.bed", sep="\t")
a
a = pd.read_csv("windows.bed", sep="\t", header=Fals)
a = pd.read_csv("windows.bed", sep="\t", header=False)
a = pd.read_csv("windows.bed", sep="\t", header=None)
a
a.columns
a
b = pd.read_csv("motifs/window.AP1.bed", header=None, sep="\t", usecols=[3,])
b
b.columns.set_names(["AP1",])
b
b.columns = ["AP1"]
b
pd.concat([a,b], axis=1)
b.shape
a.shape
b = pd.read_csv("motifs/window.AP1.bed", header=None, sep="\t")
b
a
pd.merge(a,b, left_on=[0,1,2], right_on=[0,1,2])
pd.merge(a,b, left_on=[0,1,2], right_on=[0,1,2], how="left")
b
a.columns
a.columns = ["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"]
b.columns= [ "window_chrom", "window_start", "window_end", "AP1"]
pd.merge(a,b, on=["window_chrom", "window_start", "window_end"], how="left")
import os
os.walk(".")
a = _
a
a.next()
a[0].next()
a.next
dir(a)
a.close()
os.walk(".")
a = _
a[0]
b, _, _ = os.walk('.')
b, _ = os.walk('.')
b
os.walk(".")
os.walk(".").next()
dir(os.walk("."))
os.walk(".").send
os.walk(".").send()
os.walk(".").__net__
os.walk(".").__next__
os.walk(".").__next__()
for a ,b, c in os.walk("motifs"):
    print a,b,c
    
for a ,b, c in os.walk("motifs"):
    print( a,b,c)
    
windows = pd.read_csv("windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
windows
motif_files = os.walk("motifs").__next__()[2]
motif_files
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows = pd.merge(windows, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
windows
X_train = windows.filter(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ])
X_train.shape
windows.shape
X_train = windows.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
X_train.shape
X_train.index = windows.window_name
X_train
y = windows['label']
from sklearn.cross_validation import StratifiedKFold, cross_val_score
from sklearn.ensemble import GradientBoostingClassifier
estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
cv = StratifiedKFold(y = labels, n_folds = 10, shuffle = True, random_state = 0)
cv = StratifiedKFold(y = y, n_folds = 10, shuffle = True, random_state = 0)
scores = cross_val_score(estimator, predictors_df, labels, scoring = 'f1', cv = cv, n_jobs = -1)
scores = cross_val_score(estimator, X_train, y, scoring = 'f1', cv = cv, n_jobs = -1)
X_train.fillna(0)
X_train = X_train.fillna(0)
scores = cross_val_score(estimator, X_train, y, scoring = 'f1', cv = cv, n_jobs = -1)
scores
scores.mean()
precisions = cross_val_score(estimator, X_train, y, scoring = 'precision', cv = cv, n_jobs = -1)
precisions
recalls = cross_val_score(estimator, X_train, y, scoring = 'recall', cv = cv, n_jobs = -1)
recalls
rocs = cross_val_score(estimator, X_train, y, scoring = 'roc_auc', cv = cv, n_jobs = -1)
rocs
estimator.fit(X_train, y)
importances = pd.Series(estimator.feature_importances_, index = X_train.columns).sort_values(ascending = False)
importances.head(10)
importances.head(16)
importances.head(20)
importances.head(40)
from sklearn.svm import SVC
clf_svm  = SVC(kernal="linear")
clf_svm  = SVC(kernel="linear")
cross_val_score(clf_svm, X_train, y, scoring="f1", cv=cv, njobs=-1)
cross_val_score(clf_svm, X_train, y, scoring="f1", cv=cv, n_jobs=-1)
clf_svm  = SVC()
scorss
scores
precsions
precisions
cross_val_score(clf_svm, X_train, y, scoring="f1", cv=cv, n_jobs=-1)
X_train
X_train[y==1].mean()
X_train[y==1]
X_train.ix[y==1,]
y
X_train
y.index = X_train.index
X_train[y==1]
X_train[y==1].mean()
X_train[y==1].mean()["CTCF"]
X_train[y==0].mean()["CTCF"]
get_ipython().magic('save "commands.py" 1-117')
windows_gm12878 = pd.read_csv("../../GM12878/output-epw/windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows_gm12878 = pd.merge(windows, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
windows_gm12878.shape
X_train_gm12878 = windows_gm12878.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
y_gm12878 = windows['label']
y_gm12878 = windows_gm12878['label']
X_train_gm12878.fillna(0)
X_train_gm12878 = X_train_gm12878.fillna(0)
k562_pre_gm12878 = estimator.predict(X_train_gm12878)
X_train_gm12878.shape
X_train_gm12878.columns
windows_gm12878 = pd.read_csv("../../GM12878/output-epw/windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows_gm12878 = pd.merge(windows, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
len(motif_files)
windows_gm12878.columns
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    if name == "PRDM1": print(name) ; windows_gm12878 = pd.merge(windows, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
windows_gm12878 = pd.read_csv("../../GM12878/output-epw/windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
windows_gm12878.columns
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    if name == "PRDM1": print(name) ; windows_gm12878 = pd.merge(windows_gm12878, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
X_train_gm12878 = windows_gm12878.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
X_train_gm12878.fillna(0)
X_train_gm12878 = X_train_gm12878.fillna(0)
k562_pre_gm12878 = estimator.predict(X_train_gm12878)
X_train_gm12878.shape
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows_gm12878 = pd.merge(windows_gm12878, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
windows_gm12878 = pd.read_csv("../../GM12878/output-epw/windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows_gm12878 = pd.merge(windows_gm12878, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
X_train_gm12878 = windows_gm12878.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
X_train_gm12878 = X_train_gm12878.fillna(0)
X_train_gm12878.shape
k562_pre_gm12878 = estimator.predict(X_train_gm12878)
from sklearn.metrics import f1_score, precision_score
f1_score(y_gm12878, k562_pre_gm12878)
y_gm12878 = windows_gm12878['label']
f1_score(y_gm12878, k562_pre_gm12878)
precision_score(y_gm12878, k562_pre_gm12878)
from sklearn.metrics import recall_score
recall_score(y_gm12878, k562_pre_gm12878)
X_combined = pd.concat([X_train, X_train_gm12878], ignore_index=True)
X_combined.shape
y_combined = pd.concat([y, y_gm12878], ignore_index=True)
y_combined.shape
cv_combined = StratifiedKFold(y = y_combined, n_folds = 10, shuffle = True, random_state = 0)
cross_val_score(estimator, X_combined, y_combined, scoring="f1", cv=cv_combined, n_jobs=-1)
scores
combined_estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
combined_estimator.fit(X_combined, y_combined)
a = _
a
combined_estimator.feature_importances_
estimator.feature_importances_
combined_importances = pd.Series(combined_estimator.feature_importances_, index = X_combined.columns).sort_values(ascending = False)
combined_importances.head(20)
cross_val_score(clf_svm, X_combined, y_combined, cv = cv_combined, n_jobs=-1)
cross_val_score(clf_svm, X_combined, y_combined, cv = cv_combined, scoring="f1", n_jobs=-1)
cross_val_score(clf_svm, X_combined, y_combined, cv = cv_combined, scoring="precision", n_jobs=-1)
from sklearn.svm import LinearSVC
clf_svm = LinearSVC()
cross_val_score(clf_svm, X_combined, y_combined, cv = cv_combined, scoring="precision", n_jobs=-1)
clf_svm = LinearSVC(C=10)
cross_val_score(clf_svm, X_combined, y_combined, cv = cv_combined, scoring="precision", n_jobs=-1)
clf_svm = LinearSVC(C=100)
cross_val_score(clf_svm, X_combined, y_combined, cv = cv_combined, scoring="precision", n_jobs=-1)
combined_precision = cross_val_score(estimator, X_combined, y_combined, scoring="precision", cv=cv_combined, n_jobs=-1)
combined_precision
combined_recall = cross_val_score(estimator, X_combined, y_combined, scoring="recall", cv=cv_combined, n_jobs=-1)
combined_recall
combined_recall.mean()
combined_precision.mean()
X_train.index[0]
import numpy as np
np.random.randn(X_train.shape)
X_random = np.random.randn(*(X_train.shape))
X_random.shape
X_random = pd.DataFrame(np.random.randn(*(X_train.shape)))
X_random.shape
X_random.mean()
a = np.random.shuffle(range(y_combined.shape[0]))
a = np.random.shuffle(range(y_combined.shape[0]))
idx = range(y_combined.shape[0])
np.random.shuffle(idx)
type(idx)
idx = range(y_combined.shape[0])
dir(idx)
idx[0]
idx = list(range(y_combined.shape[0]))
idx
np.random.shuffle(idx)
idx
y_combined_shuffled = y_combined[idx]
y_combined_shuffled - y_combined
y_combined_shuffled = y_combined[idx]
y_combined_shuffled
y_combined_shuffled.index = y_combined.index
y_combined_shuffled - y_combined
cross_val_score(combined_estimator, X_combined, y_combined_shuffled, cv = cv_combined, scoring="precision", n_jobs=-1)
windows_ripple_k562 = pd.read_csv("../../../ripple/data/hicfeatures/k562.windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
windows_ripple_k562.shape
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../../ripple/data/hicfeatures/motifs/k562_" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows_ripple_k562 = pd.merge(windows_ripple_k562, b, on=["window_chrom", "window_start", "window_end"], how="left")
    
windows_ripple_k562.shape
X_train_ripple_k562 = windows_ripple_k562.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
y_ripple_k562 = windows_ripple_k562['label']
X_train_ripple_k562.fillna(0)
X_train_ripple_k562 = X_train_ripple_k562.fillna(0)
scores_f1_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'f1', cv = cv, n_jobs = -1)
cv_ripple_k562 = StratifiedKFold(y = y_ripple_k562, n_folds = 10, shuffle = True, random_state = 0)
scores_f1_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'f1', cv = cv, n_jobs = -1)
scores_f1_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'f1', cv = cv, n_jobs = -1)
if __name__== "__main__":
    scores_f1_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'f1', cv = cv, n_jobs = -1)
    
scores
scores.mean()
precisions.mean()
recall.mean()
recalls.mean()
scores.std
scores.std()
precisions.std()
recalls.std()
combined_precision
combined_precision.mean()
combined_precision.std()
combined_recall.mean()
combined_recall.std()
fdsf
