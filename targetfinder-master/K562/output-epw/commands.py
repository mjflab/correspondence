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
