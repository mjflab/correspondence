# coding: utf-8
import pandas as pd
import numpy as np
import os
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.cross_validation import StratifiedKFold, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
windows_ripple_k562 = pd.read_csv("../../../ripple/data/hicfeatures/k562.windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
motif_files = os.walk("motifs").__next__()[2]


def getMotifClustersByCorrelation(X, num_motifs=79, thresh=0.2):
    link = linkage(X.T, 'single', 'correlation')
    cindexes = list(range(num_motifs))
    for l in link:
        c1 = int(l[0])
        c2 = int(l[1])
        distance = l[2]
        if distance < 0.2:
            cindexes.append([c1, c2])

    all_clusters = [[i,] for i in range(79)]
    for l in link:
        c1 = int(l[0])
        c2 = int(l[1])
        all_clusters.append(all_clusters[c1] + all_clusters[c2])

    cindexes = list(range(num_motifs))
    for i,l in enumerate(link):
        c1 = int(l[0])
        c2 = int(l[1])
        if l[2] < 0.2:
            cindexes.append(num_motifs+i)
            cindexes.remove(c1)
            cindexes.remove(c2)

    grouped_clusters = [all_clusters[i] for i in cindexes]
    motifs_to_use = []
    for g in grouped_clusters:
        motifs_to_use.append(X_train.columns[g[0]])
    return motifs_to_use, grouped_clusters



def load_data(region_file, motif_dir_prefix):
    windows = pd.read_csv(fn, names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
    for fn in motif_files:
        name = fn.split('.')[1]
        b = pd.read_csv(motif_dir_prefix + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
        b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
        windows = pd.merge(windows, b, on=["window_chrom", "window_start", "window_end"], how="left")
        predictors_df = windows.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
        labels = windows['label']
    return windows, predictors_df, labels

def evaludate_data(X, y):
    scoring = ['f1', 'roc_auc', 'precision', 'recall']
    estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 12)
    estimator = make_pipeline(StandardScaler(), estimator)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=12)
    cv = StratifiedKFold(y = y_train, n_folds = 8, shuffle = True, random_state = 12)
    scores = {}
    for s in scoring:
        scores[s] = cross_val_score(estimator, X_train, y_train, cv=cv, scoring=s, n_jobs=-1)
    estimator = estimator.fit(X_train, y_train)
    return (scores, (X_train, X_test, y_train, y_test), estimator, cv)

X_train_ripple_k562 = windows_ripple_k562.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
X_train_ripple_k562 = X_train_ripple_k562.fillna(0)
y_ripple_k562 = windows_ripple_k562['label']
cv_ripple_k562 = StratifiedKFold(y = y_ripple_k562, n_folds = 10, shuffle = True, random_state = 0)
cv_ripple_k562 = StratifiedKFold(y = y_ripple_k562, n_folds = 10, shuffle = True, random_state = 0)
scores_f1_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'f1', cv = cv_ripple_k562, n_jobs = -1)
estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
scores_f1_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'f1', cv = cv_ripple_k562, n_jobs = -1)
scores_f1_ripple_k562
precision_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'precision', cv = cv_ripple_k562, n_jobs = -1)
precision_ripple_k562
recall_ripple_k562 = cross_val_score(estimator, X_train_ripple_k562, y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
recall_ripple_k562
clf_rfc = RandomForestClassifier(n_estimators=4000, max_depth = 5, random_state = 0)
cross_val_score(clf_rfc, X_train_ripple_k562, y_ripple_k562, scoring = 'roc_auc', cv = cv_ripple_k562, n_jobs = -1)
ripple_k562_estimator = estimator.fit(X_train_ripple_k562, y_ripple_k562)
ripple_k562_estimator.feature_importances_
ripple_k562_importances = pd.Series(estimator.feature_importances_, index = predictors_df.columns).sort_values(ascending = False)
ripple_k562_importances = pd.Series(ripple_k562_estimator.feature_importances_, index = X_train_ripple_k562.columns).sort_values(ascending = False)
ripple_k562_importances.head(10)
ripple_k562_importances.head(20)
ripple_k562_importances.head(40)
ripple_k562_importances.head(50)
cross_val_score(estimator, X_train_ripple_k562["CTCF", "CTCF-ext"], y_ripple_k562, scoring = 'f1', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562[["CTCF", "CTCF-ext"]], y_ripple_k562, scoring = 'f1', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562[["CTCF", "CTCF-ext"]], y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562[["AP1", "RUNX1"]], y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562.ix[,[1,2,3]], y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562.ix[:,[1,2,3]], y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562.ix[:,[1,2,3]], y_ripple_k562, scoring = 'f1', cv = cv_ripple_k562, n_jobs = -1)
np.random.randn(10)
np.random.randn(10) * 79
np.random.rand(10) * 79
int(np.random.rand(10) * 79)
cross_val_score(estimator, X_train_ripple_k562.ix[,set(np.random.randint(79, size=10))], y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(estimator, X_train_ripple_k562.ix[:,set(np.random.randint(79, size=10))], y_ripple_k562, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
recall_ripple_k562
precision_ripple_k562
ripple_features = pd.read_csv("../../../ripple/data/hicfeatures/RaoHuntley_K562_allchr_5kb_enh_v1_examples.txt", sep="\t")
ripple_features
get_ipython().magic('pinfo pd.read_csv')
ripple_features = pd.read_csv("../../../ripple/data/hicfeatures/RaoHuntley_K562_allchr_5kb_enh_v1_examples.txt", sep="\t", index_col=1)
ripple_features
ripple_features = pd.read_csv("../../../ripple/data/hicfeatures/RaoHuntley_K562_allchr_5kb_enh_v1_examples.txt", sep="\t", index_col=0)
ripple_features.shape
ripple_features.index
y_ripple_features = ripple_features['class']
ripple_features.columns
y_ripple_features = ripple_features['Class']
X_ripple_features = ripple_features.drop(["Class",])
X_ripple_features = ripple_features.drop(["Class",], axis=1)
cross_val_score(estimator, X_ripple_features, y_ripple_features, scoring = 'roc_auc', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(clf_rfc, X_ripple_features, y_ripple_features, scoring = 'roc_auc', cv = cv_ripple_k562, n_jobs = -1)
windows = pd.read_csv("windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows = pd.merge(windows, b, on=["window_chrom", "window_start", "window_end"], how="left")

X_train = windows.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
y = windows['label']
X_train = X_train.fillna(0)
temp_predict = ripple_k562_estimator.predict(X_train)
f1_score(y, temp_predict)
sum(temp_predict == y)
y.shape
sum(temp_predict == y == 1)
precision_score(y, temp_predict)
recall_score(y, temp_predict)
roc_auc_score(y, temp_predict)
cross_val_score(clf_rfc, X_ripple_features, y_ripple_features, scoring = 'f1', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(clf_rfc, X_ripple_features, y_ripple_features, scoring = 'recall', cv = cv_ripple_k562, n_jobs = -1)
cross_val_score(clf_rfc, X_ripple_features, y_ripple_features, scoring = 'precision', cv = cv_ripple_k562, n_jobs = -1)
0.7 * 0.185 / (0.7+0.185) * 2
cv = StratifiedKFold(y = y, n_folds = 10, shuffle = True, random_state = 0)
scores_f1 cross_val_score(estimator, X_train, y, scoring = 'f1', cv = cv, n_jobs = -1)
scores_f1 = cross_val_score(estimator, X_train, y, scoring = 'f1', cv = cv, n_jobs = -1)
scores_f1
precision_k562 = cross_val_score(estimator, X_train, y, scoring = 'precision', cv = cv, n_jobs = -1)
precision_k562
precision_k562.mean()
recall_k562 = cross_val_score(estimator, X_train, y, scoring = 'recall', cv = cv, n_jobs = -1)
recall_k562
X_train.corr
X_train.corr()
plt
plt.imshow(X_train.corr())
plt.colorbar()
plt.imshow(X_train.corr(), interpolation="none")
plt.colorbar()
X_train.columns[28]
X_train.columns[19]
b = X_train.corr()
b
np.argmax(b[1])
np.argmax(b.ix[1,])
b
b.columns()
b.columns
np.argmax(b.ix[0,1:])
ordered_columns = b.columns[0]
ordered_columns
ordered_columns = [b.columns[0],]
for i in range(78):
    ordered_columns.append(np.argmax(b.ix[i, i+1:]))

ordered_columns
sorted(ordered_columns)
ordered_columns = [b.columns[0],]
for i in range(78):
    c = b.drop(ordered_columns)
    ordered_columns.append(np.argmax(c.ix[i,:]))

b
b.index
for i in range(78):
    c = b.drop(ordered_columns, axis=1)
    ordered_columns.append(np.argmax(c.ix[i,:]))

ordered_columns = [b.columns[0],]
for i in range(78):
    c = b.drop(ordered_columns, axis=1)
    ordered_columns.append(np.argmax(c.ix[i,:]))

ordered_columns
sorted(ordered_columns)
ordered_columns
b = X_train[:,ordered_columns].corr()
b = X_train.ix[:,ordered_columns].corr()
b
ordered_columns
ordered_columns = [b.columns[0],]
for i in range(78):
    c = b.drop(ordered_columns, axis=1)
    ordered_columns.append(np.argmax(c.ix[ordered_columns[-1],:]))

ordered_columns
b = X_train.ix[:,ordered_columns].corr()
b
plt.imshow(X_train.corr(), interpolation="none")
plt.colorbar()
b = X_train.corr()
ordered_columns = [b.columns[0],]
for i in range(78):
    c = b.drop(ordered_columns, axis=1)
    print( c.shape); ordered_columns.append(np.argmax(c.ix[ordered_columns[-1],:]))

ordered_columns
b = X_train.ix[:,ordered_columns].corr()
b
plt.imshow(b, interpolation="none")
plt.colorbar()
b
b[np.abs(b) < 0.5] = None
b
plt.imshow(b, interpolation="none")
plt.colorbar()
ordered_columns[:34]
ordered_columns[0,]
ordered_columns[0:1]
ordered_columns[0:1] + ordered_columns[35:36] +  ordered_columns[40:41] + ordered_columns[46:47] + ordered_columns[52:53] + ordered_columns[58:59] + ordered_columns[63:77] + ordered_columns[78:]
use_columns = ordered_columns[0:1] + ordered_columns[35:36] +  ordered_columns[40:41] + ordered_columns[46:47] + ordered_columns[52:53] + ordered_columns[58:59] + ordered_columns[63:77] + ordered_columns[78:]
cross_val_score(estimator, X[use_columns], y, cv=cv, n_jobs=-1, scoring="f1")
cross_val_score(estimator, X_train[use_columns], y, cv=cv, n_jobs=-1, scoring="f1")
scores_f1
b[np.abs(b) < 0.8] = None
plt.imshow(b, interpolation="none")
plt.colorbar()
ordered_columns
cross_val_score(estimator, X_train[ordered_columns[:10]], y, cv=cv, n_jobs=-1, scoring="f1")
cv = StratifiedKFold(y = y, n_folds = 10, shuffle = True, random_state = 123)
cross_val_score(estimator, X_train[ordered_columns[:10]], y, cv=cv, n_jobs=-1, scoring="f1")
from sklearn.feature_selection import RFE
selector = RFE(estimator, 10, step=0.1)
selector.fit(X_train, y)
selector = _
selector.support_
X.columns[selector.support_]
X_train.columns[selector.support_]
cross_val_score(estimator, X_train[selector.support_], y, cv=cv, n_jobs=-1, scoring="f1")
cross_val_score(estimator, X_train.ix[:,selector.support_], y, cv=cv, n_jobs=-1, scoring="f1")
scores_f1
selector = RFE(estimator, 20, step=1)
selector = selector.fit(X_train, y)
selector.support_
X_train.columns[selector.support_]
ordered_columns
random_index = np.random.randint(10)
random_index
random_index = np.random.randint(10, size = 10)
random_index
random_index = np.random.randint(79, size = 10)
cross_val_score(estimator, X_train.ix[:,random_index], y, cv=cv, n_jobs=-1, scoring="f1")
X_train.columns[random_index]
len(use_columns)
use_columns
cross_val_score(estimator, X_train_ripple_k562.ix[:,random_index], y_ripple_k562, cv=cv_ripple_k562, n_jobs=-1, scoring="f1")
cross_val_score(estimator, X_train_ripple_k562.ix[:,use_columns], y_ripple_k562, cv=cv_ripple_k562, n_jobs=-1, scoring="f1")
scores_f1_ripple_k562
pdist
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
distance =pdist(X_train, 'correlation')
distance =pdist(X_train.T, 'correlation')
distance
len(distance)
79*78/2
dendrogram(linkage(X.T, 'single', 'correlation'), color_threshold=0)
dendrogram(linkage(X_train.T, 'single', 'correlation'), color_threshold=0)
d = _
dendrogram(linkage(X_train.T, 'single', 'correlation'), color_threshold=0)
dendrogram(linkage(X_train.T, 'single', 'correlation'), color_threshold=0, labels=X_train.columns)
dendrogram(linkage(X_train.T, 'single', 'correlation'), color_threshold=0, labels=X_train.columns)
dendrogram(linkage(X_train.T, 'single', 'correlation'), color_threshold=0, labels=X_train.columns, leaf_font_size=12)
dendrogram(linkage(X_train.T, 'single', 'correlation'), color_threshold=0, labels=X_train.columns, leaf_font_size=12)
d['icoord'].shape
len(d['icoord'])
d['icoord'][-1]
d['icoord'][-2]
d['icoord'][-3]
sorted(d['icoord'])
(775 - 5)/5
(775 - 5)/79
(775 - 5)/78
linkage(X_train.T, 'single', 'correlation')
linkage(X_train.T, 'single', 'correlation').shape
print(linkage(X_train.T, 'single', 'correlation'))
X_train.columns[-3:]
X_train.columns[66]
X_train.columns[72]

motifs_to_use
cross_val_score(estimator, X_train_ripple_k562.ix[:,motifs_to_use], y_ripple_k562, cv=cv_ripple_k562, n_jobs=-1, scoring="f1")
cross_val_score(estimator, X_train.ix[:,motifs_to_use], y, cv=cv, n_jobs=-1, scoring="f1")
a = _
a.mean()
scores_f1.mean()
scores_f1
cross_val_score(estimator, X_train.ix[:,motifs_to_use], y, cv=cv, n_jobs=-1, scoring="precision")
cross_val_score(estimator, X_train.ix[:,motifs_to_use], y, cv=cv, n_jobs=-1, scoring="recall")
recall_k562
windows_gm12878 = pd.read_csv("../../GM12878/output-epw/windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
for fn in motif_files:
    name = fn.split('.')[1]
    b = pd.read_csv("../../GM12878/output-epw/motifs/" + fn, sep="\t", header=None, names=["window_chrom", "window_start", "window_end", name])
    b[name] = b[name] * 1.0 / (b["window_end"] - b["window_start"])
    windows_gm12878 = pd.merge(windows_gm12878, b, on=["window_chrom", "window_start", "window_end"], how="left")

X_train_gm12878 = windows_gm12878.drop(["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter" ], axis=1)
X_train_gm12878 = X_train_gm12878.fillna(0)
y_gm12878 = windows_gm12878['label']
scores_f1_selected = cross_val_score(estimator, X_train.ix[:,motifs_to_use], y, cv=cv, n_jobs=-1, scoring="f1")
cv_gm12878 = StratifiedKFold(y = y_gm12878, n_folds = 10, shuffle = True, random_state = 123)
scores_f1_gm12878_selected = cross_val_score(estimator, X_train_gm12878.ix[:,motifs_to_use], y_gm12878, cv=cv_gm12878, n_jobs=-1, scoring="f1")
scores_f1_gm12878_selected
metrics_names = ['f1', 'roc_auc', 'precision', 'recall']
scores_gm12878_selected = []
for m in metrics_names:
    scores_gm12878_selected.append( cross_val_score(estimator, X_train_gm12878.ix[:,motifs_to_use], y_gm12878, cv=cv_gm12878, n_jobs=-1, scoring=m))

scores_gm12878_selected
link_gm12878 = linkage(X_train_gm12878.T, 'single', 'correlation')
cindexes = list(range(79))
for i,l in enumerate(link_gm12878):
    c1 = int(l[0])
    c2 = int(l[1])
    if l[2] < 0.2:
        cindexes.append(79+i)
        cindexes.remove(c1)
        cindexes.remove(c2)

all_clusters = [[i,] for i in range(79)]
for l in link_gm12878:
    c1 = int(l[0])
    c2 = int(l[1])
    all_clusters.append(all_clusters[c1] + all_clusters[c2])

grouped_clusters_gm12878 = [all_clusters[i] for i in cindexes]
grouped_clusters
print( grouped_clusters)
print( grouped_clusters_gm12878)
len(grouped_clusters)
len(grouped_clusters_gm12878)
X_train.mean()
X_train.mean() - X_train_gm12878.mean()
print(X_train.mean() - X_train_gm12878.mean())
print(sorted(X_train.mean() - X_train_gm12878.mean()))
print(sorted(X_train.mean() - X_train_gm12878.mean()) / X_train.mean())
print(sorted((X_train.mean() - X_train_gm12878.mean()) / X_train.mean()))
print((X_train.mean() - X_train_gm12878.mean()) / X_train.mean())
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
def evaludate_data(X, y):
    scoring = ['f1', 'roc_auc', 'precision', 'recall']
    estimator = GradientBoostingClassifier(n_estimators = 400, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 12)
    estimator = make_pipeline(StandardScaler(), estimator)
    X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, y, test_size=0.2, random_state=12)
    cv = StratifiedKFold(y = y_train, n_folds = 8, shuffle = True, random_state = 12)
    scores = {}
    for s in scoring:
        scores[s] = cross_val_score(estimator, X_train, y_train, cv=cv, scoring=s)
    estimator = estimator.fit(X_train, y_train)
    return (scores, (X_train, X_test, y_train, y_test), estimator, cv)
from sklearn.cross_validation import train_test_split
tf_k562 = evaludate_data(X_train, y_train)
tf_k562 = evaludate_data(X_train, y)
def evaludate_data(X, y):
    scoring = ['f1', 'roc_auc', 'precision', 'recall']
    estimator = GradientBoostingClassifier(n_estimators = 400, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 12)
    estimator = make_pipeline(StandardScaler(), estimator)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=12)
    cv = StratifiedKFold(y = y_train, n_folds = 8, shuffle = True, random_state = 12)
    scores = {}
    for s in scoring:
        scores[s] = cross_val_score(estimator, X_train, y_train, cv=cv, scoring=s)
    estimator = estimator.fit(X_train, y_train)
    return (scores, (X_train, X_test, y_train, y_test), estimator, cv)
tf_k562 = evaludate_data(X_train, y)
tf_k562
len(tf_k562)
tf_k562[2]
f1_score(y_gm12878, tf_k562[2].predict(X_train_gm12878))
f1_score(tf_k562[1][3], tf_k562[2].predict(tf_k562[1][1]))
roc_auc_score_score(tf_k562[1][3], tf_k562[2].predict(tf_k562[1][1]))
roc_auc_score(tf_k562[1][3], tf_k562[2].predict(tf_k562[1][1]))
precision_score(tf_k562[1][3], tf_k562[2].predict(tf_k562[1][1]))
recall_score(tf_k562[1][3], tf_k562[2].predict(tf_k562[1][1]))
tf_k562[0]
precision_score(y_gm12878, tf_k562[2].predict(X_train_gm12878))
recall_score(y_gm12878, tf_k562[2].predict(X_train_gm12878))
recall_score(y_gm12878, tf_k562[2].predict(X_train_gm12878))
tf_k562_selected = evaludate_data(X_train.ix[:,motif_to_use], y)
tf_k562_selected = evaludate_data(X_train.ix[:,motifs_to_use], y)
tf_k562_selected[0]
recall_score(y_gm12878, tf_k562_selected[2].predict(X_train_gm12878.ix[:,motifs_to_use]))
tf_gm12878_selected = evaludate_data(X_train_gm12878.ix[:,motifs_to_use], y_gm12878)
tf_gm12878_selected[0]
def evaludate_data(X, y):
    scoring = ['f1', 'roc_auc', 'precision', 'recall']
    estimator = GradientBoostingClassifier(n_estimators = 400, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 12)
    #estimator = make_pipeline(StandardScaler(), estimator)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=12)
    cv = StratifiedKFold(y = y_train, n_folds = 8, shuffle = True, random_state = 12)
    scores = {}
    for s in scoring:
        scores[s] = cross_val_score(estimator, X_train, y_train, cv=cv, scoring=s)
    estimator = estimator.fit(X_train, y_train)
    return (scores, (X_train, X_test, y_train, y_test), estimator, cv)
tf_gm12878_selected = evaludate_data(X_train_gm12878.ix[:,motifs_to_use], y_gm12878)
tf_gm12878_selected[0]
def evaludate_data(X, y):
    scoring = ['f1', 'roc_auc', 'precision', 'recall']
    estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 12)
    estimator = make_pipeline(StandardScaler(), estimator)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=12)
    cv = StratifiedKFold(y = y_train, n_folds = 8, shuffle = True, random_state = 12)
    scores = {}
    for s in scoring:
        scores[s] = cross_val_score(estimator, X_train, y_train, cv=cv, scoring=s)
    estimator = estimator.fit(X_train, y_train)
    return (scores, (X_train, X_test, y_train, y_test), estimator, cv)
def evaludate_data(X, y):
    scoring = ['f1', 'roc_auc', 'precision', 'recall']
    estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 12)
    estimator = make_pipeline(StandardScaler(), estimator)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=12)
    cv = StratifiedKFold(y = y_train, n_folds = 8, shuffle = True, random_state = 12)
    scores = {}
    for s in scoring:
        scores[s] = cross_val_score(estimator, X_train, y_train, cv=cv, scoring=s, n_jobs=-1)
    estimator = estimator.fit(X_train, y_train)
    return (scores, (X_train, X_test, y_train, y_test), estimator, cv)
tf_gm12878_selected = evaludate_data(X_train_gm12878.ix[:,motifs_to_use], y_gm12878)
tf_gm12878_selected[0]
tf_k562_selected = evaludate_data(X_train.ix[:,motifs_to_use], y)
tf_k562_selected[0]
tf_k562 = evaludate_data(X_train, y)
tf_k562[0]
recall_score(y_gm12878, tf_k562_selected[2].predict(X_train_gm12878.ix[:,motifs_to_use]))
recall_score(tf_k562_selected[1][3], tf_k562_selected[2].predict(tf_k562_selected[1][2].ix[:,motifs_to_use]))
recall_score(tf_k562_selected[1][3], tf_k562_selected[2].predict(tf_k562_selected[1][2]))
recall_score(tf_k562_selected[1][3], tf_k562_selected[2].predict(tf_k562_selected[1][2]))
tf_k562_selected[1][2].shape
recall_score(tf_k562_selected[1][3], tf_k562_selected[2].predict(tf_k562_selected[1][1]))
f1_score(tf_k562_selected[1][3], tf_k562_selected[2].predict(tf_k562_selected[1][1]))
precision_score(tf_k562_selected[1][3], tf_k562_selected[2].predict(tf_k562_selected[1][1]))
nonpredictors = ['enhancer_chrom', 'enhancer_start', 'enhancer_end', 'promoter_chrom', 'promoter_start', 'promoter_end', 'window_chrom', 'window_start', 'window_end', 'window_name', 'active_promoters_in_window', 'interactions_in_window', 'enhancer_distance_to_promoter', 'bin', 'label']
training_df = pd.read_hdf('training.h5', 'training').set_index(['enhancer_name', 'promoter_name'])
predictors_df = training_df.drop(nonpredictors, axis = 1)
labels = training_df['label']
estimator = GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
cv = StratifiedKFold(y = labels, n_folds = 10, shuffle = True, random_state = 0)
scores = cross_val_score(estimator, predictors_df, labels, scoring = 'f1', cv = cv, n_jobs = -1)
score
scores
tf_epw_k562 = evaludate_data(predictors_df, labels)
tf_epw_k562[0]
predictors_df.columns
tf_epw_k562_window = evaludate_data(predictors_df.ix[:,range(2,predictors_df.shape[1], 3)], labels)
tf_epw_k562_window
tf_epw_k562_window[0]
tf_epw_k562[0]
training_df_gm12878 = pd.read_hdf('../../../targetfinder-master/paper/targetfinder/GM12878/output-epw/training.h5', 'training').set_index(['enhancer_name', 'promoter_name'])
predictors_df_gm12878 = training_df_gm12878.drop(nonpredictors, axis = 1)
labels_gm12878 = training_df_gm12878['label']
tf_epw_gm12878 = evaludate_data(predictors_df_gm12878, labels_gm12878)
tf_epw_gm12878[0]
tf_epw_gm12878 = evaludate_data(predictors_df_gm12878.ix[:,range(2,predictors_df_gm12878.shape[1], 3)], labels_gm12878)
tf_epw_gm12878_window = tf_epw_gm12878
tf_epw_gm12878_window[0]
tf_gm12878_selected[0]
tf_k562_selected[2].feature_importances_
tf_k562_selected[2].feature_importances_
dir(tf_k562_selected[2])
dir(tf_k562_selected[2].step[1])
dir(tf_k562_selected[2].steps[1])
tf_k562_selected[2].steps
tf_k562_selected[2].steps[1][1].feature_importances_
a = evaludate_data(X_train.ix[:,23], y)
a = evaludate_data(X_train_gm12878.ix[:,23:24], y_gm12878)
a[0]
windows = pd.read_csv("windows.bed", names=["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"], sep="\t", header=None)
a = evaludate_data(windows.ix[:,['enhancer_distance_to_promoter']], windows['label'])
a[0]
a[1][3]
a[1][3].shape
a[1][2].shape
sum(a[2].predict(a[1][1]) == a[1][3])
sum(a[1][3] == 1)
sum(a[2].predict(a[1][3][a[1][1]) == a[1][3]])
sum(a[1][3][a[2].predict([a[1][1]) == a[1][3]]])
sum(a[1][3][a[2].predict(a[1][1]) == a[1][3]]])
sum(a[1][3][a[2].predict(a[1][1]) == a[1][3]])
sum(a[2].predict(a[1][1]))
tf_k562[2].steps[1].oob_improvement_
tf_k562[2].steps[1][1].oob_improvement_
tf_k562[2].steps[1][1]['oob_improvement_']
tf_k562[2].steps[1][1]['oob_improvement_']
tf_k562[2].steps[1][1]
dir(tf_k562[2].steps[1][1])
tf_k562[2].steps[1][1].decision_function
tf_k562[2].steps[1][1].decision_function()
import scikit
import sklearn
sklearn.__version__
tf_k562[0]
scores_f1.mean()
evaludate_data(X_train, y)
evaludate_data(X_train, y)
a = evaludate_data(X_train, y)
a = evaludate_data(X_train, y)
