## Import Cell
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np

## Load Data Cell
expression_data = pd.read_csv('expression_data.csv')
gene_variances = expression_data.var(axis=1)
top_5000_genes = gene_variances.sort_values(ascending=False).head(5000).index
top_variable_expression = expression_data.loc[top_5000_genes].T

## Load Metadata Cell
target = pd.read_csv('../data/processed/group_annotation.csv', index_col=0)

target

## Data Splitting cell

from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

X_train, X_test, y_train, y_test = train_test_split(top_variable_expression, target['Group'], 
                                                    test_size=0.2, random_state=0, 
                                                    stratify=target['Group'])
## Logistic Regression: Cell 1

log_reg = Pipeline([
    ('scaler', StandardScaler()),
    ('log_reg', LogisticRegression(solver='liblinear'))
])

## Logistic Regression: Cell 2 
#hyperparameter tuning for logistic regression
param_grid = {
    'log_reg__C': [0.01, 0.1, 1, 10, 100, 1000],
    'log_reg__penalty': ['l1', 'l2']
}
grid_search_lr = GridSearchCV(log_reg, param_grid, cv=5, 
                              scoring='accuracy')
grid_search_lr.fit(X_train, y_train)
grid_search_lr.best_params_


## Logistic Regression: Cell 3

best_model_lr = grid_search_lr.best_estimator_ 
y_test_pred_lr = best_model_lr.predict(X_test)
print(classification_report(y_test, y_test_pred_lr))


# Save true and predicted labels with sample IDs to CSV for logistic regression group prediction
logreg_group_pred_df = pd.DataFrame({
    'Sample': X_test.index,
    'True_Label': y_test.values,
    'Pred_Label_LR': y_test_pred_lr
})

logreg_group_pred_df.to_csv('../results/logreg_group_predictions.csv', index=False)


## Logistic Regression: Cell 4 (Predict clusters)

from sklearn.cluster import AgglomerativeClustering
k = 3 

agglo = AgglomerativeClustering(n_clusters=k)
cluster_labels = agglo.fit_predict(top_variable_expression)


X_train, X_test, y_train_c, y_test_c = train_test_split(top_variable_expression, cluster_labels, 
                                                    test_size=0.2, random_state=0, 
                                                    stratify=cluster_labels)

## Logistic Regression: Cell 5 (Model training)
grid_search_lr.fit(X_train, y_train_c)
model_lr = grid_search_lr.best_estimator_
y_test_pred_cluster = model_lr.predict(X_test)
print(classification_report(y_test_c, y_test_pred_cluster))

# Save cluster label predictions, true cluster labels, and sample IDs to CSV
cluster_pred_df = pd.DataFrame({
    'Sample': X_test.index,
    'True_Cluster': y_test_c,
    'lr_pred_cluster': y_test_pred_cluster
})
cluster_pred_df.to_csv('../results/logreg_cluster_predictions.csv', index=False)

## Random Forest: Cell 1

#Random Forest Model

rf = Pipeline([
    ('scaler', StandardScaler()),
    ('rf', RandomForestClassifier(random_state=0))
])

## RF: Cell 2 

#random forest hyperparameter tuning 

param_grid_rf = {
    'rf__n_estimators': [50, 100, 200, 300],
    'rf__max_depth': [None, 3, 5, 10, 15, 20],
    'rf__min_samples_split': [2, 5, 10]
}

grid_search_rf = GridSearchCV(rf, param_grid_rf, cv=5,
                                scoring='accuracy')
grid_search_rf.fit(X_train, y_train)

## RF: Cell 3

best_model_rf = grid_search_rf.best_estimator_
y_test_pred_rf = best_model_rf.predict(X_test)
print(classification_report(y_test, y_test_pred_rf))


## Q4: retrain each predictive model using different numbers of genes. 

from sklearn.metrics import roc_auc_score

def training_adj_gene_count(model): 
    N_GENES = [10, 100, 1000, 10000]
    for n_genes in N_GENES:
        top_n_genes = gene_variances.sort_values(ascending=False).head(n_genes).index
        top_n_expression = expression_data.loc[top_n_genes].T
        X_train, X_test, y_train, y_test = train_test_split(top_n_expression, target['Group'], 
                                                            test_size=0.2, random_state=0, 
                                                        stratify=target['Group'])
        
        model.fit(X_train, y_train)
        y_test_pred = model.best_estimator_.predict(X_test)
        print(f"Classification report for top {n_genes} genes:")
        print(classification_report(y_test, y_test_pred))
        print("\n")

        y_test_score = model.best_estimator_.predict_proba(X_test)[:, 1]

        auc = roc_auc_score(y_test, y_test_score, multi_class='ovr')
        print(f"AUC for top {n_genes} genes: {auc}\n")

        



training_adj_gene_count(grid_search_lr)
training_adj_gene_count(grid_search_rf)










