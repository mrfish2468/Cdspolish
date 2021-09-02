import pandas as pd
import sys
#from libsvm.svmutil import *
from modules import preprocessing as p
from sklearn.metrics import confusion_matrix,classification_report
from sklearn import preprocessing
import time
from multiprocessing import Pool
import joblib
from joblib import Parallel, delayed
import numpy as np
from functools import partial
from contextlib import contextmanager
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn.metrics import confusion_matrix,classification_report, accuracy_score


from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, roc_curve, roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import confusion_matrix,classification_report, accuracy_score
from sklearn.model_selection import cross_val_score, GridSearchCV, train_test_split


def predict_class(input_data, model_FILE):
	model = joblib.load(model_FILE).set_params(n_jobs=1)
	p_label=model.predict(input_data.values.tolist())
	return p_label

def predict_proba(input_data, model_FILE):
	model = joblib.load(model_FILE).set_params(n_jobs=1)
	p_proba=model.predict_proba(input_data.values.tolist())
	return p_proba


def predict(df_feather, model_file, path):
	X = pd.read_feather(df_feather)
	Y = X['label'].values
	del X['label']
	position = X['position']
	X = p.double_ins_chaos(X)
	X = p.preprocessing(X)
	X = pd.DataFrame(X.drop(['position'], axis=1))

	print(X.columns)
	print('training shape: ', X.shape)
	size = 100000
	list_of_X = [X.loc[i:i+size-1,:] for i in range(0, len(X),size)]


	#get label
	pool = Pool(32)
	results = pool.map(partial(predict_class, model_FILE=model_file), list_of_X)
	pool.close()
	pool.join()

	result = []
	for i in results:
		result.extend(i)
	print(len(result))



	#get proba
	pool = Pool(32)
	proba = pool.map(partial(predict_class, model_FILE=model_file), list_of_X)
	pool.close()
	pool.join()

	proba_final = []
	for i in proba:
		proba_final.extend(i)
	print(len(proba_final))

	print("cm", confusion_matrix(Y.tolist(),result, labels = [0,1,2,3,4,5]), "\n")

	#confusion matrix
	confusion = path + '/confusion.png'
	cm =confusion_matrix(Y.tolist(),result, labels = [0,1,2,3,4,5])
	plt.figure(figsize=(15,15))
	ax= plt.subplot()
	sns.heatmap(cm, annot=True, linewidths=.5, square = True, cmap = 'Blues_r', fmt='g', annot_kws={"size":15}, ax=ax);
	ax.set_xticklabels( ['A','T','C','G','deletion','keep'])
	ax.set_yticklabels( ['A','T','C','G','deletion','keep'])

	plt.ylabel('Actual label');
	plt.xlabel('Predicted label');
	all_sample_title = 'Accuracy Score: {0}'.format(accuracy_score(Y.tolist(),result)*100)
	plt.title(all_sample_title, size = 15);
	plt.savefig(confusion)


	sub = X
	prediction = path + '/result.feather'
	DEBUG = path + '/debug.csv'
	RIGHT = path + '/right.csv'
	sub['position'] = position

	sub['label'] = Y

	sub['predict'] = result

	sub.to_feather(prediction)
	debug = sub[sub['predict'] != sub['label']]
	debug.to_csv(DEBUG)
	right = sub[sub['predict'] == sub['label']]
	right.to_csv(RIGHT)
	return prediction

