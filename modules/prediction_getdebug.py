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

def predict_class(input_data, model_FILE):

	model = joblib.load(model_FILE).set_params(n_jobs=1)
	p_label=model.predict(input_data.values.tolist())
	return p_label


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

	pool = Pool(32)
	results = pool.map(partial(predict_class, model_FILE=model_file), list_of_X)
	pool.close()
	pool.join()

	result = []
	for i in results:
		result.extend(i)
	print(len(result))

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
