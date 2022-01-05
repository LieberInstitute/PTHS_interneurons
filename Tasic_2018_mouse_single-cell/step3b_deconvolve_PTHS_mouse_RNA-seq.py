import pandas as pd
import numpy as np
from sklearn.svm import LinearSVR
from sklearn.datasets import make_regression
import matplotlib.pyplot as plt
import random
from matplotlib.pyplot import figure

reference = pd.read_csv("tables/Allen_Mouse_Brain_Reference_logCounts_p1-Adult.csv")
test_data = pd.read_csv("tables/PTHSmouse_samples_logCounts_p1-Adult.csv")

train = reference.drop("Unnamed: 0",1)
test_data = test_data.drop("Unnamed: 0",1)

score_adj = []
for o in range(len(test_data.columns)):
	print('Testing ' + str(o) + ' of ' + str(len(test_data.columns)) + ' runs.')
	test = test_data.loc[:,test_data.columns[o]]
	im_name = train.columns
	svr = LinearSVR(random_state=0, max_iter = 100000)
	model = svr.fit(train, test)
	score = model.coef_
	score[np.where(score<0)] = 0 
	score_adj.append((score/sum(score)))

df = pd.DataFrame(score_adj, index=test_data.columns, columns=im_name)
df.to_csv("tables/PTHSmouse_logCounts_deconvolution_score_p1-Adult_maxIter10000.csv")



