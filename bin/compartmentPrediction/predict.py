import numpy as np
from tensorflow import keras
import copy
import math
from Bio import SeqIO
import numpy as np
import sys, getopt
from comPredLib import getAAfeatures


json_file = open('bin/compartmentPrediction/deepModels/model_train1_6.json', 'r')
model_json = json_file.read()
json_file.close()
model = keras.models.model_from_json(model_json)
model.load_weights("bin/compartmentPrediction/deepModels/model_train1_6.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train1_99.json', 'r')
model_json = json_file.read()
json_file.close()
modelb = keras.models.model_from_json(model_json)
modelb.load_weights("bin/compartmentPrediction/deepModels/model_train1_99.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train2_5.json', 'r')
model_json = json_file.read()
json_file.close()
model2 = keras.models.model_from_json(model_json)
model2.load_weights("bin/compartmentPrediction/deepModels/model_train2_5.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train2_99.json', 'r')
model_json = json_file.read()
json_file.close()
model2b = keras.models.model_from_json(model_json)
model2b.load_weights("bin/compartmentPrediction/deepModels/model_train2_99.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train3_5.json', 'r')
model_json = json_file.read()
json_file.close()
model3 = keras.models.model_from_json(model_json)
model3.load_weights("bin/compartmentPrediction/deepModels/model_train3_5.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train3_99.json', 'r')
model_json = json_file.read()
json_file.close()
model3b = keras.models.model_from_json(model_json)
model3b.load_weights("bin/compartmentPrediction/deepModels/model_train3_99.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train4_55.json', 'r')
model_json = json_file.read()
json_file.close()
model4 = keras.models.model_from_json(model_json)
model4.load_weights("bin/compartmentPrediction/deepModels/model_train4_55.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train4_99.json', 'r')
model_json = json_file.read()
json_file.close()
model4b = keras.models.model_from_json(model_json)
model4b.load_weights("bin/compartmentPrediction/deepModels/model_train4_99.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train4_55.json', 'r')
model_json = json_file.read()
json_file.close()
model5 = keras.models.model_from_json(model_json)
model5.load_weights("bin/compartmentPrediction/deepModels/model_train5_55.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train5_99.json', 'r')
model_json = json_file.read()
json_file.close()
model5b = keras.models.model_from_json(model_json)
model5b.load_weights("bin/compartmentPrediction/deepModels/model_train5_99.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train4_55.json', 'r')
model_json = json_file.read()
json_file.close()
model6 = keras.models.model_from_json(model_json)
model6.load_weights("bin/compartmentPrediction/deepModels/model_train6_55.h5")

json_file = open('bin/compartmentPrediction/deepModels/model_train6_99.json', 'r')
model_json = json_file.read()
json_file.close()
model6b = keras.models.model_from_json(model_json)
model6b.load_weights("bin/compartmentPrediction/deepModels/model_train6_99.h5")



def getPrediction(datax):
    predictions =model.predict(datax)
    predictionsb =modelb.predict(datax)
    predictions2 =model2.predict(datax)
    predictions2b =model2b.predict(datax)
    predictions3 =model3.predict(datax)
    predictions3b =model3b.predict(datax)
    predictions4 =model4.predict(datax)
    predictions4b =model4b.predict(datax)
    predictions5 =model5.predict(datax)
    predictions5b =model5b.predict(datax)
    predictions6 =model6.predict(datax)
    predictions6b =model6b.predict(datax)

    average = np.average([predictions,predictions2,predictions3,predictions4,predictions5,predictions6,predictionsb,predictions2b,predictions3b,predictions4b,predictions5b,predictions6b], axis=0)
    return average

def predict(organismPath):	
	count2=0
	resultPath = organismPath.replace('fasta/','localization/')
	name = organismPath[0:]

	datax = np.array([])
	predictionsArray = np.array([])
	count = 0
	print("predicting")
	for record in SeqIO.parse(organismPath, "fasta"):
			count2=count2+1			
			if count > 1000:
					print("Predicted: "+str(count2)+" proteins")
					datax = datax.reshape(datax.shape[0],datax.shape[1],datax.shape[2],1)
					predictions = getPrediction(datax)
					if predictionsArray.shape[0]==0:
							predictionsArray = predictions
					else:
							predictionsArray = np.concatenate((predictionsArray, predictions), axis=0)
					datax = np.array([])							
					count = 0

			g = str(record.seq)
			g = g.replace("*","")
			g = g.replace("J","")
			
			try:
					features = getAAfeatures(g)
						#print(features.shape)
					if len(g) > 2000:
							features = np.concatenate((features[:,0:1000],features[:,features.shape[1]-1000:features.shape[1]]), axis=1)
							features = np.reshape(features, (1, 8, 2000))
					else:            
							size = ((2002 - len(g)) / 2)
							npad = ((int(size), int(size)), (0, 0))
							sdataPadded = np.pad(np.transpose(features), pad_width=npad, mode="constant", constant_values=0)
							features = np.transpose(sdataPadded)[:,0:2000]
							features = np.reshape(features, (1, 8, 2000))    

					if datax.shape[0] == 0:
							datax = features
					else:
							datax = np.concatenate((datax, features), axis=0)
					count = count + 1
			except:
				print(datax.shape, features.shape)
				raise


	datax = datax.reshape(datax.shape[0],datax.shape[1],datax.shape[2],1)        
	predictions = getPrediction(datax)
	if predictionsArray.shape[0] == 0:
	    predictionsArray = predictions
	else:
	    predictionsArray = np.concatenate((predictionsArray, predictions), axis=0)


	np.save(resultPath, predictionsArray)		
	
				
