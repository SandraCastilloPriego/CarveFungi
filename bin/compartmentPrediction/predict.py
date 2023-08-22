import numpy as np
import tensorflow as tf
import copy
import math
from Bio import SeqIO
import numpy as np
import sys, getopt
from comPredLib import get_chemical_features,get_binary_sequence
import pandas as pd

model = tf.keras.models.load_model('selected_models/model-all')
model1 = tf.keras.models.load_model('selected_models/model-bidireccional-8feat-59')
model2 = tf.keras.models.load_model('selected_models/model-bidireccional-8feat-gru-32-59')
model3 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-gru-32-59')
model4 = tf.keras.models.load_model('selected_models/model-bidireccional-8feat-sec-gru-32-59')
model5 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-sec-gru-32-59')
model6 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-feat7-gru-32-59')
model7 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-feat7-sec-gru-32-59')
model8 = tf.keras.models.load_model('selected_models/model-bidireccional-8feat-gru-cat-32-49')
model9 = tf.keras.models.load_model('selected_models/model-bidireccional-8feat-sec-gru-cat-32-50')
model10 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-gru-32-cat-170')
model11 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-gru-32-cat-mask160')
model12 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-feat7-gru-32-mask.59')
model13 = tf.keras.models.load_model('selected_models/model-bidireccional-onehot-feat7-gru-cat-32-mask-59')
model14 = tf.keras.models.load_model('selected_models/model-transformer-onehot-feat7-32-150')
model15 = tf.keras.models.load_model('selected_models/model-transformer-onehot-sec-32-100')
model16 = tf.keras.models.load_model('selected_models/model-transformer-onehot-32-300')
model17 = tf.keras.models.load_model('selected_models/model-transformer-onehot-32-100')
model18 = tf.keras.models.load_model('selected_models/model-transformer-8feat-sec-32-60')
model19 = tf.keras.models.load_model('selected_models/model-transformer-bidireccional-8feat-32-200')


class data_training:
    entry_id = None
    data3d = None
    aa = None
    onehot=None
    function = None
    sec_features =None

    size = 500
    def __init__(self, entry_id, path):
        self.aa, self.onehot, self.sec_features = self.get_aa_data(entry_id, path)
        self.entry_id = entry_id

    def process_secondary_structure_files(self, sec_file, length):
        try:
            matrix_s = np.empty((0, 3), int)
            secC = np.zeros([length])
            secH = np.zeros([length])
            secB = np.zeros([length])
            secfile = open(sec_file, "r")
            i = 0
            for line in secfile:
                if "Pred: " in line:
                    line = line.replace("Pred: ", "")
                    line = line.replace(" ", "")
                    for letter in line:
                        if 'C' == letter:
                            secC[i] = 1
                            secH[i] = 0
                            secB[i] = 0
                            i = i + 1
                        elif 'H' == letter:
                            secC[i] = 0
                            secH[i] = 1
                            secB[i] = 0
                            i = i + 1
                        elif 'E' == letter:
                            secC[i] = 0
                            secH[i] = 0
                            secB[i] = 1
                            i = i + 1
            matrix_s = np.append(secC, secH)
            matrix_s = np.append(matrix_s, secB).reshape(1, 3, length)
            matrix_s = np.swapaxes(matrix_s, 1, 2)
        except:
            raise
            print("No sec")
        return matrix_s

    def get_aa_data(self, entry_id, path):
        fastafile = open(str(path) + entry_id + ".fasta", "r")
        g = []
        for data in fastafile:
            if '>' in data:
                continue
            for aa in data:
                if '\n' in aa:
                    continue
                g.append(aa)

        sec_struct = self.process_secondary_structure_files(str(path) + entry_id+".horiz", len(g))
        if sec_struct.shape[1] > 500:
            partial_sec_struct = np.concatenate(
                (sec_struct[:, 1:251, :], sec_struct[:, sec_struct.shape[1] - 250:sec_struct.shape[1], :]), axis=1)
            sec_struct = partial_sec_struct
        elif sec_struct.shape[1] < 500:
            size = 500 - sec_struct.shape[1]
            missing_part = np.zeros([1, size, 3])
            sec_struct = np.concatenate((sec_struct, missing_part), axis=1)
        #print(sec_struct.shape)

        sec_struct = np.array(sec_struct.reshape(1, 500, 3), dtype=bool)


        chem_features = get_chemical_features(g)
        #print(chem_features.shape)
        if len(chem_features) > 500:
            partial_chem_features = np.concatenate(
                (chem_features[1:251], chem_features[len(chem_features) - 250:len(chem_features)]), axis=0)
            chem_features = partial_chem_features
        elif sec_struct.shape[0] < 500:
            size = 500 - chem_features.shape[0]
            missing_part = np.zeros([size, 7])
            chem_features = np.concatenate((chem_features, missing_part), axis=0)

        chem_features = np.array(chem_features.reshape(1, 500, 7), dtype=np.float16)
        #print(chem_features.shape)
        sequence = np.array(get_binary_sequence(g), dtype=bool)

        if len(sequence) > 500:
            partial_sequence = np.concatenate((sequence[1:251], sequence[len(sequence) - 250:len(sequence)]), axis=0)
            sequence = partial_sequence
        elif sequence.shape[0] < 500:
            size = 500 - sequence.shape[0]
            missing_part = np.zeros([size, 21])
            sequence = np.concatenate((sequence, missing_part), axis=0)
        sequence = np.array(sequence.reshape(1, 500, 21), dtype=bool)
        print(sequence.shape)
        return chem_features, sequence, sec_struct


    def get_id(self):
        return self.entry_id

    def get_aa(self):
        return self.aa, self.onehot

    def get_sec(self):
        return self.sec_features



def getPrediction(chem_feat, onehot, sec_feat):
    predictions1 =model1([chem_feat])
    predictions2 =model2([chem_feat])
    predictions3 =model3([onehot])
    predictions4 =model4([chem_feat,sec_feat])
    predictions5 =model5([onehot,sec_feat])
    predictions6 = model6([onehot, chem_feat])
    predictions7 = model7([onehot,chem_feat, sec_feat])
    predictions8 = model8([chem_feat])
    predictions9 = model9([chem_feat,sec_feat])
    predictions10 = model10([onehot])
    predictions11 = model11([onehot])
    predictions12 = model12([onehot, chem_feat])
    predictions13 = model13([onehot, chem_feat])
    predictions14 = model14([onehot, chem_feat])
    predictions15 = model15([onehot, sec_feat])
    predictions16 = model16([onehot])
    predictions17 = model17([onehot])
    predictions18 = model18([chem_feat,sec_feat])
    predictions19 = model19([chem_feat])

    #print(predictions1.shape)
    x = [predictions1.numpy().reshape(chem_feat.shape[0], 4),predictions2.numpy().reshape(chem_feat.shape[0], 4),predictions3.numpy().reshape(chem_feat.shape[0], 4),predictions4.numpy().reshape(chem_feat.shape[0], 4),predictions5.numpy().reshape(chem_feat.shape[0], 4),predictions6.numpy().reshape(chem_feat.shape[0], 4),predictions7.numpy().reshape(chem_feat.shape[0], 4),predictions8.numpy().reshape(chem_feat.shape[0], 4),predictions9.numpy().reshape(chem_feat.shape[0], 4),predictions10.numpy().reshape(chem_feat.shape[0], 4),predictions11.numpy().reshape(chem_feat.shape[0], 4),predictions12.numpy().reshape(chem_feat.shape[0], 4),predictions13.numpy().reshape(chem_feat.shape[0], 4),predictions14.numpy().reshape(chem_feat.shape[0], 4),predictions15.numpy().reshape(chem_feat.shape[0], 4),predictions16.numpy().reshape(chem_feat.shape[0], 4),predictions17.numpy().reshape(chem_feat.shape[0], 4),predictions18.numpy().reshape(chem_feat.shape[0], 4),predictions19.numpy().reshape(chem_feat.shape[0], 4)]
    x = np.concatenate(x, axis=1)
    y = model(x)
    return y.numpy()

def predict(organismPath, organism, suma):
    resultPath = '/scratch/project_2001899/CarveFungi/data/annotations/localization/'

    functional_annotation = pd.read_csv('/scratch/project_2001899/CarveFungi/data/annotations/functional/'+organism, skiprows=4, sep="\t",header=None)
    functional_annotation = functional_annotation[functional_annotation[7].notnull()]
    chem_features= np.array([])
    sec_features = np.array([])
    ids=[]
    for key, row in functional_annotation.iterrows():
        train = data_training(row[0],organismPath)
        ids.append(row[0])
        if chem_features.shape[0]==0:
            chem_features, onehot = train.get_aa()
            sec_features=train.get_sec()
        else:
            chem_features_o, onehot_o =train.get_aa()
            chem_features = np.concatenate([chem_features,chem_features_o])
            onehot = np.concatenate([onehot,onehot_o])
            sec_features_conc= train.get_sec()
            sec_features = np.concatenate([sec_features,sec_features_conc])

    print("predicting")

    predictions = getPrediction(chem_features, onehot, sec_features)
    
    i = 0
    for p in predictions:
        if suma[i]:
            p[2]=p[2]+0.05
        i=i+1

    df = pd.DataFrame(predictions, columns=["E", "M", "P", "O"])
    df["Ids"] = ids
    df.to_csv(resultPath+organism+".loc_pred")


