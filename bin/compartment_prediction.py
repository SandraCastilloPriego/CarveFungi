import tensorflow as tf
import numpy as np
from Bio import SeqIO
from tensorflow.keras.models import load_model
import os


def run_compartment_predictions(path, fasta_file):
    name =fasta_file.split("/")
    name = name[-1]

    if os.path.isfile(path+name+".loc_pred"):
        return path+name+".loc_pred"

    fasta_id=[]
    fasta_seq=[]
    print(fasta_file)
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_id.append(record.id)
        fasta_seq.append(str(record.seq))


    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    aa_to_int = {aa: i+1 for i, aa in enumerate(amino_acids)}  # +1 to reserve 0 for padding

    chem_feat = {"A": [0.382, 0.389, 0.400, 0.623, 0.112, 0.559, 0.000],
                "E": [0.693, 0.608, -0.778, 0.946, 0.369, 0.299, 0.335],
                "K": [0.689, 0.734, -0.867, 0.869, 0.535, 0.905, 0.120],
                "L": [0.608, 0.710, 0.844, 0.377, 0.455, 0.556, 0.000],
                "D": [0.618, 0.505, -0.778, 1.000, 0.257, 0.257, 0.502],
                "R": [0.839, 0.832, -1.000, 0.808, 0.711, 1.000, 0.236],
                "F": [0.791, 0.835, 0.622, 0.400, 0.709, 0.509, 0.000],
                "I": [0.608, 0.712, 1.000, 0.400, 0.455, 0.559, 0.000],
                "Q": [0.688, 0.645, -0.778, 0.808, 0.440, 0.525, 0.324],
                "W": [1.000, 1.000, -0.200, 0.415, 1.000, 0.547, 0.047],
                "G": [0.306, 0.275, -0.089, 0.692, 0.000, 0.555, 0.269],
                "S": [0.468, 0.407, -0.178, 0.708, 0.152, 0.528, 0.516],
                "H": [0.736, 0.688, -0.711, 0.800, 0.562, 0.705, 0.211],
                "N": [0.613, 0.550, -0.778, 0.892, 0.328, 0.503, 0.484],
                "M": [0.705, 0.724, 0.422, 0.438, 0.540, 0.533, 0.000],
                "P": [0.521, 0.531, -0.356, 0.615, 0.320, 0.602, 0.142],
                "Y": [0.876, 0.851, -0.289, 0.477, 0.729, 0.526, 0.073],
                "C": [0.554, 0.447, 0.556, 0.423, 0.313, 0.471, 1.000],
                "V": [0.532, 0.600, 0.933, 0.454, 0.342, 0.555, 0.000],
                "T": [0.543, 0.518, -0.156, 0.662, 0.264, 0.546, 0.258],
                "X": [0.639, 0.628, -0.109, 0.640, 0.430, 0.562, 0.226],
                "Z": [0.691, 0.627, -0.778, 0.877, 0.405, 0.412, 0.33],
                "B": [0.616, 0.528, -0.778, 0.946, 0.293, 0.38, 0.493], "J": [0.608, 0.711, 0.922, 0.389, 0.455, 0.558, 0],
                "U": [0.554, 0.447, 0.556, 0.423, 0.313, 0.471, 1.000],
                "O": [0.689, 0.734, -0.867, 0.869, 0.535, 0.905, 0.120], "ERROR": [0, 0, 0, 0, 0, 0, 0]}


    def get_chemical_features(aa):
        sequence = np.array([])
        for i in range(750):
            try:
                if sequence.shape[0] == 0:
                    sequence = np.array(chem_feat[aa[i]]).reshape(1, 7)
                else:
                    sequence = np.concatenate((sequence, np.array(chem_feat[aa[i]]).reshape(1, 7)))
            except:
                if sequence.shape[0] == 0:
                    sequence = np.array(chem_feat["ERROR"]).reshape(1, 7)
                else:
                    sequence = np.concatenate((sequence, np.array(chem_feat["ERROR"]).reshape(1, 7)))

        return sequence.reshape(1, 750, 7)
    
    

    def encode_sequence(seq, aa_to_int, max_len=750):
        encoded = [aa_to_int.get(aa, 0) for aa in seq[:max_len]]
        return encoded + [0] * (max_len - len(encoded))

    # Example: sequences = list of strings, labels = list of 8-element vectors
    encoded_seqs = [encode_sequence(seq, aa_to_int) for seq in fasta_seq]
    features= [get_chemical_features(aa) for aa in fasta_seq]

    encoded_seqs = tf.convert_to_tensor(encoded_seqs, dtype=tf.int32)
    features = tf.convert_to_tensor(features, dtype=tf.float32)

    features= tf.squeeze(features, axis=1)
    try:
        model = load_model('/home/scsandra/Projects/FoodID/CarveFungi/bin/model_binary_output10_18.keras')
    except Exception as e:
        print(f"Error loading model: {e}")
        return None

    # Load test data
    X_test = np.array(encoded_seqs)

    # Predict probabilities
    y_pred = model.predict([X_test,features])

    import pandas as pd
    newDF = pd.DataFrame()
    E = []
    M = []
    P = []
    O = []
    Ex = []
    N = []
    Mb = [] 
    Gg = [] 

    for i in range(len(fasta_id)):
        pred= y_pred[i]
        E.append(pred[3])
        M.append(pred[2])
        P.append(pred[5])
        Ex.append(pred[6])
        N.append(pred[1])
        Mb.append(pred[7])
        O.append(pred[0])
        Gg.append(pred[4])

    newDF["E"]=E
    newDF["M"]=M
    newDF["P"]=P
    newDF["O"]=O
    newDF["Ex"]=Ex
    newDF["N"]=N
    newDF["Mb"]=Mb
    newDF["Gg"]=Gg
    newDF["Ids"]=fasta_id
    
    newDF.to_csv(path+name+".loc_pred")
    return path+name+".loc_pred"