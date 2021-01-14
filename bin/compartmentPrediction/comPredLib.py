import numpy as np
aa = ["A", "E", "K", "L", "D", "R", "F", "I", "Q", "W", "G", "S", "H", "N", "M", "P", "Y", "C", "V", "T", "X", "Z", "B",
      "U", "O"]
Mass = [0.382, 0.693, 0.688, 0.608, 0.618, 0.839, 0.790, 0.608, 0.688, 1.000, 0.306, 0.468, 0.736, 0.613, 0.705, 0.522,
        0.876, 0.554, 0.532, 0.543, 0.591, 0.591, 0.591, 0.746, 0.608]
Hydrophobicity = [0.616, 0.043, 0.283, 0.943, 0.028, 0.00, 1.00, 0.943, 0.251, 0.878, 0.501, 0.359, 0.165, 0.236, 0.738,
                  0.711,
                  0.880, 0.680, 0.825, 0.450, 0.5, 0.5, 0.5, 0.5, 0.5]
Polarity = [8.1, 12.3, 11.3, 4.9, 13, 10.5, 5.2, 5.2, 10.5, 5.4, 9, 9.2, 10.4, 11.6, 5.7, 8, 6.2, 5.5, 5.9, 8.6, 8.3,
            8.3, 8.3, 8.3, 8.3]
COOHpka = [2.34, 2.19, 2.18, 2.36, 1.88, 2.17, 1.83, 2.36, 2.17, 2.83, 2.34, 2.21, 1.82, 2.02, 2.28, 1.99, 2.2, 1.96,
           2.32, 2.09, 2.177, 2.177, 2.177, 2.177, 2.177]
NH3pka = [9.69, 9.67, 8.95, 9.6, 9.6, 9.04, 9.13, 9.6, 9.13, 9.39, 9.6, 9.15, 9.17, 8.8, 9.21, 10.6, 9.11, 10.128, 9.62,
          9.1, 9.4144, 9.4144, 9.4144, 9.4144, 9.4144]
pI = [6, 3.22, 9.74, 5.98, 2.77, 10.76, 5.48, 6.02, 5.65, 5.89, 5.97, 5.58, 7.59, 5.41, 5.74, 6.3, 5.66, 5.07, 5.96,
      5.6, 6.0195, 6.0195, 6.0195, 6.0195, 6.0195]
SASA = [1.181, 1.862, 2.258, 1.931, 1.587, 2.56, 2.228, 1.81, 1.932, 2.663, 0.881, 1.298, 2.025, 1.655, 2.034, 1.468,
        2.368, 1.461, 1.645, 1.525, 1.8186
    , 1.8186, 1.8186, 1.8186, 1.8186]
NCISC = [0.07187, 0.006802, 0.017708, 0.051672, -0.02382, 0.043587, 0.037552, 0.21631, 0.049211, 0.037977, 0.179052,
         0.004627, -0.01069, 0.005392, 0.002683, 0.239531, 0.023599, -0.03661, 0.057004, 0.003352, 0.0488404, 0.0488404,
         0.0488404, 0.0488404, 0.0488404]
aa_dict = {}
for i in range(len(aa)):
    aa_dict[aa[i]] = Mass[i]

aa_dict2 = {}
for i in range(len(aa)):
    aa_dict2[aa[i]] = Hydrophobicity[i]

aa_dict3 = {}
for i in range(len(aa)):
    aa_dict3[aa[i]] = Polarity[i]

aa_dict4 = {}
for i in range(len(aa)):
    aa_dict4[aa[i]] = COOHpka[i]
aa_dict5 = {}
for i in range(len(aa)):
    aa_dict5[aa[i]] = NH3pka[i]
aa_dict6 = {}
for i in range(len(aa)):
    aa_dict6[aa[i]] = pI[i]
aa_dict7 = {}
for i in range(len(aa)):
    aa_dict7[aa[i]] = SASA[i]
aa_dict8 = {}
for i in range(len(aa)):
    aa_dict8[aa[i]] = NCISC[i]



def getAAfeatures(g):
    sdataIdx = [aa_dict[c] for c in list(g)]
    sdataIdx2 = [aa_dict2[c] for c in list(g)]
    sdataIdx3 = [aa_dict3[c] for c in list(g)]
    sdataIdx4 = [aa_dict4[c] for c in list(g)]
    sdataIdx5 = [aa_dict5[c] for c in list(g)]
    sdataIdx6 = [aa_dict6[c] for c in list(g)]
    sdataIdx7 = [aa_dict7[c] for c in list(g)]
    sdataIdx8 = [aa_dict8[c] for c in list(g)]
    frequencies = np.asarray(sdataIdx)
    # print("1: ", frequencies.shape)
    frequencies = np.reshape(frequencies, (1, len(g)))
    sdataIdx2 = np.reshape(np.asarray(sdataIdx2), (1, len(g)))
    sdataIdx3 = np.reshape(np.asarray(sdataIdx3), (1, len(g)))
    sdataIdx4 = np.reshape(np.asarray(sdataIdx4), (1, len(g)))
    sdataIdx5 = np.reshape(np.asarray(sdataIdx5), (1, len(g)))
    sdataIdx6 = np.reshape(np.asarray(sdataIdx6), (1, len(g)))
    sdataIdx7 = np.reshape(np.asarray(sdataIdx7), (1, len(g)))
    sdataIdx8 = np.reshape(np.asarray(sdataIdx8), (1, len(g)))
    frequencies = np.concatenate((frequencies, sdataIdx2), axis=0)
    frequencies = np.concatenate((frequencies, sdataIdx3), axis=0)
    frequencies = np.concatenate((frequencies, sdataIdx4), axis=0)
    frequencies = np.concatenate((frequencies, sdataIdx5), axis=0)
    frequencies = np.concatenate((frequencies, sdataIdx6), axis=0)
    frequencies = np.concatenate((frequencies, sdataIdx7), axis=0)
    frequencies = np.concatenate((frequencies, sdataIdx8), axis=0)
    # print(len(g),  frequencies.shape)
    # frequencies = np.reshape( frequencies, (len(g), 8))
    # print(frequencies.shape)
    return frequencies