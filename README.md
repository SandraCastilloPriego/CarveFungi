# CarveFungi

CarveFungi is a genome-scale metabolic model reconstruction pipeline able to create a compartmentalized metabolic model of any fungal specie from its protein sequences. It is implemented using python.

Requirements:
- Tensorflow
- keras
- biopython
- pandas
- numpy
- cobrapy
- Framed (https://github.com/cdanielmachado/framed)


 CarveFungi uses a deep learning model to predict the cellular localization of the proteins and this information is used to score the reactions. 
![Deep neural network](/images/CNN.png)

CarveMe (https://doi.org/10.1093/nar/gky537,https://github.com/cdanielmachado/carveme) uses the scoring of the reactions to obtain a functional metabolic model that is able to produce biomass. 

