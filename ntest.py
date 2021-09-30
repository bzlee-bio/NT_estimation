import ntox_estimator
import tensorflow as tf
from Bio import SeqIO
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Neurotoxicity estimation.')
parser.add_argument('--fasta',metavar='fasta', type=str, help='Fasta file')
parser.add_argument('--output',metavar='output', type=str, help='Output file')
args = parser.parse_args()

print('Fasta file loading')
fasta_p = SeqIO.parse(args.fasta,'fasta')

data = []
## Data filtering
count = 0
for l in fasta_p:
    count+=1    
    if len(l.seq)>=50 and len(l.seq)<=300:

        data.append(l)

print('%d of %d sequences were loaded '%(len(data),count))        
print('Sequence converting')
so = ntox_estimator.SeqtoOHE()
(arr_ohe,id_ohe,seq) = so.gen_seq_np(data)

print('Estimation')
print(arr_ohe.shape)
model = tf.keras.models.load_model('./ntox_estimator/model_file/ntox.h5')
pred_prob = model.predict(arr_ohe[:,:,:,np.newaxis])
pred_res = np.argmax(pred_prob, 1)
# predict_x=model.predict(X_test) 
# classes_x=np.argmax(predict_x,axis=1)
print('==============================================================')
print('%d of %d sequences were predicted as neurotoxin'%(np.sum(pred_res),len(data)))

## Write results
filtered_list = [[i,j] for (i,j, v) in zip(id_ohe,seq, pred_res) if v==1]
with open(args.output,'w',newline="") as f:
    f.writelines('## Predicted neurotoxins list')
    f.writelines('ID, Sequence\n')
    for l in filtered_list:
        f.writelines(l[0]+', '+l[1]+'\n')
        
print('Predicted neurotoxins were saved into %s'%args.output)
print('Done')