# Create csv averaged matrices

import pandas as pd
import numpy as np

def samples_cov(exp_mat):
    new_exp_mat = pd.DataFrame(columns=list(exp_mat.columns)+['new_index'])
    for i in exp_mat['celltype'].unique():
        n = exp_mat['celltype'].value_counts()[i]
        t = int(n/1)
        for j in range(t):
            sampled = exp_mat[exp_mat['celltype']==i].sample(1)
            indeces2remove = sampled.index
            new_sample_avg = sampled.drop(columns=['celltype','disease']).mean()
            new_sample_avg['celltype'] = i
            new_sample_avg['disease'] = 'Y'
            new_sample_avg['new_index'] = "avg3_" + str(i) + "_" + str(j)
            new_exp_mat = new_exp_mat.append(new_sample_avg,ignore_index=True)
            exp_mat = exp_mat.drop(index=indeces2remove)
        print(f"{i} listo")
    new_exp_mat = new_exp_mat.set_index(['new_index'])
    return new_exp_mat
    
def samples_control(exp_mat):
    new_exp_mat = pd.DataFrame(columns=list(exp_mat.columns)+['new_index'])
    for i in exp_mat['celltype'].unique():
        n = exp_mat['celltype'].value_counts()[i]
        t = int(n/1)
        for j in range(t):
            sampled = exp_mat[exp_mat['celltype']==i].sample(1)
            indeces2remove = sampled.index
            new_sample_avg = sampled.drop(columns=['celltype','disease']).mean()
            new_sample_avg['celltype'] = i
            new_sample_avg['disease'] = 'N'
            new_sample_avg['new_index'] = "avg3_" + str(i) + "_" + str(j)
            new_exp_mat = new_exp_mat.append(new_sample_avg,ignore_index=True)
            exp_mat = exp_mat.drop(index=indeces2remove)
        print(f"{i} listo")
    new_exp_mat = new_exp_mat.set_index(['new_index'])
    return new_exp_mat

exp_mat = pd.read_csv("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/exp_mat_both.csv", index_col=0)

annotations_df = pd.read_csv("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/all.cell.annotation.meta.txt",sep="\t")

exp_mat = exp_mat.loc[exp_mat.index.isin(annotations_df['ID'])]

exp_mat['celltype'] = list(annotations_df['celltype'])
exp_mat['disease'] = list(annotations_df['disease'])

exp_mat_control = exp_mat[exp_mat['disease']=='N']
exp_mat_cov = exp_mat[exp_mat['disease']=='Y']

# cov_ and control_ remove are the number of cells to remove of each cell type in COVID and control

cov_remove = [2,2,1,2,0,0,1,1,1,0]
j = 0
for i in list(exp_mat_cov['celltype'].value_counts().index):
	a = pd.Series(exp_mat_cov[exp_mat_cov['celltype'] == i].index).sample(cov_remove[j])
	if len(list(a)) != 0:
		exp_mat_cov = exp_mat_cov.drop(index=list(a))
	print(f"{i} :\n{a}")
	j += 1

control_remove = [2,2,1,2,1,0,1,1,0]
j = 0
for i in list(exp_mat_control['celltype'].value_counts().index):
	a = pd.Series(exp_mat_control[exp_mat_control['celltype'] == i].index).sample(control_remove[j])
	if len(list(a)) != 0:
		exp_mat_control = exp_mat_control.drop(index=list(a))
	print(f"{i} :\n{a}")
	j += 1


exp_mat_cov_avg3 = samples_cov(exp_mat_cov)

exp_mat_control_avg3 = samples_control(exp_mat_control)

exp_mat_cov_avg3.to_csv('/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/groups/cov/avg3/exp_mat_cov_avg3.csv',sep=',')

exp_mat_control_avg3.to_csv('/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/groups/control/avg3/exp_mat_control_avg3.csv',sep=',')


exp_mat_control = exp_mat_control.drop(columns = ['celltype','disease'])
exp_mat_cov = exp_mat_cov.drop(columns = ['celltype','disease'])

# Create loom files from averaged matrices

import scanpy as sc
import loompy as lp

ex_matrix = exp_mat_control_avg3

row_attrs = {"Gene": np.array(ex_matrix.var.index)}
col_attrs = {"CellID": np.array(ex_matrix.obs.index),"nGene": np.array(np.sum(ex_matrix.X.transpose()>0,axis=0)).flatten(),"nUMI": np.array(np.sum(ex_matrix.X.transpose(),axis=0)).flatten()}

lp.create("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/groups/control/avg3/exp_mat_control_avg3.loom",ex_matrix.X.transpose(), row_attrs, col_attrs)

ex_matrix = exp_mat_cov_avg3

row_attrs = {"Gene": np.array(ex_matrix.var.index)}
col_attrs = {"CellID": np.array(ex_matrix.obs.index),"nGene": np.array(np.sum(ex_matrix.X.transpose()>0,axis=0)).flatten(),"nUMI": np.array(np.sum(ex_matrix.X.transpose(),axis=0)).flatten()}

lp.create("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/groups/cov/avg3/exp_mat_cov_avg3.loom",ex_matrix.X.transpose(), row_attrs, col_attrs)

# Create integrated averaged matrix 

exp_mat_both_avg3 = pd.concat(exp_mat_control_avg3,exp_mat_cov_avg3)

exp_mat_both_avg3.to_csv('/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/exp_mat_both_avg3.csv',sep=',')


ex_matrix = exp_mat_both_avg3

row_attrs = {"Gene": np.array(ex_matrix.var.index)}
col_attrs = {"CellID": np.array(ex_matrix.obs.index),"nGene": np.array(np.sum(ex_matrix.X.transpose()>0,axis=0)).flatten(),"nUMI": np.array(np.sum(ex_matrix.X.transpose(),axis=0)).flatten()}

lp.create("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/exp_mat_avg3.loom",ex_matrix.X.transpose(), row_attrs, col_attrs)
























