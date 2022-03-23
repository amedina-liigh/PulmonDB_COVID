import pandas as pd

# SI TIENES LOS DATOS DE LOS REGULONES EN EL FORMATO PKL, PUEDES CARGARLOS ASÍ
if pkl_format == 1:
    import pickle
    import gzip

    with gzip.open(path_regulons_targets_pkl, 'rb') as f:
        regs_signatures = pickle.load(f)
    f.close()
else :
    #SI TIENES LOS DATSO EN LOS .CSV, O SEA, EL DE REG.CSV, PUEDES HACER ESTO
    from pyscenic.cli.utils import load_signatures
    regs_signatures = load_signatures(path_regulons_targets_csv)

    # ESTE PASO PUEDE TARDAR UN BASTANTE

# USANDO LAS VARIABLES QUE CREÉ EN R, OBTENGO LOS REGULONES QUE QUIERO
indexes_of_regs = {}
for i,j in enumerate(regs_signatures):
  indexes_of_regs[j.name] = i

# ESTOS SON TOMANDO EN CUENTA LOS RESULTADOS DEL WILCONXON TEST, EL KOLMOGOROV TEST Y EL RSS
TF_targets = {}
for i in r.ccmk.keys():
  TF_targets[i] = {}
  for j in r.ccmk[i]:
    TF_targets[i][j] = regs_signatures[indexes_of_regs[j]].genes

# ESTOS SON TOMANDO EN CUENTA LOS RESULTADOS DEL WILCONXON TEST, EL KOLMOGOROV TEST Y EL PSEUDO-FOLD CHANGE
TF_targets2 = {}
for i in r.lfc.keys():
  TF_targets2[i] = {}
  for j in r.lfc[i]:
    TF_targets2[i][j] = regs_signatures[indexes_of_regs[j]].genes
