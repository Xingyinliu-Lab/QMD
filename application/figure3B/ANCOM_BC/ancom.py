from ancomP.stats.ancom import ancom
from ancomP.linalg import composition
import pandas as pd


genusdata=pd.read_csv('PRJNA763023_single.csv',index_col=None,header=0)
taxalist=list(genusdata.columns)
taxalist.remove('group')

sumabundance = genusdata[taxalist].sum(axis=1)
genusdata[taxalist] = genusdata[taxalist].div(
    sumabundance, axis='rows')


table=genusdata[taxalist]
table=table.fillna(value=0)
table=composition.multiplicative_replacement(table)
table=pd.DataFrame(table)
table.index=genusdata.index
table.columns=taxalist
labels, uniques = pd.factorize(genusdata['group'])

grouping = pd.Series(labels,index=genusdata.index)
permu=500
results = ancom(table, grouping, significance_test='permutative-anova', permutations=permu)
results.to_csv('ancom.csv')
pass
