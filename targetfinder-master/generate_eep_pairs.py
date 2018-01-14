import pandas as pd

for i in ['K562', 'GM12878', 'NHEK', 'IMR90', 'HeLa-S3', 'HUVEC']:
    training_df = pd.read_hdf('%s/output-eep/training.h5' % i, 'training').set_index('window_name')
    pairs = training_df.loc[:,['enhancer_chrom', 'enhancer_start', 'enhancer_end', 'promoter_chrom', 'promoter_start', 'promoter_end', 'label', 'promoter_name', 'enhancer_name']]
    pairs.to_csv('%s/output-eep/pairs.promoter_name.bedpe'%i, sep='\t', header=False, index=False)
