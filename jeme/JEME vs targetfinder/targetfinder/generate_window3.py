import pandas as pd

def write_window(cell):
    data = pd.read_csv(cell + "/training.csv")
    #out_data = data.ix[:, ['window_chrom', 'window_start', 'window_end']]
    #out_data.to_csv(cell + "/output-epw/windows.bed3", sep='\t', header=False, index=False)
    out_data = data.ix[data['label']==1, ['window_chrom', 'window_start', 'window_end', 'promoter_name', 'enhancer_name']]
    out_data.to_csv(cell + "/windows.promoter_name.pos.bed", sep='\t', header=False, index=False)
    out_data = data.ix[data['label']==0, ['window_chrom', 'window_start', 'window_end', 'promoter_name', 'enhancer_name']]
    out_data.to_csv(cell + "/windows.promoter_name.neg.bed", sep='\t', header=False, index=False)
    #out_data = data.ix[:, ["window_chrom", "window_start", "window_end", "window_name", "label", "active_promoters_in_window", "enhancer_distance_to_promoter"]]
    #out_data.to_csv(cell + "/output-epw/windows.bed", sep='\t', header=False, index=False)

for cell in ['GM12878', 'K562']:
    write_window(cell)

