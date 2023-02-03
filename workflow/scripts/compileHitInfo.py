"""
This code adds info for species that have hits for each gene and adds it to a CSV file
^ i can't make words today
"""

import sys
import pandas as pd

file_path = sys.argv[1]
score_file = sys.argv[2]

# add header to compiled_hits csv file
# header = "description,hit_id,targ_id,target_len,hit_score,hit_evalue,hit_bias,hsp_len,hsp_score,hsp_evalue_cond,hsp_evalue_ind,hsp_bias,output_table \n"
# f = open(file_path, 'r')
# olines = f.readlines()
# olines.insert(0, header)
# f.close()
# f = open(file_path, 'w')
# f.writelines(olines)
# f.close()

# make a csv file with only the highest score value per strain per profile HMM
df = pd.read_csv(file_path)
# make new column for strain name and profile HMM
df[["strain", "gene"]] = df['description'].str.rsplit("_", 1, expand=True)


# df["strain"] = df['description'].str.rsplit("_", 2, expand=True)[0]
# df["gene"] = df['description'].str.rsplit("_", 2, expand=True)[1] + \
#     "_" + df['description'].str.rsplit("_", 2, expand=True)[2]

# score_df = df.groupby(["strain", "gene"], as_index=False).hit_score.max()
score_df = df.loc[df.groupby(["strain", "gene"])["hit_score"].idxmax()].reset_index(drop=True)
score_df = score_df[["strain", "gene", "hit_score", "hit_id", "description"]]
score_df.to_csv(score_file)
