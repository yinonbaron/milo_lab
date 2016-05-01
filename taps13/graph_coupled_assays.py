import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file

reactions = pd.DataFrame.from_csv('../res/nadh_forming_enzymes.tsv', sep='\t')
