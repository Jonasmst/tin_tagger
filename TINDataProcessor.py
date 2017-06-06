import subprocess
import pandas as pd
import os


class TINDataProcessor(object):

    def __init__(self):
        pass

    def load_dataset(self, filepath):
        """
        Reads a dataset and returns it as a pandas dataframe
        """

        # TODO: For certain splice types, start_ex and stop_ex is NaN
        # Specify datatypes
        datatypes = {
            "as_id": "int",
            "included_counts": "int",
            "excluded_counts": "int",
            "psi": "float",
            "psi_denominator": "float",
            "magnitude": "float",
            "name": "object",
            #"start_ex": "int",
            #"end_ex": "int",
            "exons": "object",
            "splice_type": "object",
            "novel_splice": "int",
            "graph_id": "int",
            "symbol": "object",
            "chr": "object",
            "strand": "object",
            "rpkm": "float",
            "exon1": "object",
            "exon2": "object",
            #"prev_exon_stop": "int",
            #"prev_exon_start": "int",
            #"next_exon_start": "int",
            #"next_exon_stop": "int",
            "first_exon_in_splice": "object",
            "last_exon_in_splice": "object",
            "first_exon_start": "int",
            "last_exon_stop": "int"
        }

        # Infer delimiter by file extension
        name, extension = os.path.splitext(filepath)
        delimiter = ","
        if extension.lower() == ".tsv":
            delimiter = "\t"

        # Read into pandas dataframe
        df = pd.read_csv(filepath, dtype=datatypes, sep=delimiter)

        # Trim .0 from exon names
        #df.exon1 = df.exon1.str.replace("\.0$", "")
        #df.exon2 = df.exon2.str.replace("\.0$", "")

        # TODO: Don't do this; return full dataset
        df = df.loc[df.splice_type == "ES"]

        return df


