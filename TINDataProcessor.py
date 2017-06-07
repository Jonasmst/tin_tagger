import subprocess
import pandas as pd
import os
import json


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
            "first_exon_in_splice": "object",
            "last_exon_in_splice": "object",
            "first_exon_start": "int",
            "last_exon_stop": "int",
            #"prev_exon_stop": "int",
            #"prev_exon_start": "int",
            #"next_exon_start": "int",
            #"next_exon_stop": "int",
            # "start_ex": "int",
            # "end_ex": "int",
        }

        # Infer delimiter by file extension
        name, extension = os.path.splitext(filepath)
        delimiter = ","
        if extension.lower() == ".tsv":
            delimiter = "\t"

        # Read into pandas dataframe
        df = pd.read_csv(filepath, dtype=datatypes, sep=delimiter)

        # TODO: Don't do this; return full dataset
        df = df.loc[df.splice_type == "ES"]

        return df

    def filter_dataset(self, dataset, filters):
        """
        Takes a Pandas Dataframe in param dataset and filter it based on criteria given in param filters.
        Returns a Pandas DataFrame if there are datapoints in the remaining filtered dataset, otherwise returns None.
        """

        min_inc_counts = filters["included_counts"]
        min_exc_counts = filters["excluded_counts"]
        min_psi = filters["psi"]
        min_rpkm = filters["rpkm"]

        keep_splice_types = []

        # Find which splice types are marked to be included
        for name, info in filters["splice_type"].items():
            if info["enabled_var"] == 1:
                keep_splice_types.append(name)

        # Filter
        filtered_dataset = dataset.loc[
            (dataset["included_counts"] >= min_inc_counts) &
            (dataset["excluded_counts"] >= min_exc_counts) &
            (dataset["psi"] >= min_psi) &
            (dataset["rpkm"] >= min_rpkm)
        ]

        filtered_dataset = filtered_dataset.loc[filtered_dataset["splice_type"].isin(keep_splice_types)]

        # Check if there are remaining entries after filtering
        if filtered_dataset.empty:
            return None

        return filtered_dataset

    def save_filters(self, filters, filepath):
        """
        Saves the dataset filters, given in param filters, to disk, at location given by param filepath.
        Returns True if everything is OK, otherwise False.
        """

        try:
            with open(filepath, "w") as outfile:
                json.dump(filters, outfile)
                return True
        except IOError as e:
            print "Error when writing filters to file: %s", e.message
            return False

    def read_filters(self, filepath):
        """
        Read and return dataset filters from filepath. Returns a dictionary if everything works, otherwise None.
        """

        try:
            with open(filepath, "r") as infile:
                filters = json.load(infile)
                return filters
        except IOError as e:
            print "Error when reading filters from file: %s", e.message
            return None




