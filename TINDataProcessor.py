import subprocess
import pandas as pd
import os
import json
import random


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

        # TODO: Count occurrences if not already done
        if "occurrences" not in list(df.columns):
            df["occurrences"] = df.groupby("as_id")["name"].transform(len)

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

    def is_event_reported_in_sample(self, sample_name, as_id, dataset):
        """
        Returns True if an event is reported for a given sample, otherwise returns False
        """
        event_sample_names = list(dataset.loc[dataset["as_id"] == as_id]["name"])
        return sample_name in event_sample_names

    def get_sample_tag_by_as_id(self, sample_name, as_id, dataset):
        """
        Returns the event tag for a given sample and as_id
        """
        row = dataset.loc[(dataset["as_id"] == as_id) & (dataset["name"] == sample_name)]
        sample_tag = row["event_tag"].iloc[0]
        return sample_tag

    def get_gene_rpkm_by_sample_name(self, sample_name, gene_symbol, dataset):
        """
        Returns the gene RPKM value for a given gene in a given sample.
        """
        try:
            gene_rpkm = dataset.loc[(dataset["symbol"] == gene_symbol) & (dataset["name"] == sample_name)]["rpkm"].iloc[0]
            return gene_rpkm
        except IndexError:
            # Data not found for this sample/gene
            return 0.0

    def get_coverage_by_coordinates(self, coordinates, bam_file_path, testing):
        """
        Runs samtools depth -r on a BAM-file to find the average number of reads covering a region. Returns the
        average coverage. The param 'testing' indicates whether or not the program is run for testing purposes
        (e.g. outside of the TSD firewall where BAM-files are not available) or not. For usual operation, 'testing' is
        false.
        """
        if testing:
            return random.randint(100, 2000)
        else:
            bam_file_path
            command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++;}END{if (cnt>0){ print sum/cnt } else print 0}'" % (coordinates, bam_file_path)
            samtools_output = subprocess.check_output(command, shell=True)

            region_coverage = -1.0
            try:
                region_coverage = float(samtools_output)
            except ValueError:
                print "ERROR: Samtools output can't be converted to float:"
                print samtools_output

            return region_coverage

    def get_bam_file_paths(self):
        """
        Returns a dictionary where keys are sample names and values are the path to that sample's BAM-file.
        """
        # TODO: Read this from a file or something
        # TODO: Should the DataProcessor take care of these paths, rather than the main app? Affects coverage-func.
        bam_directory = "/tsd/p19/data/durable/Projects/CRC-RNA-Seq/hisat2/transcriptome"
        bam_paths = {}

        for x in range(1, 11):
            bam_paths["sample%s" % str(x)] = os.path.join(bam_directory, "%s.sorted.bam" % str(x))

        return bam_paths

    def get_tag_by_sample_name_and_as_id(self, sample_name, as_id, dataset):
        """
        Returns the tag for this as_id for the given sample.
        """
        try:
            tag = dataset.loc[(dataset["as_id"] == as_id) & (dataset["name"] == sample_name)]["event_tag"].iloc[0]
            print "Found tag:", tag
            return tag
        except IndexError as e:
            print "ERROR, can't find tag for sample %s, as_id %d. Message:\n%s" % (sample_name, as_id, e.message)

    def set_tag_by_sample_name_and_as_id(self, new_tag, sample_name, as_id, dataset):
        """
        Sets a given tag for this sample and this as_id
        """

        # Get index of the row in question
        try:
            row_index = dataset.loc[(dataset["as_id"] == as_id) & (dataset["name"] == sample_name)].index.tolist()[0]
            # Assign new tag to this row
            dataset.set_value(row_index, "event_tag", new_tag)
            print "Set tag:", new_tag
        except IndexError as e:
            print "ERROR: Can't find row index for sample %s, as_id %d. Message:\n%s" % (sample_name, as_id, e.message)










