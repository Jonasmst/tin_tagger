import subprocess
import pandas as pd
import os
import json
import random


class TINDataProcessor(object):

    def __init__(self, tag_no_tag):
        self.tag_no_tag = tag_no_tag

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

        # Trim .0 from exon names
        df.exon1 = df.exon1.str.replace("\.0$", "")
        df.exon2 = df.exon2.str.replace("\.0$", "")

        # TODO: Don't do this; return full dataset
        df = df.loc[df.splice_type == "ES"]

        # Count occurrences if not already done
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
            #print "Found tag:", tag
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
            #print "Set tag:", new_tag
        except IndexError as e:
            print "ERROR: Can't find row index for sample %s, as_id %d. Message:\n%s" % (sample_name, as_id, e.message)

    def get_row_data(self, current_row_index, dataset, sample_names, bam_paths, testing):
        """
        Returns formatted data for the current row index.
        The format is as follows:

        {
            "splice_type": splice_type,
            "gene_symbol": gene_symbol,
            "sample_of_interest": sample_name,
            "location": coordinates,
            "exons": exons,
            "strand": strand,
            "as_id": as_id,
            "exon_psi": str(psi),
            "max_gene_rpkm": float,
            "samples":
                {
                    sample_name:
                    {
                        "gene_rpkm": float,
                        "event_tag": int,
                        "is_reported": bool,
                        "exons":
                        {
                            exon_name:
                            {
                                "coverage": int,
                                "max_coverage": int,
                                "psi": float
                            }
                        }
                    }
                }
        }
        """
        # TODO: Find included_counts and excluded_counts for other exons than the one in question

        # Get information
        splice_type = dataset.iloc[current_row_index]["splice_type"]
        sample_name = dataset.iloc[current_row_index]["name"]
        as_id = dataset.iloc[current_row_index]["as_id"]
        psi = dataset.iloc[current_row_index]["psi"]
        gene_symbol = dataset.iloc[current_row_index]["symbol"]
        strand = dataset.iloc[current_row_index]["strand"]
        exons = dataset.iloc[current_row_index]["exons"]
        chrom = dataset.iloc[current_row_index]["chr"]
        splice_start = dataset.iloc[current_row_index]["first_exon_start"]
        splice_stop = dataset.iloc[current_row_index]["last_exon_stop"]
        prev_exon_start = dataset.iloc[current_row_index]["prev_exon_start"]  # NaN for AT/AP
        prev_exon_stop = dataset.iloc[current_row_index]["prev_exon_stop"]  # NaN for AT/AP
        next_exon_start = dataset.iloc[current_row_index]["next_exon_start"]  # NaN for AT/AP
        next_exon_stop = dataset.iloc[current_row_index]["next_exon_stop"]  # NaN for AT/AP
        # Handle negative strand start- and stop- coordinates
        if strand == "-":
            splice_start = dataset.iloc[current_row_index]["last_exon_stop"]  # NaN for AT/AP
            splice_stop = dataset.iloc[current_row_index]["first_exon_start"]  # NaN for AT/AP
            prev_exon_start = dataset.iloc[current_row_index]["prev_exon_stop"]  # NaN for AT/AP
            prev_exon_stop = dataset.iloc[current_row_index]["prev_exon_start"]  # NaN for AT/AP
            next_exon_start = dataset.iloc[current_row_index]["next_exon_stop"]  # NaN for AT/AP
            next_exon_stop = dataset.iloc[current_row_index]["next_exon_start"]  # NaN for AT/AP
        included_counts = dataset.iloc[current_row_index]["included_counts"]
        excluded_counts = dataset.iloc[current_row_index]["excluded_counts"]
        prev_exon_name = dataset.iloc[current_row_index]["exon1"]  # NaN for AT/AP
        next_exon_name = dataset.iloc[current_row_index]["exon2"]  # NaN for AT/AP
        # prev_exon_id = dataset.iloc[current_row_index]["start_ex"]  # NaN for AT/AP
        # next_exon_id = dataset.iloc[current_row_index]["end_ex"]  # NaN for AT/AP

        # Create coordinates from chr, start and stop
        coordinates = str(chrom) + ":" + str(int(splice_start)) + "-" + str(int(splice_stop))

        # General row data
        row_data = {
            "splice_type": splice_type,
            "gene_symbol": gene_symbol,
            "sample_of_interest": sample_name,
            "location": coordinates,
            "exons": exons,
            "strand": strand,
            "as_id": as_id,
            "exon_psi": psi,
            "included_counts": included_counts,
            "excluded_counts": excluded_counts,
        }

        # Keep track of exon coverages and gene RPKMs
        all_gene_rpkms = []
        upstream_exon_coverages = [-1]
        downstream_exon_coverages = [-1]
        exon_of_interest_coverages = [-1]

        ########################
        # Sample-specific data #
        ########################
        samples_data = {}
        for s_name in sample_names:
            # Add entry for sample in the container
            if s_name not in samples_data.keys():
                samples_data[s_name] = {}

            # Find if sample is reported by SpliceSeq or not
            is_reported = self.is_event_reported_in_sample(s_name, as_id, dataset)
            samples_data[s_name]["is_reported"] = is_reported
            # Default to sample not being tagged
            sample_tag = self.tag_no_tag
            # Default to RPKM being 0 (in case it's not reported)
            gene_rpkm = self.get_gene_rpkm_by_sample_name(s_name, gene_symbol, dataset)

            if is_reported:
                sample_tag = self.get_sample_tag_by_as_id(s_name, as_id, dataset)

            samples_data[s_name]["gene_rpkm"] = gene_rpkm
            all_gene_rpkms.append(gene_rpkm)
            samples_data[s_name]["event_tag"] = sample_tag

            #########################
            # Find exon information #
            #########################
            sample_exons = {}
            # TODO: Also handle AT and AP events
            if splice_type in ["ES", "ME", "AD", "AA", "RI"]:
                # Handle upstream exon
                upstream_exon_coords = str(chrom) + ":" + str(int(prev_exon_start)) + "-" + str(int(prev_exon_stop))
                upstream_exon = {
                    "exon_name": prev_exon_name,
                    "coverage": self.get_coverage_by_coordinates(upstream_exon_coords, bam_paths[s_name], testing)
                }
                sample_exons["upstream_exon"] = upstream_exon

                # Handle downstream exon
                downstream_exon_coords = str(chrom) + ":" + str(int(next_exon_start)) + "-" + str(int(next_exon_stop))
                downstream_exon = {
                    "exon_name": next_exon_name,
                    "coverage": self.get_coverage_by_coordinates(downstream_exon_coords, bam_paths[s_name], testing)
                }
                sample_exons["downstream_exon"] = downstream_exon

                # Keep track of coverage values
                upstream_exon_coverages.append(upstream_exon["coverage"])
                downstream_exon_coverages.append(downstream_exon["coverage"])

            # Handle the main exon
            main_exon_coords = str(chrom) + ":" + str(int(splice_start)) + "-" + str(int(splice_stop))
            exon_of_interest = {
                "exon_name": exons,
                "coverage": self.get_coverage_by_coordinates(main_exon_coords, bam_paths[s_name], testing),
                "psi": psi,
                "included_counts": included_counts,
                "excluded_counts": excluded_counts
            }
            sample_exons["exon_of_interest"] = exon_of_interest

            # Keep track of exon coverage values
            exon_of_interest_coverages.append(exon_of_interest["coverage"])

            # Add exons data to this sample
            samples_data[s_name]["exons"] = sample_exons

        row_data["samples"] = samples_data
        row_data["max_gene_rpkm"] = max(all_gene_rpkms)
        row_data["max_upstream_exon_coverage"] = max(upstream_exon_coverages)
        row_data["max_downstream_exon_coverage"] = max(downstream_exon_coverages)
        row_data["max_exon_of_interest_coverage"] = max(exon_of_interest_coverages)

        return row_data









