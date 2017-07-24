import subprocess
import pandas as pd
import os
import json
import random
import MySQLdb as mysql

# TODO: Check that label "event_tag" is present in the dataset. If not, create one and set all to default (-1)
# TODO: Deprecate samtools for now
# TODO: Max 2 decimals on exon coverage values

TAG_INTERESTING = 0
TAG_NOT_INTERESTING = 1
TAG_UNCERTAIN = 2
TAG_NO_TAG = -1

class TINDataProcessor(object):

    def __init__(self, tag_no_tag):
        self.testing = True
        self.tag_no_tag = tag_no_tag
        if not self.testing:
            self.db = mysql.connect("localhost", "crc_spliceseq", "TODO:insert_password_here", "crc_spliceseq")

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

        # Add event_tag column if not present
        if "event_tag" not in list(df.columns):
            df["event_tag"] = TAG_NO_TAG  # Default to no tag

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

    def get_main_exon_psi_by_asid(self, sample_names, as_id):
        """
        Queries the SpliceSeq database for PSI values for a given as_id in all samples.
        Returns a dictionary on the form:
        {
            "<sample name>": {
                "sample_id": <sample ID>,
                "psi": <PSI>,
                "included_counts": <included counts>,
                "excluded_counts": <excluded_counts>
            },
        }
        """
        print "Finding main exon PSI, included counts, excluded counts."
        query = """
        SELECT
            a.as_id,
            a.sample_id,
            sample.name,
            a.psi,
            a.included_counts,
            a.excluded_counts
        FROM as_counts AS a
        INNER JOIN sample
            ON sample.sample_id=a.sample_id
        WHERE
            sample.name IN(%s)
            AND
             a.as_id=%d
        """ % (",".join('"' + s + '"' for s in sample_names), as_id)

        if not self.testing:
            df = pd.read_sql_query(query, self.db)
        else:
            # If test-running, create dummy-data
            df = pd.DataFrame({
                "as_id": [71] * 10,
                "sample_id": [3, 4, 7, 8, 9, 10, 11, 12, 13, 14],
                "name": ["sample%s" % str(x) for x in range(1, 11)],
                "psi": [0.044, 0.334, 0.026, 0.022, 0, 0, 0, 0, 0, 0],
                "included_counts": [6, 6, 4, 4, 0, 0, 0, 0, 0, 0],
                "excluded_counts": [66, 81, 79, 83, 75, 39, 66, 61, 67, 75]
            })

        # Create a formatted dictionary for the results
        results_dict = {}
        for sample_name in sample_names:
            results_dict[sample_name] = {}
            sample_df = df.loc[df["name"] == sample_name]

            # Add sample-specific data to results
            results_dict[sample_name]["sample_id"] = sample_df["sample_id"].iloc[0]
            results_dict[sample_name]["psi"] = sample_df["psi"].iloc[0]
            results_dict[sample_name]["included_counts"] = sample_df["included_counts"].iloc[0]
            results_dict[sample_name]["excluded_counts"] = sample_df["excluded_counts"].iloc[0]

        return results_dict

    def get_main_exon_rpkm_by_asid(self, sample_names, as_id):
        """
        Queries the SpliceSeq database for RPKM values for a given exon.
        Returns a dictionary on the form:
        {
            "<sample name>": {
                "combined_rpkm": <total rpkm value for all exons combined>,
                "max_combined_rpkm": <max combined_rpkm across all samples>,
                "exons":
                    <exon_id>: {
                        "rpkm": <rpkm value>,
                        "max_rpkm": <max rpkm value>,
                        "tot_reads": <total reads value>
                    },
            },
        }
        """

        # First, find how many (partial) exons are in this splicing events
        if not self.testing:
            # TODO: Handle this being empty
            num_exons_query = "SELECT * FROM as_ref_exon WHERE as_id=%d;" % as_id
            num_exons_result = pd.read_sql_query(num_exons_query, self.db)
        else:
            num_exons_result = pd.DataFrame({
                "as_id": [as_id] * 2,
                "exon_id": [210, 211]
            })

        affected_exon_ids = list(num_exons_result["exon_id"].unique())

        print "Finding main exon RPKM"
        query = """
        SELECT
            s.name,
            s.sample_id,
            are.as_id,
            are.exon_id,
            ec.tot_reads,
            ec.rpkm
        FROM sample AS s
        INNER JOIN exon_counts AS ec
            ON s.sample_id=ec.sample_id
        INNER JOIN as_ref_exon AS are
            ON ec.exon_id=are.exon_id
        WHERE
            are.as_id=%d
            AND
            s.name IN(%s)
        """ % (as_id, ", ".join('"' + s + '"' for s in sample_names))

        if not self.testing:
            df = pd.read_sql_query(query, self.db)
        else:
            df = pd.DataFrame({
                "name": ["sample1", "sample2", "sample3", "sample4", "sample5", "sample8", "sample9", "sample10"] + ["sample1", "sample2", "sample3", "sample4"],
                "sample_id": [3, 4, 7, 8, 9, 12, 13, 14, 3, 4, 7, 8],
                "as_id": [71] * 12,
                "exon_id": [210] * 8 + [211] * 4,
                "tot_reads": [13, 12, 7, 8, 8, 5, 7, 8, 1, 2, 4, 1],
                "rpkm": [1.09, 0.96, 0.63, 0.58, 0.61, 0.42, 0.60, 0.66, 0.09, 0.17, 0.39, 0.07]
            })

        # If resulting dataframe is empty, create an empty DF with the correct column names so the merge still completes.
        if df.empty:
            df = pd.DataFrame({
                "name": [],
                "sample_id": [],
                "as_id": [],
                "exon_id": [],
                "tot_reads": [],
                "rpkm": []
            })

        # Rename resulting dataframe
        final_dataframe = df

        # Detect and handle unreported samples
        # If there's not (number of samples * number of exons) rows, something's missing
        if len(df) < (len(affected_exon_ids) * len(sample_names)):
            # Create template DataFrame to compare to resulting dataframe, containing sample names and exon ids
            number_of_samples = len(sample_names)
            template_exon_ids = []
            for exon in affected_exon_ids:
                template_exon_ids = template_exon_ids + [exon] * len(sample_names)
            template_df = pd.DataFrame({
                "name": sample_names * 2,
                "exon_id": template_exon_ids
            })

            # Merge with resulting dataframe
            merged = template_df.merge(df, on=["name", "exon_id"], how="left")

            # Forward fill as_id column
            merged["as_id"].fillna(method="ffill", inplace=True)

            # Fill rpkm and tot reads to 0
            filled = merged.fillna(0, axis=1)

            # Update dtypes that have somehow become floats
            filled[["as_id", "sample_id", "tot_reads"]] = filled[["as_id", "sample_id", "tot_reads"]].astype(int)

            final_dataframe = filled

        # Create column for combined RPKM and max combined RPKM
        final_dataframe["combined_rpkm"] = final_dataframe.groupby(["name", "exon_id"])["rpkm"].transform(sum)
        final_dataframe["max_combined_rpkm"] = final_dataframe["combined_rpkm"].max()

        # Convert output to a dictionary
        results_dict = {}
        max_rpkms = {}
        for ex in affected_exon_ids:
            max_rpkms[ex] = final_dataframe.loc[final_dataframe["exon_id"] == ex]["rpkm"].max()

        for sample_name in sample_names:
            results_dict[sample_name] = {}

            # Get info about this sample
            sample_df = final_dataframe.loc[final_dataframe["name"] == sample_name]
            results_dict[sample_name]["combined_rpkm"] = sample_df["combined_rpkm"].iloc[0]
            results_dict[sample_name]["max_combined_rpkm"] = sample_df["max_combined_rpkm"].iloc[0]

            # Create dictionary for exons
            results_dict[sample_name]["exons"] = {}

            # Get info about all exons in this sample
            for exon_id in affected_exon_ids:
                exon_df = sample_df.loc[sample_df["exon_id"] == exon_id]
                results_dict[sample_name]["exons"][exon_id] = {}
                results_dict[sample_name]["exons"][exon_id]["rpkm"] = exon_df["rpkm"].iloc[0]
                results_dict[sample_name]["exons"][exon_id]["max_rpkm"] = max_rpkms[exon_id]
                results_dict[sample_name]["exons"][exon_id]["tot_reads"] = exon_df["tot_reads"].iloc[0]

        print "Main exon RPKM query resulting dictionary:"
        print json.dumps(results_dict, indent=4)

        return results_dict

    def get_flanking_exons_rpkm_by_exon_ids(self, sample_names, prev_exon_id, next_exon_id):
        """
        Queries the SpliceSeq database for RPKM values for the flanking exons.
        Returns a dictionary of rpkm and tot_reads values for each exon in each sample:
        {
            "<sample name>": {
                "<exon id>": {
                    "rpkm": <rpkm_value>,
                    "tot_reads": <tot_reads_value>
                }
            }
        }
        """

        print "Finding flanking exons RPKM"
        query = """
        SELECT
            s.name,
            e.sample_id,
            e.exon_id,
            e.tot_reads,
            e.rpkm
        FROM exon_counts AS e
        INNER JOIN sample AS s
            ON s.sample_id=e.sample_id
        WHERE
            s.name IN(%s)
            AND
            e.exon_id IN(%s)
        """ % (", ".join("'" + s + "'" for s in sample_names), ", ".join(str(x) for x in [prev_exon_id, next_exon_id]))

        if not self.testing:
            df = pd.read_sql_query(query, self.db)
        else:
            df = pd.DataFrame({
                "name": ["sample%s" % str(x) for x in range(1, 11)] * 2,
                "sample_id": [3, 4, 7, 8, 9, 10, 11, 12, 13, 14] * 2,
                "exon_id": [prev_exon_id] * 10 + [next_exon_id] * 10,
                "tot_reads": [118, 135, 115, 125, 133, 58, 93, 93, 106, 136, 190, 200, 181, 190, 229, 134, 133, 154, 166, 188],
                "rpkm": [13.3, 14.5, 14.0, 12.3, 13.6, 5.4, 12.2, 10.5, 12.3, 15.2, 16.9, 17.1, 17.5, 14.8, 18.5, 10.0, 13.8, 13.7, 15.3, 16.6]
            })

        # If resulting dataframe is empty, create an empty DF with the correct column names so the merge still completes.
        if df.empty:
            df = pd.DataFrame({
                "name": [],
                "exon_id": [],
                "sample_id": [],
                "tot_reads": [],
                "rpkm": []
            })

        # Rename resulting dataframe
        final_dataframe = df

        # Detect and handle unreported samples
        # If there's not (number of samples * number of exons) rows, something's missing
        if len(df) < (2 * len(sample_names)):
            # Create template DataFrame to compare to resulting dataframe, containing sample names and exon ids
            number_of_samples = len(sample_names)
            template_df = pd.DataFrame({
                "name": sample_names * 2,
                "exon_id": [prev_exon_id] * number_of_samples + [next_exon_id] * number_of_samples
            })

            # Merge template DF with result DF to produce NaNs for RPKM and tot_reads in unreported samples
            merged_dataframe = template_df.merge(df, how="left", on=["rpkm", "tot_reads"])

            # Fill NaNs in rpkm and tot_reads columns with 0s
            final_dataframe = merged_dataframe.fillna(0, axis=1)

        # Convert output to a dictionary
        results_dict = {}
        # Include max RPKM for both exons
        prev_exon_max_rpkm = final_dataframe.loc[final_dataframe["exon_id"] == prev_exon_id]["rpkm"].max()
        next_exon_max_rpkm = final_dataframe.loc[final_dataframe["exon_id"] == next_exon_id]["rpkm"].max()
        max_rpkms = {
            prev_exon_id: prev_exon_max_rpkm,
            next_exon_id: next_exon_max_rpkm
        }
        for sample_name in sample_names:
            results_dict[sample_name] = {}

            # Get info about this sample
            sample_df = final_dataframe.loc[final_dataframe["name"] == sample_name]

            # Get info about both exons in this sample
            for exon_id in [prev_exon_id, next_exon_id]:
                exon_df = sample_df.loc[sample_df["exon_id"] == exon_id]
                results_dict[sample_name][exon_id] = {}
                results_dict[sample_name][exon_id]["rpkm"] = exon_df["rpkm"].iloc[0]
                results_dict[sample_name][exon_id]["tot_reads"] = exon_df["tot_reads"].iloc[0]
                results_dict[sample_name][exon_id]["max_rpkm"] = max_rpkms[exon_id]

        # print "Flanking exons RPKM/tot_reads dict:"
        # print json.dumps(results_dict, indent=4)

        return results_dict

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
        prev_exon_id = dataset.iloc[current_row_index]["start_ex"]  # NaN for AT/AP
        next_exon_id = dataset.iloc[current_row_index]["end_ex"]  # NaN for AT/AP

        print "AS_ID:", as_id
        print "Exon name:", exons
        print "Prev exon name:", prev_exon_name
        print "Prev exon ID:", int(prev_exon_id)
        print "Next exon name:", next_exon_name
        print "Next exon ID:", int(next_exon_id)
        print "-------------------------------"

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
            "prev_exon_id": int(prev_exon_id),
            "next_exon_id": int(next_exon_id)
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










