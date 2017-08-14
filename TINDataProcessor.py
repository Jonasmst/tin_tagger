import subprocess
import pandas as pd
import os
import json
import random
import MySQLdb as mysql
import numpy as np
from TINLearner import TINLearner

# TODO: Max 2 decimals on exon coverage values
# TODO: Test the relative differences between exons when using SAMtools and DB RPKM.
# - Are the relative differences the same? Otherwise, which one should we use?

TAG_NO_TAG = -1
TAG_INTERESTING = 0
TAG_NOT_INTERESTING = 1
TAG_UNCERTAIN = 2

# Fix pandas print width
pd.set_option("display.width", 250)
pd.set_option("max.columns", 100)


class TINDataProcessor(object):

    def __init__(self, tin_tagger, tag_no_tag):
        self.testing = True
        self.tin_tagger = tin_tagger
        self.tag_no_tag = tag_no_tag
        if not self.testing:
            print "DataProcessor: Connecting to DB"
            self.tin_tagger.set_statusbar_text("Hello from the DataProcessor!")
            self.db = mysql.connect("localhost", "crc_spliceseq", "TODO:insert_password_here", "crc_spliceseq")

        self.tin_learner = TINLearner(self)
        self.enable_decision_tree = True

        self.samtools_enabled = False  # Flag to determine whether or not to use SAMtools. For testing.

    def load_dataset(self, filepath, processQueue):
        """
        Reads a dataset and returns it as a pandas dataframe
        """
        print "DataProcessor: load_dataset()"

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
            "prev_exon_name": "object",
            "next_exon_name": "object",
            "prev_exon_chr_start": "int",
            "next_exon_chr_start": "int",
            "event_tag": "int",
            "avg_rpkm": "float",
            "max_avg_rpkm": "float",
            "prev_exon_rpkm": "float",
            "prev_exon_max_rpkm": "float",
            "next_exon_rpkm": "float",
            "next_exon_max_rpkm": "float",
            "max_gene_rpkm": "float",
            "prev_exon_tot_reads": "int",
            "next_exon_tot_reads": "int"
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
        #df.exon1 = df.exon1.str.replace("\.0$", "")
        #df.exon2 = df.exon2.str.replace("\.0$", "")

        # Count occurrences if not already done
        if "occurrences" not in list(df.columns):
            df["occurrences"] = df.groupby("as_id")["name"].transform(len)

        # Add event_tag column if not present
        if "event_tag" not in list(df.columns):
            df["event_tag"] = TAG_NO_TAG  # Default to no tag

        # Calc max RPKM for prev exon
        df["prev_exon_max_rpkm"] = df.groupby("as_id")["prev_exon_rpkm"].transform(max)

        # Calc max RPKM for next exon
        df["next_exon_max_rpkm"] = df.groupby("as_id")["next_exon_rpkm"].transform(max)

        # Calc max average RPKM for main exon
        df["max_avg_rpkm"] = df.groupby("as_id")["avg_rpkm"].transform(max)

        # Calc max gene RPKM
        df["max_gene_rpkm"] = df.groupby("as_id")["rpkm"].transform(max)

        processQueue.put(df)

    def filter_dataset(self, dataset, filters):
        """
        Takes a Pandas Dataframe in param dataset and filter it based on criteria given in param filters.
        Returns a Pandas DataFrame if there are datapoints in the remaining filtered dataset, otherwise returns None.
        """

        # Work on a copy of the dataset
        df = dataset.copy()

        # Filter on ints
        int_fields = ["included_counts", "excluded_counts", "tot_reads", "occurrences"]
        for i in int_fields:
            df = df.loc[df[i] >= int(filters[i][1].get())]

        # Filter float fields
        float_fields = [
            "psi",
            "rpkm",
            "avg_rpkm",
            "prev_exon_tot_reads",
            "prev_exon_rpkm",
            "next_exon_tot_reads",
            "next_exon_rpkm",
            "avg_tot_reads",
            "prev_exon_max_rpkm",
            "next_exon_max_rpkm",
            "max_avg_rpkm",
            "max_gene_rpkm",
            "max_psi",
            "percent_of_max_psi",
            "percent_of_max_rpkm",
            "main_rpkm_to_upstream_rpkm_ratio",
            "main_rpkm_to_downstream_rpkm_ratio",
            "sum_psi_all_samples",
            "sum_psi_other_samples",
            "mean_psi_other_samples",
            "psi_diff_from_mean_other_samples",
            "sum_rpkm_all_samples",
            "sum_rpkm_other_samples",
            "mean_rpkm_other_samples",
            "rpkm_percentage_of_mean_other_samples"
        ]
        for f in float_fields:
            df = df.loc[df[f] >= float(filters[f][1].get())]

        # Filter on splice types
        include_types = []  # List of splice types to include
        for st, st_fields in filters["splice_type"].items():
            if st_fields[0]:
                include_types.append(st)

        df = df.loc[df["splice_type"].isin(include_types)]

        # Filter on tags
        include_tags = []  # List of tags to include
        if filters["event_tag"]["interesting"][0]:
            include_tags.append(TAG_INTERESTING)
        if filters["event_tag"]["not_interesting"][0]:
            include_tags.append(TAG_NOT_INTERESTING)
        if filters["event_tag"]["uncertain"][0]:
            include_tags.append(TAG_UNCERTAIN)
        if filters["event_tag"]["no_tag"][0]:
            include_tags.append(TAG_NO_TAG)

        df = df.loc[df["event_tag"].isin(include_tags)]

        # Return filtered dataset
        return df

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

    def get_next_untagged_asid(self, as_id, dataset):

        # Find all events that are not yet tagged and get their as_ids
        untagged = dataset.loc[dataset["event_tag"] == TAG_NO_TAG]
        untagged_asids = sorted(list(untagged.as_id.unique()))

        # Return the first untagged as_id, if any
        if len(untagged_asids) > 0:
            return untagged_asids[0]

        # Otherwise just return the same as_id that was passed in
        print "ERROR: No untagged events found"
        return as_id

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

    def get_row_data(self, as_id, dataset, sample_names, testing):
        """
        Returns formatted data for the current row index.
        The format is as follows:

        row_data = {
            "splice_type": <splice_type>,
            "gene_symbol": <gene_symbol>,
            "location": <coordinates>,
            "exons": <exons>,
            "strand": <strand>,
            "as_id": <as_id>,
            "max_gene_rpkm": <max_gene_rpkm>,
            "prev_exon_id": <prev_exon_id>,
            "next_exon_id": <next_exon_id>,
            "prev_exon_name": <prev_exon_name>,
            "next_exon_name": <next_exon_name>,
            "samples":
                {
                    <sample_name>:
                    {
                        "gene_rpkm": float,
                        "event_tag": int,
                        "is_reported": bool,
                    }
                }
        }
        """
        # Get information
        event_df = dataset.loc[dataset["as_id"] == as_id]
        splice_type = event_df["splice_type"].iloc[0]
        gene_symbol = event_df["symbol"].iloc[0]
        strand = event_df["strand"].iloc[0]
        exons = event_df["exons"].iloc[0]
        chrom = event_df["chr"].iloc[0]
        prev_exon_name = event_df["prev_exon_name"].iloc[0]
        next_exon_name = event_df["next_exon_name"].iloc[0]
        prev_exon_id = event_df["start_ex"].iloc[0]
        next_exon_id = event_df["end_ex"].iloc[0]

        # Try int-casting exon IDs
        try:
            prev_exon_id = int(prev_exon_id)
            next_exon_id = int(next_exon_id)
        except ValueError:
            # Don't do anything.
            pass

        # General row data
        row_data = {
            "splice_type": splice_type,
            "gene_symbol": gene_symbol,
            "exons": exons,
            "strand": strand,
            "as_id": as_id,
            "prev_exon_id": prev_exon_id,
            "next_exon_id": next_exon_id,
            "prev_exon_name": prev_exon_name,
            "next_exon_name": next_exon_name,
            "coords": event_df["coords"].iloc[0]
        }

        # TEST: Get predicted tags from decision tree
        if self.enable_decision_tree:
            decision_tree_tags = self.tin_learner.predict_tag_decision_tree(event_df)
            # decision_tree_tags is either a dataframe or a boolean. Pandas throws a ValueError when checking the truth
            # value of a dataframe, so we have to surround it with a try/catch
            try:
                if not decision_tree_tags:
                    print "ERROR: Decision tree has not been fitted, yet attempted to predict value"
            except ValueError:
                # Result is a dataframe
                event_df.is_copy = False  # Suppress the annoying SettingWithCopyWarning
                event_df.loc[:, "decision_tree_tag"] = decision_tree_tags
                # TODO: At some point include the predicted tags in the TINTagger's original dataset
        # END TEST

        ########################
        # Sample-specific data #
        ########################
        samples_data = {}
        for s_name in sample_names:

            # Add entry for sample in the container
            if s_name not in samples_data.keys():
                samples_data[s_name] = {}

            try:
                # Get row for this sample
                sample_row = dataset.loc[(dataset["as_id"] == as_id) & (dataset["name"] == s_name)].iloc[0]
                # Find if sample is reported by SpliceSeq or not
                samples_data[s_name]["is_reported"] = self.is_event_reported_in_sample(s_name, as_id, dataset)
                samples_data[s_name]["gene_rpkm"] = sample_row["rpkm"]
                samples_data[s_name]["max_gene_rpkm"] = sample_row["max_gene_rpkm"]
                samples_data[s_name]["event_tag"] = sample_row["event_tag"]
                samples_data[s_name]["prev_exon_rpkm"] = sample_row["prev_exon_rpkm"]
                samples_data[s_name]["prev_exon_max_rpkm"] = sample_row["prev_exon_max_rpkm"]
                samples_data[s_name]["next_exon_rpkm"] = sample_row["next_exon_rpkm"]
                samples_data[s_name]["next_exon_max_rpkm"] = sample_row["next_exon_max_rpkm"]
                samples_data[s_name]["avg_rpkm"] = sample_row["avg_rpkm"]
                samples_data[s_name]["max_avg_rpkm"] = sample_row["max_avg_rpkm"]
                samples_data[s_name]["psi"] = sample_row["psi"]
                samples_data[s_name]["included_counts"] = sample_row["included_counts"]
                samples_data[s_name]["excluded_counts"] = sample_row["excluded_counts"]
                samples_data[s_name]["max_gene_rpkm"] = sample_row["max_gene_rpkm"]
                samples_data[s_name]["decision_tree_prediction"] = sample_row["decision_tree_tag"]
            except IndexError:
                # Sample not present, fill with "blanks"
                samples_data[s_name]["is_reported"] = False
                samples_data[s_name]["gene_rpkm"] = 0
                samples_data[s_name]["max_gene_rpkm"] = 0
                samples_data[s_name]["event_tag"] = -1
                samples_data[s_name]["prev_exon_rpkm"] = 0
                samples_data[s_name]["prev_exon_max_rpkm"] = 0
                samples_data[s_name]["next_exon_rpkm"] = 0
                samples_data[s_name]["next_exon_max_rpkm"] = 0
                samples_data[s_name]["avg_rpkm"] = 0
                samples_data[s_name]["max_avg_rpkm"] = 0
                samples_data[s_name]["psi"] = 0
                samples_data[s_name]["included_counts"] = 0
                samples_data[s_name]["excluded_counts"] = 0
                samples_data[s_name]["max_gene_rpkm"] = 0

        row_data["samples"] = samples_data

        return row_data

    def get_rpkm_for_mutually_exclusive_exons(self, sample_names, as_id):
        """
        Finds rpkm-values for exons affected by a mutually exclusive exons (ME) event. If multiple exons present on one
        or both of the mutually exclusive exons, returns the combined RPKM of all (partial) exons. E.g., for the event
        "4\5.1:5.2", returns the RPKM for exon 4 and the combined RPKM for exon 5.1 + 5.2. Returned dictionary takes
        the form:
        {
           <sample_name>: {
                "first_exon": {
                    "joined_exon_name": <Name as reported by SpliceSeq, e.g. "5.1:5.2", or "4">,
                    "combined_rpkm": <Combined RPKM value of all exons comprising the first part of the ME event>
                },
                "second_exon": {
                    "joined_exon_name": <Same as for 'first_exon'>,
                    "combined_rpkm": <Same as for 'first_exon'>
                }
           },
           ...
        }
        """
        print "DataProcessor: get_rpkm_for_mutually_exclusive_exons"
        # Query the SpliceSeq DB
        query = """
        SELECT
            s.name,
            s.sample_id,
            ar.exons,
            are.as_id,
            are.exon_id,
            e.exon_name,
            ec.tot_reads,
            ec.rpkm
        FROM sample AS s
        INNER JOIN exon_counts AS ec
            ON s.sample_id=ec.sample_id
        INNER JOIN as_ref_exon AS are
            ON ec.exon_id=are.exon_id
        INNER JOIN exon AS e
            ON e.exon_id=are.exon_id
        INNER JOIN as_ref AS ar
            ON ar.as_id=are.as_id
        WHERE
            are.as_id=%d
            AND
            s.name IN(%s)
        """ % (as_id, ", ".join('"' + s + '"' for s in sample_names))

        # TODO: Handle some samples not being reported for
        if self.testing:
            sample_names = ["sample%s" % str(x) for x in range(1, 11)]
            df = pd.DataFrame({
                "name": sample_names * 3,
                "sample_id": [3, 4, 7, 8, 9, 10, 11, 12, 13, 14] * 3,
                "exons": ["5|6.1:6.2"] * 30,
                "as_id": [46834] * 30,
                "exon_id": [169331] * 10 + [169332] * 10 + [169333] * 10,
                "exon_name": ["5"] * 10 + ["6.1"] * 10 + ["6.2"] * 10,
                "tot_reads": random.sample(range(0, 100), 30),
                "rpkm": random.sample(range(200, 1000), 30)
            })
        else:
            # Read results into Pandas DataFrame
            df = pd.read_sql_query(query, self.db)

        # In ME events, there are always two main exons, first and second
        # TODO: Handle possibility of result being empty DF
        affected_exons = df["exons"].iloc[0]
        first_main_exons, second_main_exons = affected_exons.split("|")  # Both may contain partials, e.g.: "5.1:5.2"

        # Get all exons in first and second main exons (if multiple partial exons)
        all_first_main_exons = first_main_exons.split(":")  # Split into a list, e.g. ["5.1", "5.2"]
        all_second_main_exons = second_main_exons.split(":")  # Split into a list, e.g. ["5.1", "5.2"]

        # Grab RPKM values for the different exons for each sample
        return_data = {}
        for sample_name in sample_names:
            # TODO: Handle the possibility of this not being reported for every sample
            # Allocate an entry for this sample in the return data
            return_data[sample_name] = {}
            if sample_name not in list(df["name"].unique()):
                # Sample not reported
                return_data[sample_name]["is_reported"] = False
                # Move on to next sample
                continue
            return_data[sample_name]["is_reported"] = True
            # Get data for this sample
            sample_data = df.loc[df["name"] == sample_name]
            # Get RPKM values for both main exons
            first_exons_data = sample_data.loc[sample_data["exon_name"].isin(all_first_main_exons)]
            return_data[sample_name]["first_exon"] = {
                "joined_exon_name": ":".join(all_first_main_exons),
                "combined_rpkm": first_exons_data["rpkm"].sum()
            }
            # Do the same for the second main exon
            second_exon_data = sample_data.loc[sample_data["exon_name"].isin(all_second_main_exons)]
            return_data[sample_name]["second_exon"] = {
                "joined_exon_name": ":".join(all_second_main_exons),
                "combined_rpkm": second_exon_data["rpkm"].sum()
            }

        # Find max RPKM for first and second main exons
        all_first_exon_rpkms = []
        all_second_exon_rpkms = []
        for sample_name in sample_names:
            if return_data[sample_name]["is_reported"]:
                all_first_exon_rpkms.append(return_data[sample_name]["first_exon"]["combined_rpkm"])
                all_second_exon_rpkms.append(return_data[sample_name]["second_exon"]["combined_rpkm"])

        max_first_exon_rpkm = max(all_first_exon_rpkms)
        max_second_exon_rpkm = max(all_second_exon_rpkms)

        # Add max values to all samples' return data
        for sample_name in sample_names:
            if return_data[sample_name]["is_reported"]:
                return_data[sample_name]["first_exon"]["max_combined_rpkm"] = max_first_exon_rpkm
                return_data[sample_name]["second_exon"]["max_combined_rpkm"] = max_second_exon_rpkm

        # Finally, return the data
        return return_data

    def get_dataset_from_database(self, db_url, db_user, db_pass, db_name, queue):
        """
        Queries the database to retrieve information about main exon PSI/RPKM and flanking exons RPKM.
        """
        # Connect do database
        try:
            self.db = mysql.connect(db_url, db_user, db_pass, db_name)
        except mysql.OperationalError as e:
            print "ERROR: Unable to connect to database:"
            print e.message
            queue.put(False)
            return

        # Sample names
        sample_names = ["sample%s" % str(x) for x in range(1, 11)]

        # Find RPKM for flanking exons
        flanking_exons_query = """
        SELECT \
            sample.name, sample.sample_id, \
            ex1.exon_name AS prev_exon_name, ex1.chr_start AS prev_exon_chr_start, ex1.chr_stop AS prev_exon_chr_stop, \
            ex2.exon_name AS next_exon_name, ex2.chr_start AS next_exon_chr_start, ex2.chr_stop AS next_exon_chr_stop, \
            ar.as_id, ar.splice_type, ar.exons, \
            ec1.tot_reads AS prev_exon_tot_reads, ec1.rpkm AS prev_exon_rpkm, \
            ec2.tot_reads AS next_exon_tot_reads, ec2.rpkm AS next_exon_rpkm \
        FROM \
            sample \
            INNER JOIN as_counts ON as_counts.sample_id=sample.sample_id \
            INNER JOIN as_ref AS ar ON ar.as_id=as_counts.as_id \
            INNER JOIN exon AS ex1 ON ar.start_ex=ex1.exon_id \
            INNER JOIN exon AS ex2 ON ar.end_ex=ex2.exon_id \
            LEFT JOIN exon_counts AS ec1 ON (ar.start_ex=ec1.exon_id AND sample.sample_id=ec1.sample_id) \
            LEFT JOIN exon_counts AS ec2 ON (ar.end_ex=ec2.exon_id AND sample.sample_id=ec2.sample_id)
        WHERE \
            sample.name IN(%s)
        """ % (",".join('"' + s + '"' for s in sample_names))
        self.tin_tagger.set_statusbar_text("Querying for flanking exons RPKM")
        flanking_exons_df = pd.read_sql_query(flanking_exons_query, self.db)

        # Find PSI and included/excluded counts for main exon
        main_exon_query = """
        SELECT \
            sample.name, sample.sample_id, \
            ac.as_id, ac.psi, ac.included_counts, ac.excluded_counts, \
            ar.splice_type, ar.start_ex, ar.end_ex, ar.novel_splice, ar.exons, \
            g.graph_id, g.symbol, g.chr, g.strand, \
            gc.rpkm
        FROM
            sample \
            INNER JOIN as_counts AS ac ON ac.sample_id=sample.sample_id \
            INNER JOIN as_ref AS ar ON ar.as_id=ac.as_id \
            INNER JOIN graph AS g ON g.graph_id=ar.graph_id \
            INNER JOIN gene_counts AS gc ON gc.graph_id=ar.graph_id AND gc.sample_id=sample.sample_id \
        WHERE \
            sample.name IN(%s)
        """ % (",".join('"' + s + '"' for s in sample_names))

        print "\nQuerying for main exon"
        self.tin_tagger.set_statusbar_text("Querying for main exon PSI")
        main_exon_df = pd.read_sql_query(main_exon_query, self.db)

        # Merge flanking exons data and main exons data together into a single dataset
        merged_df = main_exon_df.merge(flanking_exons_df, on=["sample_id", "as_id", "name", "exons", "splice_type"], how="outer")

        # Find average RPKM for every exon in the main exon(s)
        main_exon_rpkm_query = """
        SELECT \
            sample.sample_id, sample.name, \
            as_counts.as_id, as_ref_exon.exon_id, \
            exon_counts.rpkm, exon_counts.tot_reads, \
            exon.exon_name, exon.chr_start, exon.chr_stop
        FROM sample \
            INNER JOIN as_counts ON as_counts.sample_id=sample.sample_id \
            INNER JOIN as_ref ON as_ref.as_id=as_counts.as_id \
            INNER JOIN as_ref_exon ON as_ref_exon.as_id=as_ref.as_id \
            LEFT JOIN exon_counts ON exon_counts.sample_id=sample.sample_id AND exon_counts.exon_id=as_ref_exon.exon_id \
            LEFT JOIN exon ON exon.exon_id=as_ref_exon.exon_id AND exon.graph_id=as_ref.graph_id
        WHERE \
            sample.name IN(%s);
        """ % (",".join('"' + s + '"' for s in sample_names))
        self.tin_tagger.set_statusbar_text("Querying for main exon RPKM")
        main_exon_rpkm_df = pd.read_sql_query(main_exon_rpkm_query, self.db)
        # Replace NaNs with 0 (IGV-lookup shows that non-reported RPKMs is due to zero read count)
        main_exon_rpkm_df["rpkm"].fillna(0, inplace=True)
        # Calc average RPKM for all main exons in each sample
        main_exon_rpkm_df["avg_rpkm"] = main_exon_rpkm_df.groupby(["as_id", "name"])["rpkm"].transform("mean")
        # Replace NaNs of tot_reads with 0
        main_exon_rpkm_df["tot_reads"].fillna(0, inplace=True)
        # Calc average total reads for all main exons in each sample
        main_exon_rpkm_df["avg_tot_reads"] = main_exon_rpkm_df.groupby(["as_id", "name"])["rpkm"].transform("mean")

        # Remove duplicates of as_id and name (otherwise we'll have one identical entry per exon name in the main exon)
        self.tin_tagger.set_statusbar_text("Removing duplicates..")
        main_exon_rpkm_deduped = main_exon_rpkm_df.drop_duplicates(["as_id", "name"])[["sample_id", "name", "as_id", "avg_rpkm", "chr_start", "chr_stop", "tot_reads", "avg_tot_reads"]]

        self.tin_tagger.set_statusbar_text("Merging datasets..")
        # Merge together datasets to include average RPKM for main exon
        unprocessed_final = merged_df.merge(main_exon_rpkm_deduped, on=["name", "sample_id", "as_id"], how="inner")

        # Drop NaNs in prev_exon_tot_reads, prev_exon_rpkm, next_exon_tot_reads, next_exon_rpkm
        final_df = unprocessed_final.copy()
        final_df["prev_exon_tot_reads"].fillna(0, inplace=True)
        final_df["prev_exon_rpkm"].fillna(0, inplace=True)
        final_df["next_exon_tot_reads"].fillna(0, inplace=True)
        final_df["next_exon_rpkm"].fillna(0, inplace=True)
        final_df["next_exon_chr_start"].fillna(-1, inplace=True)
        final_df["next_exon_chr_stop"].fillna(-1, inplace=True)
        final_df["next_exon_name"].fillna("N/A", inplace=True)
        final_df["prev_exon_chr_start"].fillna(-1, inplace=True)
        final_df["prev_exon_chr_stop"].fillna(-1, inplace=True)
        final_df["prev_exon_name"].fillna("N/A", inplace=True)
        final_df["start_ex"].fillna(-1, inplace=True)
        final_df["end_ex"].fillna(-1, inplace=True)
        final_df["avg_rpkm"].fillna(0, inplace=True)
        final_df["rpkm"].fillna(0, inplace=True)
        final_df["tot_reads"].fillna(0, inplace=True)
        final_df["avg_tot_reads"].fillna(0, inplace=True)

        # Count occurrences if not already done
        if "occurrences" not in list(final_df.columns):
            final_df["occurrences"] = final_df.groupby("as_id")["name"].transform(len)

        # Add event_tag column if not present
        if "event_tag" not in list(final_df.columns):
            final_df["event_tag"] = TAG_NO_TAG  # Default to no tag

        self.tin_tagger.set_statusbar_text("Finding RPKM max. values and ratios.")

        # Calc max RPKM for prev exon
        final_df["prev_exon_max_rpkm"] = final_df.groupby("as_id")["prev_exon_rpkm"].transform(max)

        # Calc max RPKM for next exon
        final_df["next_exon_max_rpkm"] = final_df.groupby("as_id")["next_exon_rpkm"].transform(max)

        # Calc max average RPKM for main exon
        final_df["max_avg_rpkm"] = final_df.groupby("as_id")["avg_rpkm"].transform(max)

        # Calc max gene RPKM
        final_df["max_gene_rpkm"] = final_df.groupby("as_id")["rpkm"].transform(max)

        # TEST: Fill NaNs in max_gene_rpkm
        final_df["max_gene_rpkm"] = final_df.groupby("as_id")["max_gene_rpkm"].transform(lambda s: s.loc[s.first_valid_index()])
        # END TEST

        # Calc max PSI
        self.tin_tagger.set_statusbar_text("Finding PSI max. values and ratios")
        final_df["max_psi"] = final_df.groupby("as_id")["psi"].transform(max)

        # Create coords column
        final_df["coords"] = final_df["chr"].map(str) + ":" + final_df["chr_start"].map(str) + "-" + final_df["chr_stop"].map(str)

        # Datatypes
        datatypes = {
            "name": str,
            "as_id": int,
            "psi": float,
            "included_counts": int,
            "excluded_counts": int,
            "splice_type": str,
            "start_ex": int,
            "end_ex": int,
            "novel_splice": int,
            "exons": str,
            "graph_id": int,
            "symbol": str,
            "chr": str,
            "strand": str,
            "rpkm": float,
            "prev_exon_name": str,
            "prev_exon_chr_start": int,
            "prev_exon_chr_stop": int,
            "next_exon_name": str,
            "next_exon_chr_start": int,
            "next_exon_chr_stop": int,
            "prev_exon_tot_reads": int,
        }

        # Convert datatypes
        for column_name, column_type in datatypes.items():
            final_df[column_name] = final_df[column_name].astype(column_type)

        # Variables to be used for learning purposes # TODO: Division by zero errors
        print "Assigning variables for algo learning purposes"
        # Percent of max PSI
        final_df["percent_of_max_psi"] = final_df["psi"] / final_df["max_psi"]
        final_df["percent_of_max_psi"].fillna(0, inplace=True)
        # Percent of max RPKM
        final_df["percent_of_max_rpkm"] = final_df["avg_rpkm"] / final_df["max_avg_rpkm"]
        final_df["percent_of_max_rpkm"].fillna(0, inplace=True)
        # Main exon to upstream RPKM ratio
        final_df["main_rpkm_to_upstream_rpkm_ratio"] = final_df["avg_rpkm"] / final_df["prev_exon_rpkm"]
        final_df["main_rpkm_to_upstream_rpkm_ratio"].fillna(0, inplace=True)
        final_df["main_rpkm_to_upstream_rpkm_ratio"] = final_df["main_rpkm_to_upstream_rpkm_ratio"].replace(np.inf, 1)
        # Main exon to downstream RPKM ratio
        final_df["main_rpkm_to_downstream_rpkm_ratio"] = final_df["avg_rpkm"] / final_df["next_exon_rpkm"]
        final_df["main_rpkm_to_downstream_rpkm_ratio"].fillna(0, inplace=True)
        final_df["main_rpkm_to_downstream_rpkm_ratio"] = final_df["main_rpkm_to_downstream_rpkm_ratio"].replace(np.inf, 1)

        # PSI diff (in percentage) from mean PSI in other samples
        final_df["sum_psi_all_samples"] = final_df.groupby("as_id")["psi"].transform(sum)
        final_df["sum_psi_other_samples"] = final_df["sum_psi_all_samples"] - final_df["psi"]
        final_df["mean_psi_other_samples"] = final_df["sum_psi_other_samples"] / (final_df["occurrences"].astype(float) - 1)
        # For events only occurring in 1 sample, mean_psi_other_samples will be np.inf. We'll replace them with 0s and
        # just ignore occurrences < 4 or something when training algos.
        final_df["mean_psi_other_samples"].replace(np.inf, 0.00)
        final_df["psi_diff_from_mean_other_samples"] = final_df["psi"] - final_df["mean_psi_other_samples"]

        # RPKM diff (in percentage) from mean RPKM in other samples
        final_df["sum_rpkm_all_samples"] = final_df.groupby("as_id")["avg_rpkm"].transform(sum)
        final_df["sum_rpkm_other_samples"] = final_df["sum_rpkm_all_samples"] - final_df["avg_rpkm"]
        final_df["mean_rpkm_other_samples"] = final_df["sum_rpkm_other_samples"] / (final_df["occurrences"].astype(float) - 1)
        # For events only occurring in 1 sample, mean_rpkm_other_samples will be np.inf. We'll replace them with 0s and
        # just ignore occurrences < 4 or something when training algos.
        final_df["mean_rpkm_other_samples"].replace(np.inf, 0.00)
        final_df["rpkm_percentage_of_mean_other_samples"] = final_df["avg_rpkm"] / final_df["mean_rpkm_other_samples"]
        final_df["rpkm_percentage_of_mean_other_samples"].replace(np.inf, 0.00)
        final_df["rpkm_percentage_of_mean_other_samples"].fillna(0.00)

        # One-hot encode splice_type column
        self.tin_tagger.set_statusbar_text("One-hot encoding splice type column")
        splicetype_dummies = pd.get_dummies(final_df.splice_type, prefix="splicetype", drop_first=True)
        final_df = pd.concat([final_df, splicetype_dummies], axis=1)

        # Columns for ML assigned tags
        final_df["decision_tree_tag"] = TAG_NO_TAG
        final_df["random_forest_tag"] = TAG_NO_TAG
        final_df["neural_net_tag"] = TAG_NO_TAG

        self.tin_tagger.set_statusbar_text("Done fetching and preprocessing data.")

        # Finally, add dataframe to queue
        queue.put(final_df)






