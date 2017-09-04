from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.model_selection import train_test_split
from sklearn.exceptions import NotFittedError


class TINLearner(object):

    def __init__(self, tin_dataprocessor):
        self.data_processor = tin_dataprocessor
        self.tag_column = "event_tag"
        self.decision_tree = DecisionTreeClassifier()

        # Present columns to include in training
        self.training_columns = [
            "psi",
            "included_counts",
            "excluded_counts",
            "rpkm",
            "prev_exon_tot_reads",
            "next_exon_tot_reads",
            "prev_exon_rpkm",
            "next_exon_rpkm",
            "avg_rpkm",
            "occurrences",
            "prev_exon_max_rpkm",
            "next_exon_max_rpkm",
            "max_avg_rpkm",
            "max_gene_rpkm",
            "splicetype_AD",
            "splicetype_AP",
            "splicetype_AT",
            "splicetype_ES",
            "splicetype_ME",
            "splicetype_RI",
            "percent_of_max_psi",
            "percent_of_max_rpkm",
            "main_rpkm_to_upstream_rpkm_ratio",
            "main_rpkm_to_downstream_rpkm_ratio",
            "psi_diff_from_mean_other_samples",
            "rpkm_percentage_of_mean_other_samples",
            "rpkm_abs_diff_mean_other_samples"
        ]

    def prepare_dataset(self, dataset):
        """
        Prepares dataset for training. Ensures all appropriate columns are present
        in the dataset. Extract only events that have been tagged.
        """
        # TODO: Sanitize column names and dataset size

        """
        PARAMETERS:
        -----------
        - # Main exon RPKM (combined if multiple exons) --> avg_rpkm
        - # Upstream exon RPKM --> prev_exon_rpkm
        - # Downstream exon RPKM --> next_exon_rpkm
        - # Main exon RPKM diff from upstream exon RPKM --> main_rpkm_to_upstream_rpkm_ratio
        - # Main exon RPKM diff from downstream exon RPKM --> main_rpkm_to_downstream_rpkm_ratio
        - # Main exon RPKM percent of max RPKM --> percent_of_max_rpkm 
        - # Main exon PSI percent of max PSI --> percent_of_max_psi
        - # Main exon PSI diff from mean main exon PSI in other samples --> psi_diff_from_mean_other_samples
        - # Main exon diff from mean main exon RPKM in other samples --> rpkm_percentage_of_mean_other_samples 
        - Upstream exon RPKM diff from mean upstream RPKM in other samples
        - Downstream exon RPKM diff from mean downstream RPKM in other samples
        - # Main exon included counts
        - # Main exon excluded counts
        - Diff in main exon included counts from mean main exon included counts in other samples
        - Diff in main exon excluded counts from mean main exon excluded coutns in other samples
        - Splice type mapped to numeric
        - Is_reported mapped to numeric?
        """

        # Extract events that have been tagged
        tagged_events = dataset.loc[dataset[self.tag_column] > -1]  # -1 means they're not tagged. 0, 1, 2 are valid tags.
        tagged_events = tagged_events.loc[tagged_events["occurrences"] > 1]  # NaNs are present if event's only present in a single sample
        print "Tagged events:", len(tagged_events)
        minimum_tagged_events = 10
        if len(tagged_events) < minimum_tagged_events:
            return False

        return tagged_events

    def train_decision_tree(self, dataset):
        """
        Trains a decision tree on the provided dataset.
        """

        # Get a sanitized dataset. This returns False if there are too few event tags to train the tree
        dataset = self.prepare_dataset(dataset)

        try:
            if not dataset:
                print "Cannot train decision tree: Too few tagged events"
                return False
        except ValueError:
            # It's a dataframe, not a boolean, which means training succeeded. Continue operation
            pass

        # Split data into training set and validation set
        training_set, validation_set = train_test_split(dataset, test_size=0.2)

        # Fit tree
        """
        Tips on fitting a tree:
        -----------------------
        - Control the number of samples at the leaf nodes. Few samples often means overfitting, more samples means it
          can't memorize the data. Use min_samples_leaf=5 (for example) in the DT classifier.
        - Visualize the tree as you are training by using the export function. Use max_depth=3 as an initial tree depth
          to get a feel for how the tree is fitting to your data, and then increase the depth.
          
        """
        self.decision_tree.fit(training_set[self.training_columns], training_set[self.tag_column])

        # Test accuracy of tree
        score = self.decision_tree.score(validation_set[self.training_columns], validation_set[self.tag_column])

        print "Decision tree prediction accuracy: %.2f" % score

        print_tree = True
        if print_tree:
            from sklearn.externals.six import StringIO
            import pydot
            dot_data = StringIO()
            export_graphviz(
                self.decision_tree,
                out_file=dot_data,
                class_names = ["Deviating", "Non deviating"],
                feature_names=self.training_columns,
                filled=True,
                leaves_parallel=True,
                rounded=True,
                max_depth=3
            )
            graph = pydot.graph_from_dot_data(dot_data.getvalue())
            graph[0].write_pdf("decision_tree.pdf")
            print "Tree written to file decision_tree.pdf"

        # Training successful, return accuracy
        return score

    def predict_tag_decision_tree(self, event_df):
        """
        Takes an event DataFrame and predicts the value of event_tag column

        :return: False if tree has not yet been fitted, np.array of tags if it has been fitted.
        """

        # Only predict if the event occurs in more than one sample (otherwise a lot of the params are NaNs)
        if event_df["occurrences"].iloc[0] <= 1:
            print "Occurrence == 1, cancelling prediction."
            return False

        try:
            tag_predictions = self.decision_tree.predict(event_df[self.training_columns])
            return tag_predictions
        except NotFittedError:
            return False
        except ValueError as e:
            print "ERROR when predicting: %s" % e.message
            print "--Occurrences:", event_df["occurrences"]
            return False

    def predict_all_events_decision_tree(self, dataset):
        """
        Takes a dataset and runs predicitons on every splicing event.
        Returns predicted tags.
        """
        # TODO: Or return entire dataset?

        # First, exclude any event occurring in only 1 sample
        dataset = dataset.loc[dataset["occurrences"] > 1]

        # Also, only consider exon skipping events. # TODO: This should at some point be removed.
        dataset = dataset.loc[dataset["splice_type"] == "ES"]

        # Predict tags
        try:
            tag_predictions = self.decision_tree.predict(dataset[self.training_columns])
            dataset.loc[:, "decision_tree_tag"] = tag_predictions
            return dataset
        except NotFittedError:
            return False
        except ValueError as e:
            print "ERROR when predicting: %s" % e.message

            # TEST
            print "Attempting to find NaNs:"
            mask = dataset[self.training_columns].isnull().any()
            print mask

            sub = dataset[self.training_columns]
            print "\n10 first events with NaNs"
            print sub[sub.isnull().any(axis=1)].head(10)
            # END TEST

            return False



