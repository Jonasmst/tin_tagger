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
            "rpkm_percentage_of_mean_other_samples"
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
        print "Tagged events:", len(tagged_events)
        minimum_tagged_events = 20
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

        # TEST
        print "Training set:"
        print training_set
        print "\nValidation set:"
        print validation_set
        # END TEST

        # Fit tree
        self.decision_tree.fit(training_set[self.training_columns], training_set[self.tag_column])

        # Test accuracy of tree
        score = self.decision_tree.score(validation_set[self.training_columns], validation_set[self.tag_column])

        print "Decision tree prediction accuracy: %.2f" % score

        # Training successful, return accuracy
        return score

    def predict_tag_decision_tree(self, event_df):
        """
        Takes an event DataFrame and predicts the value of event_tag column

        :return: False if tree has not yet been fitted, np.array of tags if it has been fitted.
        """
        try:
            tag_predictions = self.decision_tree.predict(event_df[self.training_columns])
            return tag_predictions
        except NotFittedError:
            return False
