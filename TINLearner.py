
class TINLearner(object):

    def __init__(self, tin_tagger):
        self.tin_tagger = tin_tagger
        self.dataset = None

    def prepare_dataset(self, dataset):
        """
        Prepares dataset for training. Maps categoricals to numericals, ensures all appropriate columns are present
        in the dataset.
        """

        """
        PARAMETERS:
        -----------
        - Main exon RPKM (combined if multiple exons)
        - Upstream exon RPKM
        - Downstream exon RPKM
        - Main exon RPKM diff from upstream exon RPKM
        - Main exon RPKM diff from downstream exon RPKM
        - Main exon diff from mean main exon RPKM in other samples
        - Main exon PSI diff from main exon PSI in other samples
        - Upstream exon RPKM diff from mean upstream RPKM in other samples
        - Downstream exon RPKM diff from mean downstream RPKM in other samples
        - Main exon included counts
        - Main exon excluded counts
        - Diff in main exon included counts from mean main exon included counts in other samples
        - Diff in main exon excluded counts from mean main exon excluded coutns in other samples
        - Splice type mapped to numeric
        - Is_reported mapped to numeric?
        """

        # Present columns to include in training
        present_keep_columns = [
            "included_counts",
            "excluded_counts",
            "psi",
            "rpkm"
        ]

        # Generate (or update) new columns
        dataset["upstream_exon_rpkm"]