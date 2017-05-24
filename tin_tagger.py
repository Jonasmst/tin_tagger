import Tkinter as tk
import ttk
import sys
import re
import pandas as pd
import random
import os
import subprocess

COLOR_BLUE = "#3498db"
COLOR_YELLOW = "#f1c40f"
COLOR_RED = "#e74c3c"
COLOR_GREEN = "#2ecc71"
COLOR_WHITE = "#ecf0f1"
COLOR_DARKBLUE = "#34495e"
COLOR_PURPLE = "#9b59b6"
COLOR_LIGHTRED = "#ffe5e2"
COLOR_DARKWHITE = "#bdc3c7"
COLOR_DARKGRAY = "#7f8c8d"
COLOR_ORANGE = "#e67e22"

BAM_PATH = "/tsd/p19/data/durable/Projects/CRC-RNA-Seq/hisat2/transcriptome"


class ResizingCanvas(tk.Canvas):
    """
    This function resizes the canvas AND the objects within the canvas
    (except text). It is bound to the <Configure> event, which is invoked
    upon window resize (I think).
    """
    def __init__(self, parent, **kwargs):
        tk.Canvas.__init__(self, parent, **kwargs)
        self.bind("<Configure>", self.on_resize)
        self.height = self.winfo_reqheight()
        self.width = self.winfo_reqwidth()

    def on_resize(self, event):
        # Determine the ratio of old width/height to new width/height
        wscale = float(event.width) / self.width
        hscale = float(event.height) / self.height
        self.width = event.width
        self.height = event.height
        # Resize the canvas
        self.config(width=self.width, height=self.height)
        # Rescale all the objects tagged with the "all" tag
        self.scale("all", 0, 0, wscale, hscale)


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """ Provides natural sort as a key to sorted() """
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


class MainApplication(tk.Tk):
    def __init__(self, dataset_path, *args, **kwargs):
        # Super init
        tk.Tk.__init__(self, *args, **kwargs)

        # Bind various keypresses
        self.bind("<Key>", self.key_pressed)
        self.bind("<Escape>", self.close_window)

        ### Variables needed ###
        # Mapping exon skipping abbreviations to full text
        self.splice_type_map = {
            "ES": "Exon skipping",
            "AD": "Alternative donor site",
            "AA": "Alternative acceptor site",
            "RI": "Retained intron",
            "ME": "Mutually exclusive exons",
            "AT": "Alternative terminator",
            "AP": "Alternative promoter"
        }
        # BAM filepaths. TODO: Read this from file or something (*samplename*:*pathtobam*\n)
        self.bam_paths = {}
        for x in range(1, 11):
            self.bam_paths["sample%s" % str(x)] = os.path.join(BAM_PATH, "%s.sorted.bam" % str(x))
        # Pandas dataset
        self.dataset_path = dataset_path
        # Controls coverage graphics
        self.wide_cov_mode = False
        # Dataset
        self.dataset = pd.read_csv(self.dataset_path, sep="\t")
        # Current row
        self.current_row = 0  # Default to 0, controlled by prev/next buttons
        # Flag to define if we're just testing or running for real
        self.testing = False

        # Window sizes
        self.window_height = 1000
        self.window_width = 1000
        self.sidebar_width = 600

        # master.geometry("%dx%d+%d+%d" % (window_width+sidebar_width, window_height, window_width+500, window_height))

        # TODO: Can I grid() stuff instead of pack()-ing them?

        # Initialize a canvas
        canvas_frame = tk.Frame(self)
        canvas_frame.pack(fill=tk.BOTH, expand=tk.YES, side=tk.LEFT)
        self.canvas = ResizingCanvas(canvas_frame, width=self.window_width, height=self.window_height, bg="white", highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=tk.YES)

        # Draw event on canvas
        # self.draw_alternative_terminator_event(None)

        # Draw a sidebar
        self.draw_sidebar()

        # self.current_row_data = self.create_dummy_data()
        # dummy_data = self.create_dummy_data()
        self.current_row_data = self.get_row_data_by_index(self.current_row)
        self.handle_row(self.current_row_data)

    def next_button_pressed(self):
        """
        Handles next button presses. Loads the next row in the dataset.
        """
        # TODO: Passing self-args to functions doesn't make sense. Change get_row_data and handle_row functions
        self.current_row += 1
        if self.current_row > len(self.dataset):
            print "Wops: No more events (at end of dataset)"
            self.current_row -= 1

        # Load data for row
        self.current_row_data = self.get_row_data_by_index(self.current_row)
        self.handle_row(self.current_row_data)

    def prev_button_pressed(self):
        """
        Handles previous button presses. Loads the previous row in the dataset.
        """
        self.current_row -= 1
        if self.current_row < 0:
            print "Wops: No more previous events (at beginning of dataset)"
            self.current_row += 1
        self.current_row_data = self.get_row_data_by_index(self.current_row)
        self.handle_row(self.current_row_data)

    def get_row_data_by_index(self, row_index):
        """
        Returns formatted data for a given row index.
        """
        # TODO: Sanitize everything
        # TODO: Find included_counts and excluded_counts for other exons than the one in question

        # Get all sample names
        sample_names = list(self.dataset["name"].unique())

        splice_type = self.dataset.iloc[row_index]["splice_type"]
        sample_name = self.dataset.iloc[row_index]["name"]
        as_id = self.dataset.iloc[row_index]["as_id"]
        psi = self.dataset.iloc[row_index]["psi"]
        gene_symbol = self.dataset.iloc[row_index]["symbol"]
        strand = self.dataset.iloc[row_index]["strand"]
        exons = str(self.dataset.iloc[row_index]["exons"])
        chrom = self.dataset.iloc[row_index]["chr"]
        splice_start = self.dataset.iloc[row_index]["first_exon_start"]
        splice_stop = self.dataset.iloc[row_index]["last_exon_stop"]
        prev_exon_start = self.dataset.iloc[row_index]["prev_exon_start"]  # NaN for AT/AP
        prev_exon_stop = self.dataset.iloc[row_index]["prev_exon_stop"]  # NaN for AT/AP
        next_exon_start = self.dataset.iloc[row_index]["next_exon_start"]  # NaN for AT/AP
        next_exon_stop = self.dataset.iloc[row_index]["next_exon_stop"]  # NaN for AT/AP
        # Handle negative strand start- and stop- coordinates
        if strand == "-":
            splice_start = self.dataset.iloc[row_index]["last_exon_stop"]  # NaN for AT/AP
            splice_stop = self.dataset.iloc[row_index]["first_exon_start"]  # NaN for AT/AP
            prev_exon_start = self.dataset.iloc[row_index]["prev_exon_stop"]  # NaN for AT/AP
            prev_exon_stop = self.dataset.iloc[row_index]["prev_exon_start"]  # NaN for AT/AP
            next_exon_start = self.dataset.iloc[row_index]["next_exon_stop"]  # NaN for AT/AP
            next_exon_stop = self.dataset.iloc[row_index]["next_exon_start"]  # NaN for AT/AP
        included_counts = self.dataset.iloc[row_index]["included_counts"]
        excluded_counts = self.dataset.iloc[row_index]["excluded_counts"]
        prev_exon_name = str(self.dataset.iloc[row_index]["exon1"])  # NaN for AT/AP
        next_exon_name = str(self.dataset.iloc[row_index]["exon2"])  # NaN for AT/AP
        # prev_exon_id = self.dataset.iloc[row_index]["start_ex"]  # NaN for AT/AP
        # next_exon_id = self.dataset.iloc[row_index]["end_ex"]  # NaN for AT/AP

        # Create coordinates from chr, start and stop
        coordinates = str(chrom) + ":" + str(int(splice_start)) + "-" + str(int(splice_stop))

        # Find coverage for all exons in every sample
        exon_coverages = {}
        for s_name in sample_names:
            # First find for exon in question
            coverage = self.get_coverage_by_coordinates(s_name, coordinates)
            # included_counts, excluded_counts = self.get_included_and_excluded_counts
            if exons not in exon_coverages.keys():
                exon_coverages[exons] = {}
            exon_coverages[exons][s_name] = coverage

            # Then previous exon (if applicable)
            if not pd.isnull(prev_exon_name):
                if pd.isnull(prev_exon_start) or pd.isnull(prev_exon_stop):
                    print "Error: prev_exon_name is fine, but there's no prev_exon_start or prev_exon_stop"
                # Find coords
                prev_coords = str(chrom) + ":" + str(int(prev_exon_start)) + "-" + str(int(prev_exon_stop))
                prev_coverage = self.get_coverage_by_coordinates(s_name, prev_coords)
                if prev_exon_name not in exon_coverages.keys():
                    exon_coverages[prev_exon_name] = {}
                exon_coverages[prev_exon_name][s_name] = prev_coverage

            # Lastly, the next exon (if applicable)
            if not pd.isnull(next_exon_name):
                if pd.isnull(next_exon_start) or pd.isnull(next_exon_stop):
                    print "Error: next_exon_name is fine, but there's no next_exon_start or next_exon_stop"
                # Find coords
                next_coords = str(chrom) + ":" + str(int(next_exon_start)) + "-" + str(int(next_exon_stop))
                next_coverage = self.get_coverage_by_coordinates(s_name, next_coords)
                if next_exon_name not in exon_coverages.keys():
                    exon_coverages[next_exon_name] = {}
                exon_coverages[next_exon_name][s_name] = next_coverage

        # General row data
        row_data = {
            "splice_type": splice_type,
            "gene_symbol": gene_symbol,
            "sample_of_interest": sample_name,
            "location": coordinates,
            "exons": exons,
            "strand": strand,
            "as_id": as_id,
            "exon_psi": str(psi)
        }

        # Sample-specific data
        samples_data = {}
        for s_name in sample_names:
            if s_name not in samples_data.keys():
                samples_data[s_name] = {}
            samples_data[s_name]["gene_rpkm"] = self.get_gene_rpkm_by_sample(s_name, gene_symbol)
            # Find exons in this sample
            sample_exons = {}
            for e_name in exon_coverages.keys():
                if e_name not in sample_exons.keys():
                    sample_exons[e_name] = {}
                sample_exons[e_name]["coverage"] = exon_coverages[e_name][s_name]
                if s_name == sample_name and e_name == exons:
                    sample_exons[e_name]["psi"] = psi
                else:
                    sample_exons[e_name]["psi"] = "N/A"
                sample_exons[e_name]["max_coverage"] = max(exon_coverages[e_name].values())

            # Add exon data to sample
            samples_data[s_name]["exons"] = sample_exons

        row_data["samples"] = samples_data

        return row_data

    def get_coverage_by_coordinates(self, sample_name, coordinates):
        """
        Runs samtools to find coverage for region in a sample.
        """
        # TODO: Get path to file, run samtools via subprocess, calc avg. coverage for region and return
        # For now, just return a random int.
        # return random.randint(100, 2000)

        if self.testing:
            return random.randint(100, 2000)
        else:
            path_to_bam = self.bam_paths[sample_name]
            command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++;}END{if (cnt>0){ print sum/cnt } else print 0}'" % (coordinates, path_to_bam)
            # TODO: Something about security issue with shell=True, but it allows piping in awk commands, so I'll keep it for now
            samtools_output = subprocess.check_output(command, shell=True)

            region_coverage = -1.0
            try:
                region_coverage = float(samtools_output)
            except:
                print "ERROR: Samtools output can't be converted to float:"
                print samtools_output

            return region_coverage

    def get_gene_rpkm_by_sample(self, sample_name, gene_name):
        """
        Returns the RPKM of a gene for a given sample
        """
        df = self.dataset.loc[(self.dataset["symbol"] == gene_name) & (self.dataset["name"] == sample_name)]
        # Handle genes that are not present in every sample
        rpkm = -1.0
        if len(df) > 0:
            rpkm = df["rpkm"].iloc[0]
        return rpkm

    def handle_row(self, data):
        """
        Controls the drawing of events. Populates sidebar and how the canvas is drawn.
        """
        # Populate the sidebar with general information
        self.location_text["text"] = data["location"]
        self.splice_type_text["text"] = "%s (%s)" % (self.splice_type_map[data["splice_type"]], data["splice_type"])
        self.exons_text["text"] = data["exons"]
        self.sample_text["text"] = data["sample_of_interest"]
        self.gene_text["text"] = data["gene_symbol"]
        self.strand_text["text"] = data["strand"]
        self.psi_text["text"] = data["exon_psi"]

        # Populate canvas
        if data["splice_type"] == "AT":
            self.draw_alternative_terminator_event(data)
        elif data["splice_type"] == "ES":
            self.draw_exon_skipping_event(data)

    def create_dummy_data(self):
        # Create a dummy json-esque data container
        dummy_data = {
            "splice_type": "AT",
            "gene_symbol": "TP53",
            "sample_of_interest": "sample5",
            "location": "1:12345678-123456789",
            "exons": "10",
            "strand": "+",
            "samples": {
                "sample1": {
                    "gene_rpkm": 2132,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2315,
                        "psi": 0.97,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 123,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }
                    ],
                },
                "sample2": {
                    "gene_rpkm": 2300,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2120,
                        "psi": 0.99,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 93,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }
                    ]
                },
                "sample3": {
                    "gene_rpkm": 2310,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2202,
                        "psi": 0.97,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 130,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }]
                },
                "sample4": {
                    "gene_rpkm": 2420,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2211,
                        "psi": 0.98,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 123,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }]
                },
                "sample5": {
                    "gene_rpkm": 2123,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 122,
                        "psi": 0.01,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 2221,
                        "psi": 0.99,
                        "max_coverage": 2221
                    }],
                },
                "sample6": {
                    "gene_rpkm": 2032,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2202,
                        "psi": 0.97,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 123,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }]
                },
                "sample7": {
                    "gene_rpkm": 1588,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2310,
                        "psi": 0.99,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 211,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }]
                },
                "sample8": {
                    "gene_rpkm": 2299,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2120,
                        "psi": 0.99,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 123,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }]
                },
                "sample9": {
                    "gene_rpkm": 2009,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 2010,
                        "psi": 0.99,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 231,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }
                    ]
                },
                "sample10": {
                    "gene_rpkm": 2100,
                    "exons": [{
                        "exon_name": "10",
                        "coverage": 3102,
                        "psi": 0.99,
                        "max_coverage": 3102,
                    }, {
                        "exon_name": "11",
                        "coverage": 220,
                        "psi": 0.01,
                        "max_coverage": 2221
                    }]
                }
            },
        }

        return dummy_data

    def draw_alternative_terminator_event(self, data):
        """
        Draws exons on the canvas, one row per sample, where each row contains
        PSI and coverage metrics.

        wide_cov_mode controls whether to draw wide coverage columns or not.
        """

        # First, wipe the canvas
        self.canvas.delete("all")

        # Get list containing sample data
        samples_data = data["samples"]
        sample_of_interest = data["sample_of_interest"]
        exon_of_interest = data["exons"]

        # Keep track of number of samples and adjust vertical space for each samples accordingly
        number_of_samples = len(samples_data)

        # Specify dimensions
        window_height = self.canvas.height
        window_width = self.canvas.width

        # Calc how much vertical space we have per sample
        height_per_sample = window_height / number_of_samples

        # Keep track of sample index
        sample_index = 0

        # Loop the samples
        # for x in range(len(samples_data)):
        for sample_name in sorted(samples_data.keys(), key=natural_sort_key):
            # Get sample data
            data = samples_data[sample_name]
            # Define a start y-coordinate for this sample
            sample_start_y = sample_index * height_per_sample

            # Define some dimesions
            poly_start = height_per_sample / 4
            poly_height = height_per_sample / 2
            x_offset = 200
            x_offset_increment = 400

            # Draw a separator at the very bottom
            self.canvas.create_line(10, sample_start_y + height_per_sample, window_width - 10, sample_start_y + height_per_sample)

            # Highlight background if this is the sample in question
            background_color = COLOR_DARKWHITE
            if sample_name == sample_of_interest:
                background_color = COLOR_WHITE
            self.canvas.create_rectangle(0, sample_start_y, window_width, sample_start_y + height_per_sample, fill=background_color)

            # Print sample name on the left
            self.canvas.create_text(100, sample_start_y + (height_per_sample * 0.25), text=sample_name, width=100)
            # Print gene RPKM on the left too.
            # TODO: Find out a good placement for RPKM values
            rpkm_x = 45
            rpkm_y = sample_start_y + height_per_sample - 10
            self.canvas.create_text(rpkm_x, rpkm_y, text="RPKM: %s" % data["gene_rpkm"], width=100)

            # Draw exons in a loop
            for exon in data["exons"]:
                # TODO: Draw exons in order
                exon_name = exon["exon_name"]
                exon_coverage = str(exon["coverage"])
                exon_psi = str(exon["psi"])

                # Draw exon
                outline_color = COLOR_DARKBLUE
                outline_width = 1
                if exon_name == exon_of_interest and sample_name == sample_of_interest:
                    outline_color = COLOR_DARKBLUE
                    outline_width = 3
                self.canvas.create_polygon(
                    [
                        50 + x_offset, poly_start + sample_start_y,  # Start point upper left
                        50 + x_offset, poly_start + sample_start_y + poly_height,  # Down to the lower left corner
                        150 + x_offset, poly_start + sample_start_y + poly_height,  # Right to the lower right corner
                        150 + x_offset, poly_start + sample_start_y + ((poly_height / 4) * 3),  # Up to where the tail begins
                        200 + x_offset, poly_start + sample_start_y + ((poly_height / 4) * 3),  # Right to where the tail ends
                        200 + x_offset, poly_start + sample_start_y + (poly_height / 4),  # Up to the upper right corner of the tail
                        150 + x_offset, poly_start + sample_start_y + (poly_height / 4),  # Left to where the tail begins
                        150 + x_offset, poly_start + sample_start_y,  # Up to the upper right corner
                        50 + x_offset, poly_start + sample_start_y  # Back to where it started
                    ],
                    fill=COLOR_BLUE,
                    # outline=COLOR_RED
                    outline=outline_color,
                    width=outline_width
                )

                # Calc percentage of max coverage for this exon
                percent_of_max_coverage = float(exon["coverage"]) / float(exon["max_coverage"])
                # Draw coverage frame
                frame_y_padding = 5
                frame_border_width = 1
                frame_width = 10
                if self.wide_cov_mode:
                    frame_width = 50
                frame_start_x = 30 + x_offset
                if self.wide_cov_mode:
                    frame_start_x = x_offset - 10
                frame_stop_x = frame_start_x + frame_width

                # Draw coverage frame
                frame_start_y = sample_start_y + frame_y_padding
                frame_stop_y = sample_start_y + height_per_sample - frame_y_padding
                frame_height = frame_stop_y - frame_start_y - (2 * frame_border_width)
                height_of_coverage = float(frame_height) * float(percent_of_max_coverage)
                self.canvas.create_rectangle(frame_start_x, frame_start_y, frame_stop_x, frame_stop_y, width=frame_border_width, outline="black")

                # Draw coverage fill
                start_x = frame_start_x + frame_border_width
                start_y = frame_stop_y - frame_border_width - height_of_coverage
                stop_x = frame_stop_x - frame_border_width
                stop_y = frame_stop_y - frame_border_width
                self.canvas.create_rectangle(start_x, start_y, stop_x, stop_y, outline=COLOR_RED, fill=COLOR_RED)

                # Write coverage to the left of exon
                if self.wide_cov_mode:
                    coverage_text_y_pos = start_y - 7
                    if percent_of_max_coverage >= 0.9:
                        # Draw text inside the frame, otherwise it'll overflow
                        coverage_text_y_pos = start_y + 10
                        # Also draw a white shadow, the background is all red
                    self.canvas.create_text(frame_start_x + (frame_width / 2), coverage_text_y_pos, text=exon_coverage, width=150, fill="black", font="Helvetica 16")       # This is for text on top of coverage column
                else:
                    self.canvas.create_text(100 + x_offset - (85 + frame_width), poly_start + sample_start_y + (poly_height / 2), text=exon_coverage, width=150, fill=COLOR_DARKBLUE, font="Helvetica 16")

                # Write name of exon in the middle
                self.canvas.create_text(101 + x_offset, poly_start + sample_start_y + (poly_height / 2) + 1, text=exon_name, width=100, fill="black")  # Text shadow
                self.canvas.create_text(100 + x_offset, poly_start + sample_start_y + (poly_height / 2), text=exon_name, width=100, fill=COLOR_WHITE)
                # Write PSI below exon
                self.canvas.create_text(100 + x_offset + 125, poly_start + sample_start_y + (poly_height / 2), text=exon_psi, fill=COLOR_DARKBLUE, font="Helvetica 16")

                x_offset += x_offset_increment

            sample_index += 1

        # Finally, add "all" tag to all drawn objects
        self.canvas.addtag_all("all")

    def draw_exon_skipping_event(self, data):
        """
        Draws exon skipping event on canvas
        """

        # First, wipe the canvas
        self.canvas.delete("all")

        # Get list containing sample data
        samples_data = data["samples"]
        sample_of_interest = data["sample_of_interest"]
        exon_of_interest = data["exons"]

        # Keep track of number of samples and adjust vertical space for each samples accordingly
        number_of_samples = len(samples_data)

        # Specify dimensions
        window_height = self.canvas.height
        window_width = self.canvas.width

        # Calc how much vertical space we have per sample
        height_per_sample = window_height / number_of_samples

        # Keep track of sample index
        sample_index = 0

        # Loop the samples
        for sample_name in sorted(samples_data.keys(), key=natural_sort_key):
            # Get sample data
            data = samples_data[sample_name]
            # Define a start y-coordinate for this sample
            sample_start_y = sample_index * height_per_sample

            # Define some dimesions
            poly_start_y = height_per_sample / 4
            poly_height = height_per_sample / 2
            poly_width = 100
            poly_start_x = 45
            x_offset = 200
            x_offset_increment = 200
            left_text_x = 10
            left_text_width = 100

            # Draw a separator at the very bottom
            self.canvas.create_line(10, sample_start_y + height_per_sample, window_width - 10, sample_start_y + height_per_sample)

            # Highlight background if this is the sample in question
            background_color = COLOR_DARKWHITE
            if sample_name == sample_of_interest:
                background_color = COLOR_WHITE
            self.canvas.create_rectangle(0, sample_start_y, window_width, sample_start_y + height_per_sample, fill=background_color)

            # Print sample name on the left
            self.canvas.create_text(left_text_x, sample_start_y + (height_per_sample * 0.25), text=sample_name, width=left_text_width, anchor=tk.W)
            # Print gene RPKM on the left too. Only 2 decimals
            rpkm = data["gene_rpkm"]
            rpkm_text = "RPKM: N/A"
            if rpkm > 0:
                rpkm_text = "RPKM: %.2f" % rpkm
            self.canvas.create_text(left_text_x, sample_start_y + (height_per_sample * 0.75), text=rpkm_text, width=left_text_width, anchor=tk.W)

            #################################################
            #################### DRAW EXONS #################
            #################################################
            for exon_name in sorted(data["exons"].keys(), key=natural_sort_key):
                exon = data["exons"][exon_name]
                exon_coverage = str(exon["coverage"])
                # exon_psi = str(exon["psi"])

                # Define exon color
                outline_color = COLOR_DARKBLUE
                outline_width = 1
                if exon_name == exon_of_interest and sample_name == sample_of_interest:
                    outline_color = COLOR_DARKBLUE
                    outline_width = 3

                # Draw exon
                self.canvas.create_rectangle(
                    poly_start_x + x_offset, sample_start_y + poly_start_y,
                    poly_start_x + poly_width + x_offset, sample_start_y + poly_start_y + poly_height,
                    outline=outline_color,
                    width=outline_width,
                    fill=COLOR_BLUE
                )

                ##############################################
                ############## COVERAGE STUFF ################
                ##############################################
                # Calc percentage of max coverage for this exon
                if not float(exon["max_coverage"]) > 0.0:
                    percent_of_max_coverage = 0.00
                else:
                    percent_of_max_coverage = float(exon["coverage"]) / float(exon["max_coverage"])
                # Draw coverage frame
                frame_y_padding = 5
                frame_border_width = 1
                frame_width = 10
                if self.wide_cov_mode:
                    frame_width = 50
                frame_start_x = 30 + x_offset
                if self.wide_cov_mode:
                    frame_start_x = x_offset - 10
                frame_stop_x = frame_start_x + frame_width

                # Draw coverage frame
                frame_start_y = sample_start_y + frame_y_padding
                frame_stop_y = sample_start_y + height_per_sample - frame_y_padding
                frame_height = frame_stop_y - frame_start_y - (2 * frame_border_width)
                height_of_coverage = float(frame_height) * float(percent_of_max_coverage)
                self.canvas.create_rectangle(frame_start_x, frame_start_y, frame_stop_x, frame_stop_y, width=frame_border_width, outline="black")

                # Draw coverage fill
                start_x = frame_start_x + frame_border_width
                start_y = frame_stop_y - frame_border_width - height_of_coverage
                stop_x = frame_stop_x - frame_border_width
                stop_y = frame_stop_y - frame_border_width
                self.canvas.create_rectangle(start_x, start_y, stop_x, stop_y, outline=COLOR_RED, fill=COLOR_RED)

                # Write coverage to the left of exon
                if self.wide_cov_mode:
                    coverage_text_y_pos = start_y - 7
                    if percent_of_max_coverage >= 0.9:
                        # Draw text inside the frame, otherwise it'll overflow
                        coverage_text_y_pos = start_y + 10
                        # Also draw a white shadow, the background is all red
                    self.canvas.create_text(frame_start_x + (frame_width / 2), coverage_text_y_pos, text=exon_coverage, width=150, fill="black", font="Helvetica 16")       # This is for text on top of coverage column
                else:
                    self.canvas.create_text(100 + x_offset - (85 + frame_width), poly_start_y + sample_start_y + (poly_height / 2), text=exon_coverage, width=150, fill=COLOR_DARKBLUE, font="Helvetica 16")

                # Write name of exon in the middle
                self.canvas.create_text(poly_start_x + x_offset + (poly_width / 2) + 1, poly_start_y + sample_start_y + (poly_height / 2) + 1, text=exon_name, width=100, fill="black")  # Text shadow
                self.canvas.create_text(poly_start_x + x_offset + (poly_width / 2), poly_start_y + sample_start_y + (poly_height / 2), text=exon_name, width=100, fill=COLOR_WHITE)
                # Write PSI below exon
                # self.canvas.create_text(100 + x_offset + 125, poly_start_y+ sample_start_y + (poly_height / 2), text=exon_psi, fill=COLOR_DARKBLUE, font="Helvetica 16")

                x_offset += x_offset_increment

            sample_index += 1

        # Finally, add "all" tag to all drawn objects
        self.canvas.addtag_all("all")

    def draw_sidebar(self):
        # Add a frame to keep everything in. Note that this does not autoscale very well
        right_frame = tk.Frame(self, bg=COLOR_WHITE, width=self.sidebar_width, borderwidth=2, relief="groove")
        right_frame.pack(fill=tk.BOTH, expand=tk.YES, side=tk.RIGHT)

        # Width for all description labels
        labels_width = 15
        # Width for all text fields
        text_width = 35

        # Keep track of rows, to make shuffling rows easier
        current_row = 0

        # Add a top label displaying information
        info_label = tk.Label(right_frame, text="Event information", font="Helvetica 20", anchor=tk.CENTER, background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        info_label.grid(row=current_row, column=0, columnspan=3, sticky="ew")
        current_row += 1
        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=current_row, column=0, columnspan=3, sticky="ew")
        current_row += 1

        # Label and text fonts
        label_font = "Helvetica 16 bold"
        text_font = "Helvetica 16"

        # Location / coordinates
        location_label = tk.Label(right_frame, text="Location:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        location_label.grid(row=current_row, column=0, columnspan=1)
        coordinates_text = "coordinates"
        # location_text = tk.Label(right_frame, text=coordinates_text, font="Helvetica 16", anchor=tk.W, background=COLOR_WHITE, width=text_width, borderwidth=2, relief="groove")
        self.location_text = tk.Label(right_frame, text=coordinates_text, font=text_font, anchor=tk.W, background=COLOR_WHITE)
        self.location_text.grid(row=current_row, column=1, sticky=tk.W)
        copy_button = ttk.Button(right_frame, text="Copy", command=lambda: self.copy_coords_to_clipboard(coordinates_text))
        copy_button.grid(row=current_row, column=2, sticky=tk.E)
        current_row += 1

        # Splice type
        splice_type_label = tk.Label(right_frame, text="Splice type:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        splice_type_label.grid(row=current_row, column=0, columnspan=1)
        splice_type = "Splice type"
        self.splice_type_text = tk.Label(right_frame, text=splice_type, font=text_font, anchor=tk.W, background=COLOR_WHITE)
        self.splice_type_text.grid(row=current_row, column=1, sticky=tk.W)
        current_row += 1

        # Exons
        exons_label = tk.Label(right_frame, text="Exons:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        exons_label.grid(row=current_row, column=0, columnspan=1)
        exons = "exons"
        self.exons_text = tk.Label(right_frame, text=exons, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.exons_text.grid(row=current_row, column=1)
        current_row += 1

        # Gene symbol
        gene_label = tk.Label(right_frame, text="Gene:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        gene_label.grid(row=current_row, column=0, columnspan=1)
        gene = "gene"
        self.gene_text = tk.Label(right_frame, text=gene, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.gene_text.grid(row=current_row, column=1)
        current_row += 1

        # Strand
        strand_label = tk.Label(right_frame, text="Strand:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        strand_label.grid(row=current_row, column=0, columnspan=1)
        strand = "strand"
        self.strand_text = tk.Label(right_frame, text=strand, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.strand_text.grid(row=current_row, column=1)
        current_row += 1

        # Sample name
        sample_label = tk.Label(right_frame, text="Sample:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        sample_label.grid(row=current_row, column=0, columnspan=1)
        sample = "sample"
        self.sample_text = tk.Label(right_frame, text=sample, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.sample_text.grid(row=current_row, column=1)
        current_row += 1

        # PSI
        psi_label = tk.Label(right_frame, text="PSI:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        psi_label.grid(row=current_row, column=0, columnspan=1)
        psi = "PSI"
        self.psi_text = tk.Label(right_frame, text=psi, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.psi_text.grid(row=current_row, column=1)
        current_row += 1

        # Progress pane
        progress_label = tk.Label(right_frame, text="Progress", font="Helvetica 20", anchor=tk.CENTER, background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        progress_label.grid(row=current_row, column=0, columnspan=3, sticky="ew", pady=(100, 0))
        current_row += 1
        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=current_row, column=0, columnspan=3, sticky="ew")
        current_row += 1

        # Tagged yes/no/don't know
        tagged_width = 10
        tagged_labelframe = tk.LabelFrame(right_frame, text="Tagged events", padx=20, pady=20, background=COLOR_WHITE)
        tagged_labelframe.grid(row=current_row, column=0, columnspan=3, rowspan=2)
        yes_label = tk.Label(tagged_labelframe, text="0", font=text_font, bg=COLOR_WHITE, fg=COLOR_GREEN, borderwidth=1, relief="sunken", width=tagged_width, padx=10, pady=10)
        yes_label.grid(row=current_row, column=0, columnspan=1, sticky="nw")
        no_label = tk.Label(tagged_labelframe, text="0", font=text_font, width=tagged_width, padx=10, pady=10, bg=COLOR_WHITE, fg=COLOR_RED, borderwidth=1, relief="sunken")
        no_label.grid(row=current_row, column=1, columnspan=1, sticky="nw")
        uncertain_label = tk.Label(tagged_labelframe, text="0", font=text_font, width=tagged_width, padx=10, pady=10, bg=COLOR_WHITE, fg=COLOR_ORANGE, borderwidth=1, relief="sunken")
        uncertain_label.grid(row=current_row, column=2, columnspan=1, sticky="nw")
        current_row += 1
        yes_text = tk.Label(tagged_labelframe, text="Yes", font=text_font, bg=COLOR_WHITE, width=tagged_width, padx=10, pady=10)
        yes_text.grid(row=current_row, column=0, columnspan=1, sticky="nw")
        no_text = tk.Label(tagged_labelframe, text="No", font=text_font, bg=COLOR_WHITE, width=tagged_width, padx=10, pady=10)
        no_text.grid(row=current_row, column=1, columnspan=1, sticky="nw")
        uncertain_text = tk.Label(tagged_labelframe, text="Don't know", font=text_font, bg=COLOR_WHITE, width=tagged_width, padx=10, pady=10)
        uncertain_text.grid(row=current_row, column=2, columnspan=1, sticky="nw")
        current_row += 1

        # Saved information
        unsaved_label = tk.Label(right_frame, text="Unsaved tags:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        unsaved_label.grid(row=current_row, column=0, columnspan=1)
        num_unsaved = 10
        unsaved_text = tk.Label(right_frame, text=str(num_unsaved), font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        unsaved_text.grid(row=current_row, column=1, columnspan=1, sticky="ew")
        save_button = ttk.Button(right_frame, text="Save")
        save_button.grid(row=current_row, column=2, sticky=tk.E)
        current_row += 1

        # Event tagging
        tagging_label = tk.Label(right_frame, text="Tag event", font="Helvetica 20", anchor=tk.CENTER, background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        tagging_label.grid(row=current_row, column=0, columnspan=3, sticky="ew", pady=(100, 0))
        current_row += 1
        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=current_row, column=0, columnspan=3, sticky="ew")
        button_interesting = ttk.Button(right_frame, text="Interesting")
        button_interesting.grid(row=current_row, column=0, sticky="ew")
        button_not_intersting = ttk.Button(right_frame, text="Not interesting")
        button_not_intersting.grid(row=current_row, column=1, sticky="ew")
        button_uncertain = ttk.Button(right_frame, text="Don't know")
        button_uncertain.grid(row=current_row, column=2, sticky="ew")
        current_row += 1

        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=current_row, column=0, columnspan=3, sticky="ew", pady=(300, 0))
        current_row += 1
        # Previous tag + next and prev buttons
        previous_button = ttk.Button(right_frame, text="< previous", command=self.prev_button_pressed)
        previous_button.grid(row=current_row, column=0, sticky="w")
        previous_tag = "N/A"
        previous_label = tk.Label(right_frame, text="Previous tag: %s" % previous_tag, font=text_font, anchor=tk.CENTER, background=COLOR_WHITE)
        previous_label.grid(row=current_row, column=1, sticky="sew")
        next_button = ttk.Button(right_frame, text="next >", command=self.next_button_pressed)
        next_button.grid(row=current_row, column=2, sticky="e")
        current_row += 1

    def key_pressed(self, event):
        """
        Handles various keypresses
        """
        print "Pressed:", repr(event.char)
        if event.char == '1':
            # Redraw canvas with new color
            self.wide_cov_mode = not self.wide_cov_mode
            self.handle_row(self.current_row_data)
        if event.char == '2':
            self.canvas.delete("all")

    def copy_coords_to_clipboard(self, coords):
        """
        Copies the coordinates from the coordinates-label to the
        system clipboard
        """
        self.clipboard_clear()
        self.clipboard_append(coords)

    def close_window(self, event):
        sys.exit()


#############################
########### RUN #############
#############################
if __name__ == "__main__":
    # Get path to dataset from arguments
    # dataset_filepath = sys.argv[1]
    dataset_filepath = "/Users/jonas/Dropbox/phd/code/tin_tagger/datasets/es_only.tsv"
    app = MainApplication(dataset_filepath)
    app.wm_title("Hello world, look at meeee")
    app.mainloop()
