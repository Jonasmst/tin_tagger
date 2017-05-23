import Tkinter as tk
import ttk
import sys
import re

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
    def __init__(self, *args, **kwargs):
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

        dummy_data = self.create_dummy_data()
        self.handle_row(dummy_data)

    def handle_row(self, data):
        # Get general data
        self.location_text["text"] = data["location"]
        self.splice_type_text["text"] = "%s (%s)" % (self.splice_type_map[data["splice_type"]], data["splice_type"])
        self.exons_text["text"] = data["exons"]
        self.sample_text = data["sample_of_interest"]
        self.gene_text["text"] = data["gene_symbol"]
        # TODO: self.strand_text["text"] = data["strand"]

        if data["splice_type"] == "AT":
            self.draw_alternative_terminator_event(data)

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

    def draw_alternative_terminator_event(self, data, wide_cov_mode=True):
        """
        Draws exons on the canvas, one row per sample, where each row contains
        PSI and coverage metrics.

        wide_cov_mode controls wether to draw wide coverage columns or not.
        """

        # Get list containing sample data
        samples_data = data["samples"]
        sample_of_interest = data["sample_of_interest"]
        exon_of_interest = data["exons"]

        # Keep track of number of samples and adjust vertical space for each samples accordingly
        number_of_samples = len(samples_data)

        # Specify dimensions
        window_height = self.window_height
        window_width = self.window_width

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
            # second_x_offset = 300

            # Draw a separator at the very bottom
            self.canvas.create_line(10, sample_start_y + height_per_sample, window_width - 10, sample_start_y + height_per_sample)

            # Highlight background if this is the sample in question
            background_color = COLOR_DARKWHITE
            if sample_name == sample_of_interest:
                # self.canvas.create_rectangle(0, sample_start_y, window_width, sample_start_y + height_per_sample, fill=COLOR_WHITE, outline=COLOR_RED)
                # self.canvas.create_rectangle(0, sample_start_y, window_width, sample_start_y + height_per_sample, fill=COLOR_LIGHTRED)
                background_color = COLOR_WHITE
            self.canvas.create_rectangle(0, sample_start_y, window_width, sample_start_y + height_per_sample, fill=background_color)

            # Print sample name on the left
            self.canvas.create_text(100, sample_start_y + (height_per_sample / 2), text=sample_name, width=100)

            #### BEGIN TEST #####
            # Draw exons in a loop
            for exon in data["exons"]:
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
                if wide_cov_mode:
                    frame_width = 50
                frame_start_x = 30 + x_offset
                if wide_cov_mode:
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
                if wide_cov_mode:
                    # self.canvas.create_text(frame_start_x + (frame_width / 2) + 1, start_y + 1, text=exon_coverage, width=150, fill=COLOR_WHITE, font="Helvetica 16")  # This is for text on top of coverage column
                    # self.canvas.create_text(frame_start_x + (frame_width / 2), start_y, text=exon_coverage, width=150, fill=COLOR_DARKBLUE, font="Helvetica 16")       # This is for text on top of coverage column
                    coverage_text_y_pos = start_y - 7
                    if percent_of_max_coverage >= 0.9:
                        # Draw text inside the frame, otherwise it'll overflow
                        coverage_text_y_pos = start_y + 10
                        # Also draw a white shadow, the background is all red
                        # self.canvas.create_text(frame_start_x + (frame_width / 2) + 1, coverage_text_y_pos + 1, text=exon_coverage, width=150, fill=COLOR_WHITE, font="Helvetica 16")       # This is for text on top of coverage column
                    self.canvas.create_text(frame_start_x + (frame_width / 2), coverage_text_y_pos, text=exon_coverage, width=150, fill="black", font="Helvetica 16")       # This is for text on top of coverage column
                else:
                    self.canvas.create_text(100 + x_offset - (85 + frame_width), poly_start + sample_start_y + (poly_height / 2), text=exon_coverage, width=150, fill=COLOR_DARKBLUE, font="Helvetica 16")

                # Write name of exon in the middle
                self.canvas.create_text(101 + x_offset, poly_start + sample_start_y + (poly_height / 2) + 1, text=exon_name, width=100, fill="black")  # Text shadow
                self.canvas.create_text(100 + x_offset, poly_start + sample_start_y + (poly_height / 2), text=exon_name, width=100, fill=COLOR_WHITE)
                # Write PSI below exon
                # self.canvas.create_text(100 + x_offset, poly_start + sample_start_y + poly_height + 10, text="%s" % exon_psi, fill=COLOR_DARKBLUE, font="Helvetica 16")
                self.canvas.create_text(100 + x_offset + 125, poly_start + sample_start_y + (poly_height / 2), text=exon_psi, fill=COLOR_DARKBLUE, font="Helvetica 16")

                x_offset += x_offset_increment

            sample_index += 1
            #### END TEST #####

            # # Create the left exon
            # self.canvas.create_polygon(
                # [
                    # 50 + x_offset, poly_start + sample_start_y,  # Start point upper left
                    # 50 + x_offset, poly_start + sample_start_y + poly_height,  # Down to the lower left corner
                    # 150 + x_offset, poly_start + sample_start_y + poly_height,  # Right to the lower right corner
                    # 150 + x_offset, poly_start + sample_start_y + ((poly_height / 4) * 3),  # Up to where the tail begins
                    # 200 + x_offset, poly_start + sample_start_y + ((poly_height / 4) * 3),  # Right to where the tail ends
                    # 200 + x_offset, poly_start + sample_start_y + (poly_height / 4),  # Up to the upper right corner of the tail
                    # 150 + x_offset, poly_start + sample_start_y + (poly_height / 4),  # Left to where the tail begins
                    # 150 + x_offset, poly_start + sample_start_y,  # Up to the upper right corner
                    # 50 + x_offset, poly_start + sample_start_y  # Back to where it started
                # ],
                # fill=COLOR_BLUE,
                # # outline=COLOR_RED
                # outline=COLOR_DARKBLUE
            # )

            # # Write coverage above exon
            # self.canvas.create_text(100 + x_offset, poly_start + sample_start_y - 10, text="Hello, this text", width=150, fill=COLOR_DARKBLUE, font="Helvetica 16 bold")
            # # Write name of exon in the middle
            # self.canvas.create_text(101 + x_offset, poly_start + sample_start_y + (poly_height / 2) + 1, text="Exon 1", width=100, fill="black")  # Text shadow
            # self.canvas.create_text(100 + x_offset, poly_start + sample_start_y + (poly_height / 2), text="Exon 1", width=100, fill=COLOR_WHITE)

            # # Create the right exon
            # total_x_offset = x_offset + second_x_offset
            # self.canvas.create_polygon(
                # [
                    # 50 + total_x_offset, poly_start + sample_start_y,  # Start point upper left
                    # 50 + total_x_offset, poly_start + sample_start_y + poly_height,  # Down to the lower left corner
                    # 150 + total_x_offset, poly_start + sample_start_y + poly_height,  # Right to the lower right corner
                    # 150 + total_x_offset, poly_start + sample_start_y + ((poly_height / 4) * 3),  # Up to where the tail begins
                    # 200 + total_x_offset, poly_start + sample_start_y + ((poly_height / 4) * 3),  # Right to where the tail ends
                    # 200 + total_x_offset, poly_start + sample_start_y + (poly_height / 4),  # Up to the upper right corner of the tail
                    # 150 + total_x_offset, poly_start + sample_start_y + (poly_height / 4),  # Left to where the tail begins
                    # 150 + total_x_offset, poly_start + sample_start_y,  # Up to the upper right corner
                    # 50 + total_x_offset, poly_start + sample_start_y  # Back to where it started
                # ],
                # fill="gray",
                # outline=COLOR_DARKBLUE
            # )

            # # Write text
            # self.canvas.create_text(100 + total_x_offset, poly_start + sample_start_y - 10, text="Coverage", width=150)
            # self.canvas.create_text(100 + total_x_offset, poly_start + sample_start_y + (poly_height / 2), text="Exon 2", width=100)

        # Finally, add "all" tag to all drawn objects
        self.canvas.addtag_all("all")

    def draw_sidebar(self):
        # Add a frame to keep everything in. Note that this does not autoscale very well
        right_frame = tk.Frame(self, bg=COLOR_WHITE, width=self.sidebar_width, borderwidth=2, relief="groove")
        right_frame.pack(fill=tk.BOTH, expand=tk.YES, side=tk.RIGHT)

        # Add a top label displaying information
        info_label = tk.Label(right_frame, text="Event information", font="Helvetica 20", anchor=tk.CENTER, background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        info_label.grid(row=0, column=0, columnspan=3, sticky="ew")
        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=1, column=0, columnspan=3, sticky="ew")

        # Width for all description labels
        labels_width = 15
        # Width for all text fields
        text_width = 35

        # Label and text fonts
        label_font = "Helvetica 16 bold"
        text_font = "Helvetica 16"

        # Location / coordinates
        location_label = tk.Label(right_frame, text="Location:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        location_label.grid(row=2, column=0, columnspan=1)
        coordinates_text = "1:12345678-12345679"
        # location_text = tk.Label(right_frame, text=coordinates_text, font="Helvetica 16", anchor=tk.W, background=COLOR_WHITE, width=text_width, borderwidth=2, relief="groove")
        self.location_text = tk.Label(right_frame, text=coordinates_text, font=text_font, anchor=tk.W, background=COLOR_WHITE)
        self.location_text.grid(row=2, column=1, sticky=tk.W)
        copy_button = ttk.Button(right_frame, text="Copy", command=lambda: self.copy_coords_to_clipboard(coordinates_text))
        copy_button.grid(row=2, column=2, sticky=tk.E)

        # Splice type
        splice_type_label = tk.Label(right_frame, text="Splice type:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        splice_type_label.grid(row=3, column=0, columnspan=1)
        splice_type = "ES"
        splice_type_description = "%s (%s)" % (self.splice_type_map[splice_type], splice_type)
        # splice_type_text = tk.Label(right_frame, text=splice_type_description, font="Helvetica 16", anchor=tk.W, background=COLOR_WHITE, width=text_width, borderwidth=2, relief="groove")
        self.splice_type_text = tk.Label(right_frame, text=splice_type_description, font=text_font, anchor=tk.W, background=COLOR_WHITE)
        self.splice_type_text.grid(row=3, column=1, sticky=tk.W)

        # Exons
        exons_label = tk.Label(right_frame, text="Exons:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        exons_label.grid(row=4, column=0, columnspan=1)
        exons = "1:2.1:3.2:4.3:4:5:6:7:8:11.1:12:1:2:3:4:123:6:5:3:123:2345:456:678"
        self.exons_text = tk.Label(right_frame, text=exons, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.exons_text.grid(row=4, column=1)

        # Gene symbol
        gene_label = tk.Label(right_frame, text="Gene:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        gene_label.grid(row=5, column=0, columnspan=1)
        gene = "TP53"
        self.gene_text = tk.Label(right_frame, text=gene, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.gene_text.grid(row=5, column=1)

        # Sample name
        sample_label = tk.Label(right_frame, text="Sample:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        sample_label.grid(row=6, column=0, columnspan=1)
        sample = "sample5"
        self.sample_text = tk.Label(right_frame, text=sample, font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        self.sample_text.grid(row=6, column=1)

        # Progress pane
        progress_label = tk.Label(right_frame, text="Progress", font="Helvetica 20", anchor=tk.CENTER, background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        progress_label.grid(row=7, column=0, columnspan=3, sticky="ew", pady=(100, 0))
        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=8, column=0, columnspan=3, sticky="ew")

        # Tagged yes/no/don't know
        tagged_width = 10
        tagged_labelframe = tk.LabelFrame(right_frame, text="Tagged events", padx=20, pady=20, background=COLOR_WHITE)
        tagged_labelframe.grid(row=9, column=0, columnspan=3, rowspan=2)
        yes_label = tk.Label(tagged_labelframe, text="0", font=text_font, bg=COLOR_WHITE, fg=COLOR_GREEN, borderwidth=1, relief="sunken", width=tagged_width, padx=10, pady=10)
        yes_label.grid(row=9, column=0, columnspan=1, sticky="nw")
        no_label = tk.Label(tagged_labelframe, text="0", font=text_font, width=tagged_width, padx=10, pady=10, bg=COLOR_WHITE, fg=COLOR_RED, borderwidth=1, relief="sunken")
        no_label.grid(row=9, column=1, columnspan=1, sticky="nw")
        uncertain_label = tk.Label(tagged_labelframe, text="0", font=text_font, width=tagged_width, padx=10, pady=10, bg=COLOR_WHITE, fg=COLOR_ORANGE, borderwidth=1, relief="sunken")
        uncertain_label.grid(row=9, column=2, columnspan=1, sticky="nw")
        yes_text = tk.Label(tagged_labelframe, text="Yes", font=text_font, bg=COLOR_WHITE, width=tagged_width, padx=10, pady=10)
        yes_text.grid(row=10, column=0, columnspan=1, sticky="nw")
        no_text = tk.Label(tagged_labelframe, text="No", font=text_font, bg=COLOR_WHITE, width=tagged_width, padx=10, pady=10)
        no_text.grid(row=10, column=1, columnspan=1, sticky="nw")
        uncertain_text = tk.Label(tagged_labelframe, text="Don't know", font=text_font, bg=COLOR_WHITE, width=tagged_width, padx=10, pady=10)
        uncertain_text.grid(row=10, column=2, columnspan=1, sticky="nw")

        # Saved information
        unsaved_label = tk.Label(right_frame, text="Unsaved tags:", font=label_font, anchor=tk.W, background=COLOR_WHITE, width=labels_width)
        unsaved_label.grid(row=11, column=0, columnspan=1)
        num_unsaved = 10
        unsaved_text = tk.Label(right_frame, text=str(num_unsaved), font=text_font, anchor=tk.W, background=COLOR_WHITE, width=text_width)
        unsaved_text.grid(row=11, column=1, columnspan=1, sticky="ew")
        save_button = ttk.Button(right_frame, text="Save")
        save_button.grid(row=11, column=2, sticky=tk.E)

        # Event tagging
        tagging_label = tk.Label(right_frame, text="Tag event", font="Helvetica 20", anchor=tk.CENTER, background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        tagging_label.grid(row=12, column=0, columnspan=3, sticky="ew", pady=(100, 0))
        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=13, column=0, columnspan=3, sticky="ew")
        button_interesting = ttk.Button(right_frame, text="Interesting")
        button_interesting.grid(row=13, column=0, sticky="ew")
        button_not_intersting = ttk.Button(right_frame, text="Not interesting")
        button_not_intersting.grid(row=13, column=1, sticky="ew")
        button_uncertain = ttk.Button(right_frame, text="Don't know")
        button_uncertain.grid(row=13, column=2, sticky="ew")

        ttk.Separator(right_frame, orient=tk.HORIZONTAL).grid(row=14, column=0, columnspan=3, sticky="ew", pady=(300, 0))
        # Previous tag + next and prev buttons
        previous_button = ttk.Button(right_frame, text="< previous")
        previous_button.grid(row=15, column=0, sticky="w")
        previous_tag = "N/A"
        previous_label = tk.Label(right_frame, text="Previous tag: %s" % previous_tag, font=text_font, anchor=tk.CENTER, background=COLOR_WHITE)
        previous_label.grid(row=15, column=1, sticky="sew")
        next_button = ttk.Button(right_frame, text="next >")
        next_button.grid(row=15, column=2, sticky="e")

    def key_pressed(self, event):
        """
        Handles various keypresses
        """
        print "Pressed:", repr(event.char)

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
    app = MainApplication()
    app.wm_title("Hello world, look at meeee")
    app.mainloop()
