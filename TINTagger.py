import Tkinter as tk
import ttk
import sys
import re
from tkFileDialog import askopenfilename, asksaveasfilename
import tkFont
import pandas as pd
import numpy as np
import os
import copy
import random
import subprocess
import tkMessageBox
from TINDataProcessor import TINDataProcessor
from multiprocessing import Process, Queue
from Queue import Empty

# TODO: Cache rows so we don't need to run samtools etc when pressing previous-button
# TODO: Find PSI, included counts, excluded counts for other exons and display if it's available (gonna be an sql-call).
# TODO: Pre-fetch a number of rows in the background (e.g. fetch samtools-coverage/spliceseq DB-stuff for e.g. 100 rows)
# TODO: Instead of highlighting the sample canvas, we should highlight the exon column (not the unreported ones, though).
# TODO: Sort dataset by gene, then "exons" to show similar events after each other?
# TODO: Text is huge on linux
# TODO: Add feedback as text in the statusbar when connecting to DB and reading datasets. Blink/animate (non-static)
# TODO: Support both reading from database and reading from local file

"""
################################################
Main frame:
                    Column 0        Column 1
    Row 0           left_frame      sidebar
################################################
Left frame
                    Column 0        Column 1        Column 2?
    Row 0           sample_frame    exon_frame      buttons_frame?
################################################
Exon frame
                    Column 0
    Row 0       exon_name_frame
    Row 1         exon_canvas
                      ...
    Row n         exon_canvas
################################################
"""


COLOR_BLUE = "#3498db"
COLOR_YELLOW = "#f1c40f"
COLOR_RED = "#e74c3c"
COLOR_GREEN = "#2ecc71"
COLOR_WHITE = "#ecf0f1"
COLOR_DARKBLUE = "#34495e"
COLOR_PURPLE = "#9b59b6"
COLOR_LIGHTRED = "#ffe5e2"
COLOR_DARKWHITE = "#bdc3c7"
COLOR_GRAY = "#95a5a6"
COLOR_DARKGRAY = "#7f8c8d"
COLOR_ORANGE = "#e67e22"
COLOR_CANVAS_TEXT = COLOR_WHITE
COLOR_CANVAS_TEXT_SHADOW = "black"
COLOR_SAMPLE_HIGHLIGHT = COLOR_LIGHTRED
COLOR_EXON_NAME = "black"

COLOR_INTERESTING = COLOR_GREEN
COLOR_NOT_INTERESTING = COLOR_RED
COLOR_UNCERTAIN = COLOR_PURPLE
COLOR_NO_TAG = COLOR_DARKGRAY

TAG_INTERESTING = 0
TAG_NOT_INTERESTING = 1
TAG_UNCERTAIN = 2
TAG_NO_TAG = -1

TEXTTAG_COVERAGE = "exon_rpkm_text"
TEXTTAG_SHADOW = "text_shadow"

STYLE_BUTTON_INTERESTING_ON = "InterestingOn.TButton"
STYLE_BUTTON_NOT_INTERESTING_ON = "Not_interestingOn.TButton"
STYLE_BUTTON_UNCERTAIN_ON = "UncertainOn.TButton"
STYLE_BUTTON_INTERESTING_OFF = "InterestingOff.TButton"
STYLE_BUTTON_NOT_INTERESTING_OFF = "Not_interestingOff.TButton"
STYLE_BUTTON_UNCERTAIN_OFF = "UncertainOff.TButton"


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """ Provides natural sort as a key to sorted() """
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


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

        # Redraw text and text shadow to keep the 1-pixel difference upon resize
        # Delete all text shadows
        self.delete(TEXTTAG_SHADOW)
        for item in self.find_withtag(TEXTTAG_COVERAGE):
            if self.type(item) == "text":
                # Get current font
                font = self.itemcget(item, "font")

                # Get the current text
                current_text = self.itemcget(item, "text")

                # Get coordinates of the text object
                text_x, text_y = self.coords(item)

                # Draw shadow first, then a new instance of the text
                self.create_text(text_x + 1, text_y + 1, text=current_text, tags=TEXTTAG_SHADOW, font=font, fill=COLOR_CANVAS_TEXT_SHADOW)
                self.create_text(text_x, text_y, text=current_text, font=font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

                # Finally, delete this item
                self.delete(item)


class TINTagger(tk.Tk):

    def __init__(self, *args, **kwargs):
        # Super init
        tk.Tk.__init__(self, *args, **kwargs)

        # Different themes
        self.style = ttk.Style()
        self.available_themes = self.style.theme_names()
        self.current_theme = tk.StringVar(self)
        self.current_theme.set(self.style.theme_use())
        self.current_theme.trace("w", self.change_theme)

        # TEST
        self.all_uninteresting_buttons = []
        # END TEST

        # Text size for all canvas text items
        self.canvas_text_size = 16
        self.canvas_font = ("tkDefaultFont", self.canvas_text_size)

        # Set testing or not
        self.testing = True

        if self.testing:
            # Mac ghetto-fix for bringing GUI to front when running.
            #os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
            pass

        #####################
        # Various variables #
        #####################
        # Keep track of all of the canvases in the center window
        self.canvases = []

        # No dataset loaded by default
        self.dataset = None
        # Keep a copy of the original dataset to use when filtering
        self.original_dataset = None

        # Keep track of whether we're currently waiting for a dataset to be read
        self.reading_dataset = False

        # Names of all samples in the dataset
        self.sample_names = []

        # Start on first asid. # TODO: 0 is not necessarily a valid as_id
        self.all_asids = []
        self.current_asid = 0

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

        # Default filter options
        self.filters = self.get_default_filters()

        # Flags to indicate if current filters are valid
        self.valid_filters = {
            "psi_filter": 1,
            "included_counts_filter": 1,
            "excluded_counts_filter": 1,
            "rpkm_filter": 1
        }

        ###########################
        # Bind various keypresses #
        ###########################
        self.bind("<Escape>", self.close_window)
        self.bind("<Control-o>", lambda event=None: self.open_file())
        self.bind("<Control-s>", lambda event=None: self.save_file())
        self.bind("<Control-f>", lambda event=None: self.show_filter_dataset_window())
        self.bind("<Control-a>", lambda event=None: self.show_options())
        self.bind("<Left>", self.left_arrow_clicked)
        self.bind("<Right>", self.right_arrow_clicked)
        self.bind("<Up>", self.up_arrow_clicked)
        self.bind("<Down>", self.down_arrow_clicked)

        #################################
        # Draw graphical user interface #
        #################################

        # Prepare visual styles
        self.setup_ttk_styles()

        # Create a menu
        self.main_menu = self.create_menu()

        # Create main frame
        self.main_frame = self.create_main_frame()
        self.main_frame.columnconfigure(0, weight=1)
        self.main_frame.rowconfigure(0, weight=1)

        # Create left frame
        self.left_frame = self.create_left_frame()
        self.left_frame.rowconfigure(0, weight=1)
        self.left_frame.columnconfigure(1, weight=1)

        # Create sample frame
        self.sample_frame = self.create_sample_frame()
        self.create_sample_header()

        # Create exon frame
        self.exon_frame = self.create_exon_frame()

        # Create sidebar
        self.right_sidebar = self.create_right_sidebar()

        # Create load frame
        self.draw_animation = False  # Controls whether to show a loading animation
        self.load_frame, self.spinner_label, self.spinner_frames = self.create_load_frame()

        # Create statusbar at the bottom
        self.statusbar = self.create_statusbar()

        # A processor to handle I/O and system calls
        self.data_processor = TINDataProcessor(self, TAG_NO_TAG)

        # Paths to bam files
        # self.bam_paths = self.data_processor.get_bam_file_paths()

        ######################################
        # All set to handle user interaction #
        ######################################
        self.update_information()

    def create_main_frame(self):
        """
        Creates the main window of the application
        """
        main_frame = ttk.Frame(self, padding=0)
        main_frame.grid(column=0, row=0, sticky="NEWS")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        return main_frame

    def create_left_frame(self):
        """
        Left frame contains everything except the right sidebar.
        """
        left_frame = ttk.Frame(self.main_frame, padding=0)
        left_frame.grid(column=0, row=0, sticky="NEWS")
        return left_frame

    def create_sample_frame(self):
        """
        The sample frame is a part of the left frame and contains the sample header and sample info + RPKMs
        """
        # padding = (WEST, NORTH, EAST, SOUTH)
        sample_frame = ttk.Frame(self.left_frame, padding=(5, 0, 5, 0))
        sample_frame.grid(column=0, row=0, sticky="WNS")

        # Vertical separator between samples frame and exons frame
        sep = ttk.Separator(self.left_frame, orient=tk.VERTICAL)
        sep.grid(column=0, row=0, sticky="NSE")

        return sample_frame

    def create_sample_header(self):
        """
        The sample header is the top part of the sample frame and only displays a label saying "Sample".
        """
        sample_header = ttk.Frame(self.sample_frame, padding=1)  # The 1 padding is to compensate for the bold text in exon name increasing the exon name frame's height. The padding is to align the two frames better.
        sample_header.grid(column=0, row=0, sticky="NEW")
        header_label = ttk.Label(sample_header, text="Samples", font="TkDefaultFont 16", anchor="center")
        header_label.grid(column=0, row=0, sticky="NEWS")
        sample_header.columnconfigure(0, weight=1)

    def create_exon_frame(self):
        """
        The exon frame is part of the left frame and contains the exon names and the exon canvases
        """
        exon_frame = ttk.Frame(self.left_frame, padding=0, borderwidth=0, relief=tk.SOLID)
        exon_frame.grid(column=1, row=0, sticky="NEWS")

        return exon_frame

    def create_right_sidebar(self):
        """ Creates the sidebar to the right"""

        right_sidebar = ttk.Frame(self.main_frame, padding=5)
        right_sidebar.grid(column=1, row=0, sticky="ENS")

        # Information frame
        sidebar_information = ttk.Frame(right_sidebar)
        sidebar_information.grid(column=0, row=0, sticky="NEWS")
        right_sidebar.rowconfigure(0, weight=1)

        # Write header
        info_header = ttk.Label(sidebar_information, text="Information", font="TkDefaultFont 20", anchor=tk.CENTER)
        info_header.grid(column=0, row=0, columnspan=3, sticky="NEWS")

        # Define fonts for labels and text
        label_font = "TkDefaultFont 16 bold"
        text_font = "TkDefaultFont 16"

        # Keep track of which row we're working on
        current_row = 1

        # Alternative splicing ID
        asid_label = ttk.Label(sidebar_information, text="AS ID:", font=label_font)
        asid_label.grid(column=0, row=current_row, sticky="W")
        self.asid_text = ttk.Label(sidebar_information, text="<AS ID>", font=text_font)
        self.asid_text.grid(column=1, row=current_row, sticky="W", columnspan=2)
        current_row += 1

        # Gene symbol
        gene_label = ttk.Label(sidebar_information, text="Gene symbol:", font=label_font)
        gene_label.grid(column=0, row=current_row, sticky=tk.W)
        self.gene_text = ttk.Label(sidebar_information, text="<Gene symbol>", font=text_font)
        self.gene_text.grid(column=1, row=current_row, sticky=tk.W, columnspan=2)
        current_row += 1

        # Splice type
        splice_label = ttk.Label(sidebar_information, text="Splicing type:", font=label_font)
        splice_label.grid(column=0, row=current_row, sticky=tk.W)
        self.splice_type_text = ttk.Label(sidebar_information, text="<Splice Type>", font=text_font)
        self.splice_type_text.grid(column=1, row=current_row, sticky=tk.W, columnspan=2)
        current_row += 1

        # Exons
        exons_label = ttk.Label(sidebar_information, text="Exons:", font=label_font)
        exons_label.grid(column=0, row=current_row, sticky="W")
        self.exons_text = ttk.Label(sidebar_information, text="<Exons>", font=text_font, wraplength=250)
        self.exons_text.grid(column=1, row=current_row, sticky="W", columnspan=2)
        current_row += 1

        # Strand
        strand_label = ttk.Label(sidebar_information, text="Strand:", font=label_font)
        strand_label.grid(column=0, row=current_row, sticky="W")
        self.strand_text = ttk.Label(sidebar_information, text="<Strand>", font=text_font)
        self.strand_text.grid(column=1, row=current_row, sticky="W", columnspan=2)
        current_row += 1

        # Location / coordinates
        location_label = ttk.Label(sidebar_information, text="Location:", font=label_font)
        location_label.grid(column=0, row=current_row, sticky=tk.W)
        self.location_text = ttk.Label(sidebar_information, text="<Coordinates>", font=text_font)
        self.location_text.grid(column=1, row=current_row, sticky=tk.W)
        location_button = ttk.Button(sidebar_information, text="Copy", command=self.copy_coordinates_to_clipboard)
        location_button.grid(column=2, row=current_row, sticky="E")
        # location_text.bind("<Button-1>", location_clicked)
        #self.location_text.bind("<Enter>", lambda e: self.set_statusbar_text("Left click to copy coordinates to clipboard"))

        current_row += 1

        # TEST: Button for training decision tree
        decision_tree_training_button = ttk.Button(sidebar_information, text="Train decision tree", command=self.train_decision_tree)
        decision_tree_training_button.grid(column=0, row=current_row, sticky="W")
        current_row += 1

        next_untagged_button = ttk.Button(sidebar_information, text="Jump to next untagged event", command=self.next_untagged_event_button_clicked)
        next_untagged_button.grid(column=0, row=current_row, sticky="W")
        current_row += 1
        # END TEST

        # Buttons frame
        sidebar_buttons = ttk.Frame(right_sidebar, padding=5)
        sidebar_buttons.grid(column=0, row=1, sticky="EWS")

        # Previous button
        prev_button = ttk.Button(sidebar_buttons, text="< Previous", width=len("< Previous"), command=self.previous_button_clicked)
        prev_button.grid(column=0, row=0, sticky=tk.W)

        # Random button
        random_button = ttk.Button(sidebar_buttons, text="Random", command=self.random_button_clicked)
        random_button.grid(column=1, row=0)

        # Next button
        next_button = ttk.Button(sidebar_buttons, text="Next >", width=len("< Previous"), command=self.next_button_clicked)
        next_button.grid(column=2, row=0, sticky=tk.E)
        sidebar_buttons.columnconfigure(0, weight=1)
        #sidebar_buttons.columnconfigure(1, weight=1)
        sidebar_buttons.columnconfigure(2, weight=1)

        return right_sidebar

    def train_decision_tree(self):
        """
        Tells the dataprocessor to tell its tin_learner to train its decision tree.
        """
        tree_accuracy = self.data_processor.tin_learner.train_decision_tree(self.original_dataset)
        if not tree_accuracy:
            self.set_statusbar_text("ERROR: Decision tree could not be trained: Too few tagged events.")
        else:
            self.set_statusbar_text("Decision tree trained. Accuracy: %.2f" % tree_accuracy)

    def create_statusbar(self):
        """
        Creates the statusbar at the bottom of the main window
        """
        statusbar = ttk.Frame(self, borderwidth=1, relief=tk.SUNKEN)
        statusbar.grid(column=0, row=1, sticky="NEWS")

        # Text far left of status bar
        self.statusbar_text = ttk.Label(statusbar, text="This is the statusbar", font="TkDefaultFont")
        self.statusbar_text.grid(column=0, row=0, sticky="W")

        ttk.Separator(statusbar, orient=tk.VERTICAL).grid(column=1, row=0, sticky="NS")

        # Frame for tag counts
        tags_padding = (5, 0, 5, 0)
        tags_frame = ttk.Frame(statusbar, padding=tags_padding)
        tags_frame.grid(column=2, row=0, sticky="E")
        # Interesting count
        self.statusbar_text_interesting = ttk.Label(tags_frame, padding=tags_padding, text="0", font="TkDefaultFont", foreground=COLOR_INTERESTING)
        self.statusbar_text_interesting.grid(column=0, row=0, sticky="NEWS")
        self.statusbar_text_interesting.bind("<Motion>", lambda event=None: self.set_statusbar_text("Tagged interesting: %s" % self.statusbar_text_interesting["text"]))
        # Not interesting count
        self.statusbar_text_not_interesting = ttk.Label(tags_frame, padding=tags_padding, text="0", font="TkDefaultFont", foreground=COLOR_NOT_INTERESTING)
        self.statusbar_text_not_interesting.grid(column=1, row=0, sticky="NEWS")
        self.statusbar_text_not_interesting.bind("<Enter>", lambda event=None: self.set_statusbar_text("Tagged not interesting: %s" % self.statusbar_text_not_interesting["text"]))
        # Uncertain
        self.statusbar_text_uncertain = ttk.Label(tags_frame, padding=tags_padding, text="0", font="tkDefaultFont", foreground=COLOR_UNCERTAIN)
        self.statusbar_text_uncertain.grid(column=2, row=0, sticky="NEWS")
        self.statusbar_text_uncertain.bind("<Enter>", lambda event=None: self.set_statusbar_text("Tagged uncertain: %s" % self.statusbar_text_uncertain["text"]))

        ttk.Separator(statusbar, orient=tk.VERTICAL).grid(column=3, row=0, sticky="NS")

        # Tagged / total events text
        self.statusbar_text_progress = ttk.Label(statusbar, padding=tags_padding, text="0/0", font="TkDefaultFont")
        self.statusbar_text_progress.grid(column=4, row=0, sticky="E")
        self.statusbar_text_progress.bind("<Enter>", lambda event=None: self.set_statusbar_text("Events tagged: %s" % self.statusbar_text_progress["text"]))

        ttk.Separator(statusbar, orient=tk.VERTICAL).grid(column=5, row=0, sticky="NS")

        # Unsaved text (far right)
        self.statusbar_text_unsaved = ttk.Label(statusbar, padding=tags_padding, text="0", font="TkDefaultFont")
        self.statusbar_text_unsaved.grid(column=6, row=0, sticky="E")
        self.statusbar_text_unsaved.bind("<Enter>", lambda event=None: self.set_statusbar_text("Unsaved tags: %s" % self.statusbar_text_unsaved["text"]))

        statusbar.columnconfigure(0, weight=1)

        return statusbar

    def create_load_frame(self):
        """
        Creates an empty frame with "open file.." button.
        """

        # TODO: Align DB labels and entries.

        # Display prompts to load dataset
        load_frame = ttk.Frame(self.main_frame)
        load_frame.grid(column=0, row=0, sticky="NEWS")

        # Create a top frame for file and DB widgets
        top_frame = ttk.Frame(load_frame)
        top_frame.grid(column=0, row=0, sticky="NEWS")
        top_frame.columnconfigure(0, weight=1)
        top_frame.rowconfigure(0, weight=1)

        # And a bottom frame for the spinner
        bottom_frame = ttk.Frame(load_frame)
        bottom_frame.grid(column=0, row=1, sticky="NEWS")
        bottom_frame.columnconfigure(0, weight=1)

        # Divide even weight for top and bottom frames
        load_frame.rowconfigure(0, weight=1)
        load_frame.rowconfigure(1, weight=1)
        load_frame.columnconfigure(0, weight=1)

        #######################
        ###### Top frame ######
        #######################

        # Title
        title_label = ttk.Label(top_frame, text="No dataset currently loaded", font="tkDefaultFont 26")
        title_label.grid(column=0, row=0, sticky="N")
        # Frame for open file text and button
        file_frame = ttk.Frame(top_frame)
        file_frame.grid(column=0, row=1)
        # Load file label and button
        file_label = ttk.Label(file_frame, text="Load from file:")
        file_label.grid(column=0, row=0)
        file_button = ttk.Button(file_frame, text="Open file..", command=self.open_file)
        file_button.grid(column=1, row=0)
        # - OR - text
        or_text = ttk.Label(top_frame, text="- OR -", padding=(10, 100, 10, 100))
        or_text.grid(column=0, row=2)
        # Database frame
        db_labelframe = ttk.LabelFrame(top_frame, text="Load from database")
        db_labelframe.grid(column=0, row=3)
        # Left DB frame (for text lables)
        db_left_frame = ttk.Frame(db_labelframe)
        db_left_frame.grid(column=0, row=0)
        # Right DB frame (for text entries)
        db_right_frame = ttk.Frame(db_labelframe)
        db_right_frame.grid(column=1, row=0)
        # URL label
        url_label = ttk.Label(db_left_frame, text="Database URL:")
        url_label.grid(column=0, row=0, sticky="W")
        # User label
        user_label = ttk.Label(db_left_frame, text="Username:")
        user_label.grid(column=0, row=1, sticky="W")
        # Password label
        password_label = ttk.Label(db_left_frame, text="Password:")
        password_label.grid(column=0, row=2, sticky="W")
        # URL entry
        url_entry = ttk.Entry(db_right_frame)
        url_entry.grid(column=0, row=0)
        # User entry
        user_entry = ttk.Entry(db_right_frame)
        user_entry.grid(column=0, row=1)
        # Password entry
        password_entry = ttk.Entry(db_right_frame, show="*")
        password_entry.grid(column=0, row=2)
        # Name label
        name_label = ttk.Label(db_left_frame, text="Database name:")
        name_label.grid(column=0, row=3)
        # Name entry
        name_entry = ttk.Entry(db_right_frame)
        name_entry.grid(column=0, row=3)
        # Connect button
        connect_button = ttk.Button(
            db_labelframe,
            text="Connect",
            command=lambda: self.get_dataset_from_database(
                db_url=url_entry.get(),
                db_username=user_entry.get(),
                db_password=password_entry.get(),
                db_name=name_entry.get()
            )
        )
        connect_button.grid(column=0, row=4, columnspan=2, sticky="E")

        ##########################
        ###### Bottom frame ######
        ##########################
        # Label for the spinner
        spinner_label = ttk.Label(bottom_frame)
        spinner_label.grid(column=0, row=0)
        num_frames = 30
        spinner_frames = [tk.PhotoImage(file="spinner.gif", format="gif -index %i" % i) for i in range(num_frames)]

        return load_frame, spinner_label, spinner_frames

    def update_spinner_animation(self, index):
        """
        Handles the animation of a loading GIF, essentially just iterating frames in a GIF-file.
        """
        if not self.draw_animation:
            # Stop drawing this animation
            return
        frame = self.spinner_frames[index]
        index += 1
        if index >= 30:
            index = 0
        self.spinner_label.configure(image=frame)
        self.after(10, self.update_spinner_animation, index)

    def get_dataset_from_database(self, db_url, db_username, db_password, db_name):
        """
        Queries the SpliceSeq database to retrieve dataset.
        :param db_url: URL to database
        :param db_username: Username
        :param db_password: Password
        :param db_name: Name of database
        """

        # Spin loading animation while dataset is retrieved
        self.draw_animation = True
        self.after(0, self.update_spinner_animation, 0)

        # Get dataset from DB via data processor
        process_queue = Queue()
        self.reading_dataset = True
        query_db_process = Process(
            target=self.data_processor.get_dataset_from_database,
            args=(
                db_url,
                db_username,
                db_password,
                db_name,
                process_queue
            )
        )
        query_db_process.start()
        self.check_io_queue(process_queue)

    def set_statusbar_text(self, text):
        self.statusbar_text["text"] = text

    def copy_coordinates_to_clipboard(self):
        """
        Copies the genomic coordinates of the current event to the system clipboard.
        """
        # Get coordinates text
        coordinates = self.location_text["text"]

        if coordinates != "<Coordinates>":
            self.clipboard_clear()
            self.clipboard_append(coordinates)
            self.set_statusbar_text("Copied to clipboard: %s" % coordinates)
        else:
            self.set_statusbar_text("No coordinates to copy. Please load a dataset.")

    def close_window(self, event):
        # TODO: Clean up, prompt to save progress etc. before quitting.
        sys.exit()

    def create_menu(self):
        """
        Creates a menu on the top of the window.
        """
        main_menu = tk.Menu(self)
        self.config(menu=main_menu)

        sub_menu = tk.Menu(main_menu)
        main_menu.add_cascade(label="File", menu=sub_menu)
        sub_menu.add_command(label="Open file...", command=self.open_file, accelerator="Ctrl+O")
        sub_menu.add_command(label="Save", command=self.save_file, accelerator="Ctrl+S")
        sub_menu.add_separator()
        sub_menu.add_command(label="Options", command=self.show_options, accelerator="Ctrl+A")
        sub_menu.add_separator()
        sub_menu.add_command(label="Quit", accelerator="Escape", command=lambda: self.close_window(event=None))

        # Add dataset menu
        # TODO: Add an "information" dataset entry that shows the number of events, number of samples, and column names
        dataset_menu = tk.Menu(main_menu)
        main_menu.add_cascade(label="Filters", menu=dataset_menu)
        dataset_menu.add_command(label="Filter dataset", command=self.show_filter_dataset_window, accelerator="Ctrl+F")
        dataset_menu.add_separator()
        dataset_menu.add_command(label="Open filters..", command=self.read_dataset_filters)
        dataset_menu.add_command(label="Save current filters", command=self.save_dataset_filters)

        # Return main menu
        return main_menu

    def show_options(self):
        """
        Displays a window with options.
        """
        ############################################
        # Create a window to display everything in #
        ############################################
        options_window = tk.Toplevel()
        options_window.bind("<Escape>", lambda event=None: options_window.destroy())
        options_window.wm_title("Options")

        # A frame to hold the different options
        options_frame = ttk.Frame(options_window, padding=(20, 20, 20, 20))
        options_frame.grid(column=0, row=0, sticky="NEWS")
        options_frame.focus_force()  # Force focus to a widget in this window so that binds work

        # Add a label
        theme_label = ttk.Label(options_frame, text="Visual theme:", font="tkDefaultFont")
        theme_label.grid(column=0, row=0, sticky="NEWS")

        # Add an options menu with theme options
        theme_menu = ttk.OptionMenu(options_frame, self.current_theme, self.current_theme.get(), *self.available_themes)
        theme_menu.grid(column=1, row=0, sticky="NEWS")

    def change_theme(self, *args):
        """
        Called when a theme is changed in the options menu. Traces the self.current_theme variable. Changes theme.
        """
        self.style.theme_use(self.current_theme.get())

    def show_filter_dataset_window(self):

        # TODO: Handle overfloating horizontal space. Scrollbar or something.

        # Create a new window to display filters in
        window = tk.Toplevel()
        window.bind("<Escape>", lambda event=None: window.destroy())  # Close window on escape

        # Frame for filters
        filters_frame = ttk.Frame(window)
        filters_frame.grid(column=0, row=0, sticky="NEWS")
        window.rowconfigure(0, weight=1)

        ##########################
        ##### Normal filters #####
        ##########################
        current_row = 0
        normal_header = ttk.Label(filters_frame, text="Filters", font="TkDefaultFont 16")
        normal_filters_labelframe = ttk.LabelFrame(filters_frame, padding=(0, 10, 0, 10), labelanchor="n", labelwidget=normal_header)
        normal_filters_labelframe.grid(column=0, row=current_row, sticky="NEWS")
        filters_frame.rowconfigure(current_row, weight=1)
        current_row += 1
        normal_current_row = 0  # Keep track of which row we're at
        # Add filters to the labelframe
        for filter_name in ["psi", "included_counts", "excluded_counts", "rpkm", "avg_rpkm", "tot_reads", "occurrences"]:

            filter_label = ttk.Label(normal_filters_labelframe, text="%s:" % self.filters[filter_name][2])
            filter_label.grid(column=0, row=normal_current_row, sticky="NEWS")

            self.filters[filter_name][1].set(self.filters[filter_name][0])
            filter_entry = ttk.Entry(normal_filters_labelframe, textvariable=self.filters[filter_name][1])
            filter_entry.grid(column=1, row=normal_current_row, sticky="NEWS")

            normal_current_row += 1

        # Make sure labels and entries follow on window rescale
        normal_filters_labelframe.columnconfigure(0, weight=1)
        normal_filters_labelframe.columnconfigure(1, weight=1)

        ###############################
        ##### Splice type filters #####
        ###############################
        splice_header = ttk.Label(filters_frame, text="Splice type filters", font="TkDefaultFont 16")
        splice_type_labelframe = ttk.LabelFrame(filters_frame, padding=(0, 10, 0, 10), labelanchor="n", labelwidget=splice_header)
        splice_type_labelframe.grid(column=0, row=current_row, sticky="NEWS")
        filters_frame.rowconfigure(current_row, weight=1)
        current_row += 1
        current_splice_row = 0
        for splice_type in self.filters["splice_type"].keys():
            # Set current value to control variable
            self.filters["splice_type"][splice_type][1].set(self.filters["splice_type"][splice_type][0])

            splice_label = ttk.Label(splice_type_labelframe, text=self.filters["splice_type"][splice_type][2])
            splice_label.grid(column=0, row=current_splice_row, sticky="W")

            checkbutton = ttk.Checkbutton(
                splice_type_labelframe,
                variable=self.filters["splice_type"][splice_type][1],
                onvalue=True,
                offvalue=False
            )
            checkbutton.grid(column=1, row=current_splice_row, sticky="NEWS")
            current_splice_row += 1

        #############################
        ##### Event tag filters #####
        #############################
        tag_header = ttk.Label(filters_frame, text="Event tag filters", font="TkDefaultFont 16")
        tag_labelframe = ttk.LabelFrame(filters_frame, padding=(0, 10, 0, 10), labelanchor="n", labelwidget=tag_header)
        tag_labelframe.grid(column=0, row=current_row, sticky="NEWS")
        filters_frame.rowconfigure(current_row, weight=1)
        current_row += 1
        current_tag_row = 0
        for tag in self.filters["event_tag"].keys():
            # Set current value to control variable
            self.filters["event_tag"][tag][1].set(self.filters["event_tag"][tag][0])

            # Add label
            tag_label = ttk.Label(tag_labelframe, text=self.filters["event_tag"][tag][2])
            tag_label.grid(column=0, row=current_tag_row, sticky="W")

            checkbutton = ttk.Checkbutton(
                tag_labelframe,
                variable=self.filters["event_tag"][tag][1],
                onvalue=True,
                offvalue=False
            )
            checkbutton.grid(column=1, row=current_tag_row, sticky="W")
            current_tag_row += 1

        ############################
        ##### Advanced filters #####
        ############################
        avd_header = ttk.Label(filters_frame, text="Advanced filters", font="TkDefaultFont 16")
        advanced_labelframe = ttk.LabelFrame(filters_frame, padding=(0, 10, 0, 10), labelanchor="n", labelwidget=avd_header)
        advanced_labelframe.grid(column=0, row=current_row, sticky="NEWS")
        filters_frame.rowconfigure(current_row, weight=1)
        current_row += 1
        current_advanced_row = 0

        for adv_filter in [
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
        ]:
            # Set current value to control variable
            self.filters[adv_filter][1].set(self.filters[adv_filter][0])

            # Add label
            adv_label = ttk.Label(advanced_labelframe, text=self.filters[adv_filter][2])
            adv_label.grid(column=0, row=current_advanced_row, sticky="NEWS")

            # Add entry
            adv_entry = ttk.Entry(advanced_labelframe, textvariable=self.filters[adv_filter][1])
            adv_entry.grid(column=1, row=current_advanced_row, sticky="NEWS")
            current_advanced_row += 1

        # Distribute equal weight to both columns of labelframe
        advanced_labelframe.columnconfigure(0, weight=1)
        advanced_labelframe.columnconfigure(1, weight=1)

        # Frame for buttons
        buttons_frame = ttk.Frame(window)
        buttons_frame.grid(column=0, row=1, sticky="NEWS")
        buttons_frame.columnconfigure(0, weight=1)
        buttons_frame.columnconfigure(1, weight=1)
        # Cancel button
        cancel_button = ttk.Button(buttons_frame, text="Cancel", command=lambda: window.destroy())
        cancel_button.grid(column=0, row=0, sticky="W")
        # Apply button
        apply_button = ttk.Button(buttons_frame, text="Apply", command=self.apply_filters)
        apply_button.grid(column=1, row=0, sticky="E")
        # Bind apply button to enter key
        window.bind("<Return>", lambda event=None: apply_button.invoke())

        # Make sure window content follow on window rescale
        window.columnconfigure(0, weight=1)
        filters_frame.columnconfigure(0, weight=1)

    def get_default_filters(self):
        """
        Creates a default filter for the current original_dataset.
        """

        default_filters = {
            "psi": ["0.00", tk.StringVar(), "Min. PSI"],
            "included_counts": ["0", tk.StringVar(), "Min. main exon included counts"],
            "excluded_counts": ["0", tk.StringVar(), "Min. main exon excluded counts"],
            "rpkm": ["0.00", tk.StringVar(), "Min. gene RPKM"],
            "prev_exon_tot_reads": ["0", tk.StringVar(), "Min. upstream exon total reads"],  # Be aware that these do not apply to all splice types
            "prev_exon_rpkm": ["0.00", tk.StringVar(), "Min. upstream exon RPKM"],  # Be aware that these do not apply for all splice types
            "next_exon_tot_reads": ["0", tk.StringVar(), "Min. downstream exon total reads"],
            "next_exon_rpkm": ["0.00", tk.StringVar(), "Min. downstream exon RPKM"],
            "avg_rpkm": ["0.00", tk.StringVar(), "Min. main exon RPKM"],
            "tot_reads": ["0", tk.StringVar(), "Min. main exon total reads"],
            "avg_tot_reads": ["0.00", tk.StringVar(), "Min. average main exon total reads (all samples)"],
            "occurrences": ["0", tk.StringVar(), "Min. number of samples event is present in"],
            "prev_exon_max_rpkm": ["0.00", tk.StringVar(), "Min. upstream exon max RPKM"],
            "next_exon_max_rpkm": ["0.00", tk.StringVar(), "Min. downstream exon max RPKM"],
            "max_avg_rpkm": ["0.00", tk.StringVar(), "Min. main exon max RPKM"],
            "max_gene_rpkm": ["0.00", tk.StringVar(), "Min. gene max RPKM"],
            "max_psi": ["0.00", tk.StringVar(), "Min. main exon max PSI"],
            "percent_of_max_psi": ["0.00", tk.StringVar(), "Min. main exon PSI percent of max. PSI"],
            "percent_of_max_rpkm": ["0.00", tk.StringVar(), "Min. main exon RPKM percent of max in all samples"],
            "main_rpkm_to_upstream_rpkm_ratio": ["0.00", tk.StringVar(), "Min, main exon RPKM to upstream exon RPKM ratio"],
            "main_rpkm_to_downstream_rpkm_ratio": ["0.00", tk.StringVar(), "Min main exon RPKM to downstream exon RPKM ratio"],
            "sum_psi_all_samples": ["0.00", tk.StringVar(), "Min. total PSI in all samples"],
            "sum_psi_other_samples": ["0.00", tk.StringVar(), "Min. total PSI in other samples (excluding self)"],
            "mean_psi_other_samples": ["0.00", tk.StringVar(), "Min. average PSI in other sample (excluding self)"],
            "psi_diff_from_mean_other_samples": ["0.00", tk.StringVar(), "Min. difference in PSI from other samples"],
            "sum_rpkm_all_samples": ["0.00", tk.StringVar(), "Min. sum main exon RPKM in all samples"],
            "sum_rpkm_other_samples": ["0.00", tk.StringVar(), "Min. sum main exon RPKM in other samples"],
            "mean_rpkm_other_samples": ["0.00", tk.StringVar(), "Min. average main exon RPKM in other samples"],
            "rpkm_percentage_of_mean_other_samples": ["0.00", tk.StringVar(), "Min. main exon RPKM percentage of average in other samples"],
            "splice_type": {
                "AA": [True, tk.BooleanVar(), "Alternate Acceptor"],
                "AD": [True, tk.BooleanVar(), "Alternate Donor"],
                "ES": [True, tk.BooleanVar(), "Exon Skipping"],
                "RI": [True, tk.BooleanVar(), "Retained Intron"],
                "AP": [True, tk.BooleanVar(), "Alternate Promotor"],
                "AT": [True, tk.BooleanVar(), "Alternate Terminator"],
                "ME": [True, tk.BooleanVar(), "Mutually Exclusive Exon"],
            },
            "event_tag": {
                "interesting": [True, tk.BooleanVar(), "Interesting event"],
                "not_interesting": [True, tk.BooleanVar(), "Not interesting event"],
                "uncertain": [True, tk.BooleanVar(), "Uncertain"],
                "no_tag": [True, tk.BooleanVar(), "No tag"]
            }
        }

        print "Merge new filter approach onto master branch"
        return default_filters

    def apply_filters(self):
        """
        Sanitize filter values and pass them on to the dataprocessor to be applied to the dataset.
        """

        # Check if dataset is present
        try:
            if not self.original_dataset:
                tkMessageBox.showerror("Filter error", "Cannot apply filters, as no dataset is loaded. Please load one.")
                return
        except ValueError:
            # It's a valid dataframe, continue
            pass

        int_fields = ["included_counts", "excluded_counts", "tot_reads", "occurrences"]
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

        # Sanitize int fields
        int_errors = []
        for i in int_fields:
            try:
                int_value = int(self.filters[i][1].get())
                self.filters[i][0] = self.filters[i][1].get()
            except ValueError:
                print "Error sanitizing int field %s: <%s>" % (i, self.filters[i][1].get())
                int_errors.append(self.filters[i][2])

        # Sanitize float fields
        float_errors = []
        for f in float_fields:
            try:
                float_value = float(self.filters[f][1].get())
                self.filters[f][0] = self.filters[f][1].get()
            except ValueError:
                print "Error sanitizing float field %s: <%s>" % (f, self.filters[f][1].get())
                float_errors.append(self.filters[f][2])

        # Check for errors
        error_message = ""
        if len(int_errors) > 0:
            error_message += "ERROR\n\nThe following fields must be integers (whole numbers, not floats or text. E.g. 1, not 1.00 or 'one'):\n\n"
            error_message += ", ".join(["<" + e + ">" for e in int_errors])

        if len(float_errors) > 0:
            error_message += "\n\nAnd the following fields must be floats:\n\n"
            error_message += ", ".join(["<" + f + ">" for f in float_errors])

        if len(int_errors) + len(float_errors) > 0:
            tkMessageBox.showerror("Filter error", error_message)
            return

        # Store selected values for splice types (not sure if this is necessary)
        for splice_type in self.filters["splice_type"].keys():
            self.filters["splice_type"][splice_type][0] = self.filters["splice_type"][splice_type][1].get()

        # And the same for event tags (again, this might be unnecessary)
        for tag in self.filters["event_tag"].keys():
            self.filters["event_tag"][tag][0] = self.filters["event_tag"][tag][1].get()

        # Fields are sanitized
        filtered_dataset = self.data_processor.filter_dataset(self.original_dataset, self.filters)

        # Get dataset
        if filtered_dataset.empty:
            tkMessageBox.showinfo("Filter error", "No events matching the current filter.")
            return

        # New dataset is fine, update UI
        self.dataset = filtered_dataset
        self.all_asids = sorted(list(self.dataset["as_id"].unique()))
        self.current_asid = self.all_asids[0]
        self.update_information()

    def save_dataset_filters(self):
        """
        Calls the TINDataProcessor to save the current dataset filters to disk, at the param filepath.
        """

        # Open file chooser dialog
        filepath = asksaveasfilename(title="Save filters", defaultextension=".json", filetypes=[("JSON", "*.json")])

        # Check if file was chosen
        if len(filepath) > 0:
            # Convert filters to only contain native datatypes before saving to disk
            native_filters = self.convert_filters_to_native_types()
            success = self.data_processor.save_filters(native_filters, filepath)
            if not success:
                print "Error: Unable to save filters."
                self.set_statusbar_text("Something went wrong when saving filters.")
            else:
                self.set_statusbar_text("Filters saved to file: %s" % filepath)

    def read_dataset_filters(self):
        """
        Reads datasets filters from file, then calls apply_filters() to update the UI.
        """

        # Open file chooser dialog
        filepath = askopenfilename(title="Open filters", defaultextension=".json", filetypes=[("JSON", "*.json")])

        # Check if file was chosen
        if len(filepath) > 0:
            # Read filters form file
            raw_filters = self.data_processor.read_filters(filepath)

            # Check if everything went OK (i.e. that None was not returned)
            if raw_filters:
                # Convert native datatypes to tk-specific datatypes
                converted_filters = self.convert_filters_to_tk_specific_datatypes(raw_filters)
                self.filters = converted_filters
                # Apply current filters
                self.apply_filters()
            else:
                print "Error when loading filters from file."
                self.set_statusbar_text("Filters not loaded: Something went wrong when loading filters from file.")

    def center_window(self, window):
        """
        Centers a window in the screen.
        """
        window.update_idletasks()
        width = window.winfo_screenwidth()
        height = window.winfo_screenheight()
        size = tuple(int(x) for x in window.geometry().split("+")[0].split("x"))
        x_coord = width/2 - size[0]/2
        y_coord = height/2 - size[1]/2

        window.geometry("%dx%d+%d+%d" % (size + (x_coord, y_coord)))

    def setup_ttk_styles(self):
        """
        Prepares ttk.Style() attributes for buttons, frames, etc.
        """
        style = ttk.Style()
        style.configure(STYLE_BUTTON_INTERESTING_OFF, foreground=COLOR_DARKBLUE)
        style.configure(STYLE_BUTTON_NOT_INTERESTING_OFF, foreground=COLOR_DARKBLUE)
        style.configure(STYLE_BUTTON_UNCERTAIN_OFF, foreground=COLOR_DARKBLUE, font="tkDefaultFont 16 bold")
        style.configure(STYLE_BUTTON_INTERESTING_ON, foreground=COLOR_INTERESTING)
        style.configure(STYLE_BUTTON_NOT_INTERESTING_ON, foreground=COLOR_NOT_INTERESTING)
        style.configure(STYLE_BUTTON_UNCERTAIN_ON, foreground=COLOR_UNCERTAIN, font="tkDefaultFont 16 bold")
        style.map(
            STYLE_BUTTON_INTERESTING_OFF,
            foreground=[("active", COLOR_DARKBLUE), ("disabled", COLOR_DARKWHITE)]
        )
        style.map(
            STYLE_BUTTON_NOT_INTERESTING_OFF,
            foreground=[("active", COLOR_DARKBLUE), ("disabled", COLOR_DARKWHITE)]
        )
        style.map(
            STYLE_BUTTON_UNCERTAIN_OFF,
            foreground=[("active", COLOR_DARKBLUE), ("disabled", COLOR_DARKWHITE)]
        )
        style.map(
            STYLE_BUTTON_INTERESTING_ON,
            foreground=[("active", COLOR_INTERESTING)]
        )
        style.map(
            STYLE_BUTTON_NOT_INTERESTING_ON,
            foreground=[("active", COLOR_NOT_INTERESTING)]
        )
        style.map(
            STYLE_BUTTON_UNCERTAIN_ON,
            foreground=[("active", COLOR_UNCERTAIN)]
        )

        style.configure("Graytext.TLabel", foreground=COLOR_DARKWHITE)

    def open_file(self):
        """
        Presents a filedialog for the user to choose a file. If a file is selected, the dataset is read
        into a Pandas.DataFrame
        """
        # Open file chooser dialog
        filename = askopenfilename()
        self.update()  # Without this, the file dialog remains present while the file is being read

        # Check if file was chosen
        if len(filename) > 0:
            self.read_dataset(filename)

    def read_dataset(self, filepath):
        """
        Calls the data processor to read file by path.
        """
        # Draw load animation while file is being read
        self.set_statusbar_text("Loading dataset..")
        self.draw_animation = True
        self.after(0, self.update_spinner_animation, 0)

        # Read dataset in a background process
        process_queue = Queue()
        self.reading_dataset = True  # True means we should keep checking for a dataset being read
        readfile_process = Process(target=self.data_processor.load_dataset, args=(filepath, process_queue))
        readfile_process.start()
        # Listen for changes to the queue, i.e. wether or not dataset has been read
        self.check_io_queue(process_queue)

    def check_io_queue(self, queue):
        """
        Checks whether dataset is done loading.
        :param queue: The Queue used for reading in the dataset
        :return:
        """
        # TODO: Handle empty dataset / non-compliant data formats
        try:
            self.original_dataset = queue.get_nowait()
            try:
                if not self.original_dataset:
                    # Reading functions returning False means something went wrong.
                    self.reading_dataset = False
                    self.draw_animation = False
                    self.set_statusbar_text("IO ERROR: something went wrong when reading dataset.")
                    return
            except ValueError:
                # It's a dataframe
                # TODO: Check if df is empty
                # For now, do nothing
                pass

            # Result is an actual dataset
            self.all_asids = sorted(list(self.original_dataset["as_id"].unique()))
            self.dataset = self.original_dataset.copy()
            self.reading_dataset = False
            # Find and store unique sample names
            self.sample_names = list(self.dataset["name"].unique())
            # Default to the first as_id in the file
            self.current_asid = self.all_asids[0]
            self.draw_animation = False

            # Handle data
            self.update_information()
        except Empty:
            # Do nothing, just wait for data to be read
            pass

        if self.reading_dataset:
            self.after(500, self.check_io_queue, queue)

    def update_information(self):
        """
        Handles top-level stuff for each event (row):
            - Collects row data
            - Orchestrates UI update
        """
        print ""
        print "Tagger: update_information()"

        # Set waiting cursor
        #self.config(cursor="watch")

        # TEST: Right-click on uncertain buttons
        self.all_uninteresting_buttons = []
        # END TEST

        # If no dataset is loaded, prompt user to load one
        if self.is_dataset_loaded():
            # TODO: Does this have to be called for every as_id? (The hide_load_frame, that is)
            # Dataset is loaded, hide load frame and display sidebars
            self.hide_load_frame()
        else:
            # No dataset loaded, hide sidebars and display load frame
            self.show_load_frame()
            self.set_statusbar_text("No dataset loaded.")
            # Dodge the open-file-dialog when testing/debugging # TODO: Remove this
            if self.testing:
                #self.read_dataset("/Users/jonas/Dropbox/phd/code/tin_tagger/datasets/mikes_query_with_tags.tsv")
                pass
            return

        # Display which event we're at in the statusbar
        self.set_statusbar_text("Splicing event %d/%d" % (self.all_asids.index(self.current_asid) + 1, len(self.all_asids)))

        # Update tagging progress information in statusbar
        self.update_tag_information()

        # Get data for this row
        data = self.data_processor.get_row_data(self.current_asid, self.original_dataset, self.sample_names, self.testing)

        # Populate the sidebar with general information
        self.asid_text["text"] = data["as_id"]
        self.location_text["text"] = data["coords"]
        self.splice_type_text["text"] = "%s (%s)" % (self.splice_type_map[data["splice_type"]], data["splice_type"])
        self.exons_text["text"] = data["exons"]
        self.gene_text["text"] = data["gene_symbol"]
        self.strand_text["text"] = data["strand"]

        # Clear all canvases before drawing new ones
        self.clear_all_canvases()

        # Also clear exon frame before creating new buttons
        self.exon_frame.destroy()
        self.exon_frame = self.create_exon_frame()

        # Draw events
        splice_type = data["splice_type"]
        if splice_type == "AT":
            self.draw_alternative_terminator_event(data)
        elif splice_type == "ES":
            self.draw_exon_skipping_event(data)
        elif splice_type == "AD":
            self.draw_alternative_donor_events(data)
        elif splice_type == "RI":
            self.draw_retained_intron_event(data)
        elif splice_type == "AA":
            self.draw_alternative_acceptor_events(data)
        elif splice_type == "AP":
            self.draw_alternative_promotor_event(data)
        elif splice_type == "ME":
            self.draw_mutually_exclusive_exons_event(data)

        # Reset cursor now that we're done with loading everything
        #self.config(cursor="")

    def is_dataset_loaded(self):
        """
        Checks whether there's a dataset currently loaded and returns answer as a boolean
        """
        try:
            if len(self.dataset) > 0:
                return True
            else:
                self.set_statusbar_text("Current dataset is empty.")
                return False
        except TypeError:
            return False

    def show_load_frame(self):
        """
        Shows a blank UI with a prompt to open a file, also hides the sidebars
        """
        self.right_sidebar.grid_remove()
        self.main_menu.entryconfig("Filters", state="disabled")  # Disable filters menu
        self.left_frame.grid_remove()
        self.load_frame.grid()

    def hide_load_frame(self):
        """
        Hides the open-file prompt and displays the sidebar
        """
        self.load_frame.destroy()
        self.right_sidebar.grid()
        self.left_frame.grid()
        self.main_menu.entryconfig("Filters", state="normal")  # Enable filters menu

    def update_tag_information(self):
        """
        Reads the number of tagged events in total and populates the statusbar with text about how many
        positive, negative, and neutral tags are in the dataset.
        """

        # Find events that are tagged as interesting
        positive_tags = len(self.original_dataset.loc[self.original_dataset["event_tag"] == TAG_INTERESTING])
        self.statusbar_text_interesting["text"] = "%d" % positive_tags

        # Find events that are tagged as not interesting
        negative_tags = len(self.original_dataset.loc[self.original_dataset["event_tag"] == TAG_NOT_INTERESTING])
        self.statusbar_text_not_interesting["text"] = "%d" % negative_tags

        # Find events that are tagged as uncertain
        neutral_tags = len(self.original_dataset.loc[self.original_dataset["event_tag"] == TAG_UNCERTAIN])
        self.statusbar_text_uncertain["text"] = "%d" % neutral_tags

        # Find how many events are tagged, in total
        total_events = len(self.original_dataset)
        total_tags = positive_tags + negative_tags + neutral_tags
        self.statusbar_text_progress["text"] = "%d/%d" % (total_tags, total_events)

    def save_file(self):
        print "Bleep, blop, saving file."
        filepath = asksaveasfilename(title="Save filters", defaultextension=".tsv")

        if not filepath:
            print "Unable to save to file: Filepath is None"
            return
        else:
            print "Saving file to %s" % filepath
            self.original_dataset.to_csv(filepath, sep="\t", index=False)

    def next_untagged_event_button_clicked(self):
        """
        Find and display next untagged event.
        """
        # Note: self.dataset is correct, because we want to respect current filter
        untagged_asid = self.data_processor.get_next_untagged_asid(self.current_asid, self.dataset)

        if untagged_asid == self.current_asid:
            # No untagged event found
            self.set_statusbar_text("ERROR: No untagged events found.")
            return
        else:
            self.current_asid = untagged_asid
            self.update_information()

    def left_arrow_clicked(self, event):
        """
        Handles when left-arrowkey is clicked. Essentially just passes through to
        previous_button_clicked.
        """
        self.previous_button_clicked()

    def right_arrow_clicked(self, event):
        """
        Handles when right-arrowkey is clicked. Essentially just passes through to
        next_button_clicked.
        """
        self.next_button_clicked()

    def up_arrow_clicked(self, event):
        """
        Called when the up-arrow is clicked. Increases canvas text size.
        """
        self.canvas_text_size += 1
        self.canvas_font = ("tkDefaultFont", self.canvas_text_size)
        print "UP! Current size:", self.canvas_text_size
        # TODO: HERE, redraw canvas?

    def down_arrow_clicked(self, event):
        """
        Called when the down-arrow is clicked. Decreases canvas text size.
        """
        self.canvas_text_size -= 1
        self.canvas_font = ("tkDefaultFont", self.canvas_text_size)
        print "DOWN! Current size:", self.canvas_text_size
        # TODO: HERE, redraw canvas?

    def next_button_clicked(self):
        """
        Handles next-button presses: Update as_id and initiate reading of new row.
        """
        next_asid_index = self.all_asids.index(self.current_asid) + 1

        if next_asid_index >= len(self.all_asids):
            print "Woops, no more rows (reached end of dataset)"
            next_asid_index -= 1
            self.set_statusbar_text("Reached end of dataset.")
        else:
            self.current_asid = self.all_asids[next_asid_index]
            self.update_information()

    def previous_button_clicked(self):
        """
        Handles previous-button presses: Update as_id and initiate reading of new row.
        """

        previous_asid_index = self.all_asids.index(self.current_asid) - 1
        if previous_asid_index < 0:
            print "Wops, no more rows (reached start of dataset)"
            self.set_statusbar_text("Reached begnning of dataset.")
        else:
            self.current_asid = self.all_asids[previous_asid_index]
            self.update_information()

    def random_button_clicked(self):
        """
        Handles random-button presses: Pick and display a random event.
        """

        random_asid = random.choice(self.all_asids)
        if random_asid == self.current_asid:
            # If the randomly chosen event is the same as it was before, just do nothing instead of loading everything
            # again. Also, be careful with the while-loops; it may be only 1 unique as_id in the dataset.
            return

        self.current_asid = random_asid
        self.update_information()

    def print_test(self, text):
        print "TEXT:", text
        self.set_statusbar_text(text)

    def clear_all_canvases(self):
        """
        Clear and delete canvases before drawing new ones.
        """

        # Clear all objects on canvas
        for row_canvas in self.canvases:
            row_canvas.delete("all")

        # Empty list of canvases
        self.canvases = []

    def populate_samples_frame(self, data):
        """
        Takes a dataset containing information about all samples and populates the left samples bar.
        """

        # Lambda to prevent lazy-interpretation of RPKMs (for on-hover display of RPKM)
        store_rpkm_text = lambda x: (lambda p: self.print_test(x))

        # Keep track of which row we're at
        row_number = 1  # Row number 0 is the sample header
        for sample_name in sorted(data["samples"].keys(), key=natural_sort_key):
            # Get sample information
            sample_data = data["samples"][sample_name]
            gene_rpkm = sample_data["gene_rpkm"]
            try:
                gene_rpkm_percent_of_max = (float(gene_rpkm) / float(sample_data["max_gene_rpkm"])) * 100
                try:
                    int_percent = int(gene_rpkm_percent_of_max)
                except ValueError:
                    print "Percent of max:", gene_rpkm_percent_of_max
                    gene_rpkm_percent_of_max = 0
            except ZeroDivisionError:
                gene_rpkm_percent_of_max = 0
            except ValueError:
                gene_rpkm_percent_of_max = 0

            # Create frame for this sample
            frame_for_sample = ttk.Frame(self.sample_frame, padding=5)
            frame_for_sample.grid(column=0, row=row_number, sticky="NEWS")
            self.sample_frame.rowconfigure(row_number, weight=1)

            # Display RPKM info in statusbar on hover
            display_text = "%s gene RPKM: %s/%s (%.0f%%)" % (
                sample_name, str(gene_rpkm), str(sample_data["max_gene_rpkm"]), gene_rpkm_percent_of_max)
            frame_for_sample.bind("<Enter>", store_rpkm_text(display_text))
            rpkm_color = COLOR_RED

            # Create text for sample name
            text_style = "TLabel"
            if not sample_data["is_reported"]:
                text_style = "Graytext.TLabel"
                rpkm_color = COLOR_DARKWHITE
            sample_text = ttk.Label(frame_for_sample, text=sample_name, font="TkDefaultFont", anchor=tk.CENTER, style=text_style)
            sample_text.grid(column=0, row=0, sticky="NEWS")
            frame_for_sample.columnconfigure(0, weight=9)

            # Create RPKM-indicator to the right
            rpkm_frame = ttk.Frame(frame_for_sample, borderwidth=1, relief=tk.SUNKEN)
            rpkm_frame.grid(column=1, row=0, sticky="ENS")

            frame_for_sample.columnconfigure(1, weight=1)
            frame_for_sample.rowconfigure(0, weight=1)
            # Fill remaining (what's left of 100% after this RPKM)
            remainder_frame = ttk.Frame(rpkm_frame, width=10)
            if gene_rpkm_percent_of_max == 100:
                # Hack to make it look like it's actually 100% and not just 99%
                remainder_frame = tk.Frame(rpkm_frame, width=10, background=rpkm_color)
            remainder_frame.grid(column=0, row=0, sticky="NEWS")
            # Fill for RPKM
            fill_frame = tk.Frame(rpkm_frame, bg=rpkm_color)
            fill_frame.grid(column=0, row=1, sticky="NEWS")
            rpkm_frame.rowconfigure(0, weight=100 - int(gene_rpkm_percent_of_max))
            rpkm_frame.rowconfigure(1, weight=int(gene_rpkm_percent_of_max))

            row_number += 1

            # TEST: Draw separator
            #sep = tk.Frame(self.sample_frame, bg=COLOR_DARKWHITE, height=2)
            #sep.grid(row=row_number, column=0, sticky="NEWS")
            #row_number+=1
            # For this to look nice, the sample_frame needs a 1-wide padding at the bottom.
            # END TEST

    def tag_button_clicked(self, sample_name, as_id, new_tag, up_button, down_button, uncertain_button):
        """
        Triggered when a tag button is clicked. Updates button styles and stores the new tag in the dataset.
        """

        # Get current tag for this sample and as_id
        current_tag = self.data_processor.get_tag_by_sample_name_and_as_id(sample_name, as_id, self.original_dataset)

        # This is the tag to be changed in the dataset
        set_tag = new_tag

        # If the current tag is the same as it was, then we'll "toggle" the tag back to not being tagged.
        if current_tag == new_tag:
            set_tag = TAG_NO_TAG

        # Update the (original) dataset
        self.data_processor.set_tag_by_sample_name_and_as_id(set_tag, sample_name, as_id, self.original_dataset)

        # Clear all button styles
        up_button.configure(style=STYLE_BUTTON_INTERESTING_OFF)
        down_button.configure(style=STYLE_BUTTON_NOT_INTERESTING_OFF)
        uncertain_button.configure(style=STYLE_BUTTON_UNCERTAIN_OFF)

        # If it's a toggle (the old and the new are the same), we don't activate any buttons.
        if current_tag != new_tag:
            if new_tag == TAG_INTERESTING:
                up_button.configure(style=STYLE_BUTTON_INTERESTING_ON)
            elif new_tag == TAG_NOT_INTERESTING:
                down_button.configure(style=STYLE_BUTTON_NOT_INTERESTING_ON)
            elif new_tag == TAG_UNCERTAIN:
                uncertain_button.configure(style=STYLE_BUTTON_UNCERTAIN_ON)

        # Finally, update tag information in status-bar
        self.update_tag_information()

    def add_tagging_buttons(self, row_number, sample_name, is_reported, sample_tag, as_id):
        """
        Add buttons
        """
        #########################
        # Setup tagging buttons #
        #########################
        # TEST: Algo tag
        test = random.randint(0, 2)
        algo_tag_color = COLOR_INTERESTING
        if test == 1:
            algo_tag_color = COLOR_NOT_INTERESTING
        if test == 0:
            algo_tag_color = COLOR_UNCERTAIN

        if not is_reported:
            algo_tag_color = COLOR_DARKWHITE
        # END TEST: Algo tag

        # Create a frame in which to store algorithm tag indicator + frame for buttons
        tagging_container = ttk.Frame(self.exon_frame)
        algo_tag_frame_column = 1
        button_frame_column = 0
        algo_tag_indicator_width = 5
        # Make the button frame expand to fill the remaining space
        tagging_container.grid(row=row_number, column=1, sticky="NEWS")
        # Make both the algo-tag indicator and the button frame fill vertical space
        tagging_container.columnconfigure(button_frame_column, weight=1)
        tagging_container.rowconfigure(0, weight=1)

        # Create a frame for algorithm tag indicator
        algo_tag_frame = tk.Frame(tagging_container, bg=algo_tag_color, width=algo_tag_indicator_width)
        algo_tag_frame.grid(row=0, column=algo_tag_frame_column, sticky="NEWS")

        button_frame = ttk.Frame(tagging_container)
        button_frame.grid(row=0, column=button_frame_column, sticky="NEWS")
        # Tag interesting button
        up_button_style = STYLE_BUTTON_INTERESTING_ON if sample_tag == TAG_INTERESTING else STYLE_BUTTON_INTERESTING_OFF
        up_button = ttk.Button(
            button_frame,
            text=u"\u25B2", style=up_button_style,
        )
        up_button.grid(column=0, row=0, sticky="NEWS")
        # Tag not interesting button
        down_button_style = STYLE_BUTTON_NOT_INTERESTING_ON if sample_tag == TAG_NOT_INTERESTING else STYLE_BUTTON_NOT_INTERESTING_OFF
        down_button = ttk.Button(
            button_frame,
            text=u"\u25BC",
            style=down_button_style,
        )
        down_button.grid(column=0, row=1, sticky="NEWS")
        # Tag uncertain button
        uncertain_button_style = STYLE_BUTTON_UNCERTAIN_ON if sample_tag == TAG_UNCERTAIN else STYLE_BUTTON_UNCERTAIN_OFF
        uncertain_button = ttk.Button(
            button_frame,
            text="?",
            style=uncertain_button_style,
        )
        uncertain_button.grid(column=0, row=2, sticky="NEWS")

        # Set callback functions for all buttons
        up_button.config(command=lambda name=sample_name, splice_id=as_id, new_tag=TAG_INTERESTING, button_up=up_button, button_down=down_button, button_uncertain=uncertain_button: self.tag_button_clicked(name, splice_id, new_tag, button_up, button_down, button_uncertain))
        down_button.config(command=lambda name=sample_name, splice_id=as_id, new_tag=TAG_NOT_INTERESTING, button_up=up_button, button_down=down_button, button_uncertain=uncertain_button: self.tag_button_clicked(name, splice_id, new_tag, button_up, button_down, button_uncertain))
        uncertain_button.config(command=lambda name=sample_name, splice_id=as_id, new_tag=TAG_UNCERTAIN, button_up=up_button, button_down=down_button, button_uncertain=uncertain_button: self.tag_button_clicked(name, splice_id, new_tag, button_up, button_down, button_uncertain))

        # TEST
        #uncertain_button.bind("<Button-3>", self.right_click)
        down_button.bind("<Button-2>", self.right_click)
        self.all_uninteresting_buttons.append(down_button)
        #uncertain_button.bind("<Button-1>", self.right_click)
        # END TEST

        # Disable buttons if the event is not reported for this sample
        if not is_reported:
            up_button.state(["disabled"])
            down_button.state(["disabled"])
            uncertain_button.state(["disabled"])

        button_frame.rowconfigure(0, weight=1)
        button_frame.rowconfigure(1, weight=1)
        button_frame.rowconfigure(2, weight=1)
        # END tagging buttons

    def right_click(self, event):
        print "RIGHT CLICK"
        """
        # THIS WORKS, BUT IS EXTREMELY HACKY (AND SLOW)
        # Mark all as uncertain
        for child in self.exon_frame.winfo_children():
            for subchild in child.winfo_children():
                for subsubchild in subchild.winfo_children():
                    print "Class:", subsubchild.winfo_class
                    #if subsubchild.winfo_class == "Button":
                    if isinstance(subsubchild, ttk.Button):
                        print "It's a button"
                        if subsubchild.cget("text") == "?":
                            subsubchild.invoke()

        """

        for uninteresting_button in self.all_uninteresting_buttons:
            uninteresting_button.invoke()

    def draw_exon_skipping_event(self, data):
        """
        Draw exon skipping events
        """
        ######################################
        # Fill in samples information column #
        ######################################
        self.populate_samples_frame(data)

        # Keep track of row number
        row_number = 0

        as_id = data["as_id"]

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        width_per_exon_container = canvas_width / 3
        exon_width = 60
        exon_height = 80
        exon_start_y = (canvas_height - exon_height) / 2
        top_canvas_height = 20

        ################################
        # Draw exon names in top frame #
        ################################
        canvas_background = COLOR_WHITE
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Add upstream exon name
        upstream_exon_name = data["prev_exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."

        # Add exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."

        # Add downstream exon name
        downstream_exon_name = data["next_exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + ".."

        ##############
        # Draw exons #
        ##############
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        prev_exon_id = data["prev_exon_id"]
        next_exon_id = data["next_exon_id"]

        # Get rpkm/tot_reads for flanking exons
        # flanking_exons_data = self.data_processor.get_flanking_exons_rpkm_by_exon_ids(sample_names_sorted, prev_exon_id, next_exon_id)
        #main_exon_rpkm_data = self.data_processor.get_main_exon_rpkm_by_asid(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        # Iterate samples
        for sample_name in sample_names_sorted:
            sample_data = data["samples"][sample_name]

            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]
            # TODO: Update colors depending on tag
            try:
                decision_tree_prediction = sample_data["decision_tree_prediction"]
                print "TINTagger got prediction <%d> for sample <%s>" % (decision_tree_prediction, sample_name)
            except KeyError:
                # Decision tree tag is not set yet. Default to no tag
                decision_tree_prediction = TAG_NO_TAG

            # Setup colors
            canvas_background = "white"
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE

            # If event is not reported, draw monochrome canvas
            if not is_reported:
                #canvas_background = COLOR_DARKWHITE
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY

            # Initialize canvas and grid to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            # Keep track of canvases used
            self.canvases.append(row_canvas)

            # Setup tagging buttons #
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)

            # Set even weight for every row in the exon frame
            self.exon_frame.rowconfigure(row_number, weight=1)

            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ######################
            # Draw upstream exon #
            ######################
            #upstream_exon_rpkm = flanking_exons_data[sample_name][prev_exon_id]["rpkm"]
            #upstream_exon_max_rpkm = flanking_exons_data[sample_name][prev_exon_id]["max_rpkm"]
            upstream_exon_rpkm = sample_data["prev_exon_rpkm"]
            upstream_exon_max_rpkm = sample_data["prev_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(upstream_exon_rpkm) / float(upstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            if np.isnan(percent_of_max_rpkm):
                percent_of_max_rpkm = 0

            upstream_exon_start_x = (width_per_exon_container - exon_width) / 2

            # Draw exon background
            row_canvas.create_rectangle(upstream_exon_start_x, exon_start_y, upstream_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(upstream_exon_start_x, fill_start_y, upstream_exon_start_x + exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)

            # Draw rpkm text. NOTE: Text colors are controlled in ResizingCanvas.on_resize()
            upstream_exon_text_start_x = upstream_exon_start_x + (exon_width / 2)
            upstream_exon_text_start_y = exon_height - 10
            upstream_textshadow = row_canvas.create_text(upstream_exon_text_start_x + 1, upstream_exon_text_start_y + 1, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            upstream_text = row_canvas.create_text(upstream_exon_text_start_x, upstream_exon_text_start_y, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)
            # TEST: Draw a box around the RPKM values
            #upstream_bbox = row_canvas.bbox(upstream_text)
            #rpkm_box = row_canvas.create_rectangle(upstream_bbox, outline="black", fill="white")
            #rpkm_box = row_canvas.create_rectangle(text_start_x - 5, text_start_y - 10, text_start_x+5, text_start_y+10, outline="black", fill="white")
            # Raise text above rpkm-box
            #row_canvas.tag_raise(upstream_textshadow, rpkm_box)
            #row_canvas.tag_raise(upstream_text, rpkm_box)
            # END TEST

            ##################
            # Draw main exon #
            ##################
            main_exon_start_x = upstream_exon_start_x + width_per_exon_container  # Not sure if this is correct
            main_exon_text_start_x = main_exon_start_x + (exon_width / 2)

            # Get RPKM values
            #sample_rpkm_data = main_exon_rpkm_data[sample_name]
            #combined_rpkm = sample_rpkm_data["combined_rpkm"]
            #max_combined_rpkm = sample_rpkm_data["max_combined_rpkm"]
            avg_rpkm = sample_data["avg_rpkm"]
            max_avg_rpkm = sample_data["max_avg_rpkm"]
            try:
                percent_of_max_rpkm = (float(avg_rpkm) / float(max_avg_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            # Get PSI values
            #if not main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = upstream_exon_text_start_y - 20
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_COVERAGE)

            # Draw exon background
            border_width = 1
            row_canvas.create_rectangle(main_exon_start_x, exon_start_y, main_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor, width=border_width)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            row_canvas.create_rectangle(main_exon_start_x, fill_start_y, main_exon_start_x + exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor, width=border_width)

            # Draw rpkm text. NOTE: Text colors are controlled in ResizingCanvas.on_resize()
            row_canvas.create_text(main_exon_text_start_x + 1, upstream_exon_text_start_y + 1, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(main_exon_text_start_x, upstream_exon_text_start_y, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ########################
            # Draw downstream exon #
            ########################
            #downstream_exon_rpkm = flanking_exons_data[sample_name][next_exon_id]["rpkm"]
            #downstream_exon_max_rpkm = flanking_exons_data[sample_name][next_exon_id]["max_rpkm"]
            downstream_exon_rpkm = sample_data["next_exon_rpkm"]
            downstream_exon_max_rpkm = sample_data["next_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(downstream_exon_rpkm) / float(downstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            downstream_exon_start_x = main_exon_start_x + width_per_exon_container

            # Draw exon background
            row_canvas.create_rectangle(downstream_exon_start_x, exon_start_y, downstream_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            row_canvas.create_rectangle(downstream_exon_start_x, fill_start_y, downstream_exon_start_x + exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)

            # Draw rpkm text
            downstream_exon_text_start_x = downstream_exon_start_x + (exon_width / 2)
            row_canvas.create_text(downstream_exon_text_start_x + 1, upstream_exon_text_start_y + 1, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(downstream_exon_text_start_x, upstream_exon_text_start_y, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            # Update row index for next sample
            row_number += 1

        # Draw the names of the exons above the drawings
        exon_name_canvas.create_text(upstream_exon_text_start_x, top_canvas_height / 2, text=upstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(main_exon_text_start_x, top_canvas_height / 2, text=exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(downstream_exon_text_start_x, top_canvas_height / 2, text=downstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def draw_retained_intron_event(self, data):
        """
        Draws retained introns
        """
        ######################################
        # Fill in samples information column #
        ######################################
        self.populate_samples_frame(data)

        as_id = data["as_id"]
        # Keep track of row number
        row_number = 0

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        intron_container_width = canvas_width / 2
        exon_container_width = canvas_width / 4
        exon_width = exon_container_width * 0.8
        upstream_exon_start_x = (exon_container_width / 2) - (exon_width / 2)
        upstream_exon_stop_x = upstream_exon_start_x + exon_width
        exon_to_intron_distance = upstream_exon_start_x
        intron_stop_x = exon_container_width + intron_container_width + exon_to_intron_distance
        downstream_exon_start_x = intron_stop_x
        downstream_exon_stop_x = downstream_exon_start_x + exon_width
        exon_height = 80
        exon_start_y = (canvas_height - exon_height) / 2
        top_canvas_height = 20

        ################################
        # Draw exon names in top frame #
        ################################
        canvas_background = COLOR_WHITE
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Add upstream exon name
        upstream_exon_name = data["prev_exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."

        # Add exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."
        # Add downstream exon nmae
        downstream_exon_name = data["next_exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + ".."

        ##############
        # Draw exons #
        ##############
        # Get data for this event
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        prev_exon_id = data["prev_exon_id"]
        next_exon_id = data["next_exon_id"]

        # Get rpkm/tot_reads for flanking exons
        #flanking_exons_data = self.data_processor.get_flanking_exons_rpkm_by_exon_ids(sample_names_sorted, prev_exon_id, next_exon_id)
        #main_exon_rpkm_data = self.data_processor.get_main_exon_rpkm_by_asid(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        # Iterate samples
        for sample_name in sorted(data["samples"].keys(), key=natural_sort_key):
            sample_data = data["samples"][sample_name]
            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]

            # Setup colors
            canvas_background = "white"
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE
            intron_color = COLOR_DARKBLUE
            intron_border_color = COLOR_DARKBLUE
            # If event is not reported, draw monochrome canvas
            if not is_reported:
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY
                intron_color = COLOR_DARKGRAY
                intron_border_color = canvas_background
            # Initialize canvas and grid to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            # Keep track of canvases used
            self.canvases.append(row_canvas)
            # Setup tagging buttons
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)
            # Set even weight for every row in the exon frame
            self.exon_frame.rowconfigure(row_number, weight=1)
            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ######################
            # Draw upstream exon #
            ######################
            # Get RPKM values
            #upstream_exon_rpkm = flanking_exons_data[sample_name][prev_exon_id]["rpkm"]
            #upstream_exon_max_rpkm = flanking_exons_data[sample_name][prev_exon_id]["max_rpkm"]
            upstream_exon_rpkm = sample_data["prev_exon_rpkm"]
            upstream_exon_max_rpkm = sample_data["prev_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(upstream_exon_rpkm) / float(upstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            # Draw exon background
            row_canvas.create_rectangle(
                upstream_exon_start_x,
                exon_start_y,
                upstream_exon_stop_x,
                exon_start_y + exon_height,
                fill=canvas_background
            )
            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(
                upstream_exon_start_x,
                fill_start_y,
                upstream_exon_stop_x,
                fill_end_y,
                fill=exon_color,
                outline=exon_bordercolor
            )
            # Draw rpkm text
            upstream_exon_text_start_x = upstream_exon_start_x + (exon_width / 2)
            upstream_exon_text_start_y = exon_height - 10
            row_canvas.create_text(
                upstream_exon_text_start_x + 1,
                upstream_exon_text_start_y + 1,
                text="%.1f" % upstream_exon_rpkm,
                font=self.canvas_font,
                fill=COLOR_CANVAS_TEXT_SHADOW,
                tags=TEXTTAG_SHADOW
            )
            row_canvas.create_text(
                upstream_exon_text_start_x,
                upstream_exon_text_start_y,
                text="%.1f" % upstream_exon_rpkm,
                font=self.canvas_font,
                fill=COLOR_CANVAS_TEXT,
                tags=TEXTTAG_COVERAGE
            )

            ########################
            # Draw downstream exon #
            ########################
            # Get RPKM values
            #downstream_exon_rpkm = flanking_exons_data[sample_name][next_exon_id]["rpkm"]
            #downstream_exon_max_rpkm = flanking_exons_data[sample_name][next_exon_id]["max_rpkm"]
            downstream_exon_rpkm = sample_data["next_exon_rpkm"]
            downstream_exon_max_rpkm = sample_data["next_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(downstream_exon_rpkm) / float(downstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            # Draw exon background
            row_canvas.create_rectangle(
                downstream_exon_start_x,
                exon_start_y,
                downstream_exon_stop_x,
                exon_start_y + exon_height,
                fill=canvas_background,
                outline=exon_bordercolor
            )
            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm / 100) * exon_height))
            row_canvas.create_rectangle(
                downstream_exon_start_x,
                fill_start_y,
                downstream_exon_start_x + exon_width,
                fill_end_y,
                fill=exon_color,
                outline=exon_bordercolor
            )
            # Draw rpkm text
            downstream_exon_text_start_x = downstream_exon_start_x + (exon_width / 2)
            row_canvas.create_text(
                downstream_exon_text_start_x + 1,
                upstream_exon_text_start_y + 1,
                text="%.1f" % (downstream_exon_rpkm),
                font=self.canvas_font,
                fill=COLOR_CANVAS_TEXT_SHADOW,
                tags=TEXTTAG_SHADOW
            )
            row_canvas.create_text(
                downstream_exon_text_start_x,
                upstream_exon_text_start_y,
                text="%.1f" % (downstream_exon_rpkm),
                font=self.canvas_font,
                fill=COLOR_CANVAS_TEXT,
                tags=TEXTTAG_COVERAGE
            )

            ###############
            # Draw intron #
            ###############
            # Get RPKM values for the intron
            #sample_rpkm_data = main_exon_rpkm_data[sample_name]
            #combined_rpkm = sample_rpkm_data["combined_rpkm"]
            #max_combined_rpkm = sample_rpkm_data["max_combined_rpkm"]
            avg_rpkm = sample_data["avg_rpkm"]
            max_avg_rpkm = sample_data["max_avg_rpkm"]
            try:
                percent_of_max_rpkm = (float(avg_rpkm) / float(max_avg_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            intron_start_x = upstream_exon_start_x + exon_width
            # Draw intron coverage fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            # Draw intron fill
            row_canvas.create_rectangle(
                intron_start_x,
                fill_start_y,
                intron_stop_x,
                fill_end_y,
                fill=intron_color,
                outline=intron_border_color
            )
            # Draw coverage text shadow and text
            main_exon_text_start_x = (canvas_width / 2)
            row_canvas.create_text(
                main_exon_text_start_x + 1,
                upstream_exon_text_start_y + 1,
                text="%.1f" % avg_rpkm,
                font=self.canvas_font,
                fill=COLOR_CANVAS_TEXT_SHADOW,
                tags=TEXTTAG_SHADOW
            )
            row_canvas.create_text(
                main_exon_text_start_x,
                upstream_exon_text_start_y,
                text="%.1f" % avg_rpkm,
                font=self.canvas_font,
                fill=COLOR_CANVAS_TEXT,
                tags=TEXTTAG_COVERAGE
            )

            # Draw PSI for main exon (retained intron)
            #if main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = exon_start_y
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ####################################
            # Update row index for next sample #
            ####################################
            row_number += 1

        # Draw the names of the exons above the drawings
        exon_name_canvas.create_text(upstream_exon_text_start_x, top_canvas_height / 2, text=upstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(main_exon_text_start_x, top_canvas_height / 2, text=exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(downstream_exon_text_start_x, top_canvas_height / 2, text=downstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def draw_alternative_donor_events(self, data):
        """
        Draws alternative donor site events.
        """

        # Fill in samples information column
        self.populate_samples_frame(data)

        # Get strand information to know which way to draw exons
        strand = data["strand"]

        # Keep track of current row
        row_number = 0

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        donor_site_width = (canvas_width / 3) * 2
        upstream_exon_width = 60
        main_exon_width = 30
        downstream_exon_width = 60
        exon_height = 80
        exon_start_y = (canvas_height - exon_height) / 2
        top_canvas_height = 20

        #########################################
        # Draw top canvas containing exon names #
        #########################################
        canvas_background = COLOR_WHITE
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Add upstream exon name
        upstream_exon_name = data["prev_exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."

        # Add exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."

        # Add downstream exon name
        downstream_exon_name = data["next_exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + ".."

        ############################################
        ################ Draw exons ################
        ############################################
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        prev_exon_id = data["prev_exon_id"]
        next_exon_id = data["next_exon_id"]
        as_id = data["as_id"]

        # Get rpkm/tot_reads for flanking exons
        #flanking_exons_data = self.data_processor.get_flanking_exons_rpkm_by_exon_ids(sample_names_sorted, prev_exon_id, next_exon_id)
        #main_exon_rpkm_data = self.data_processor.get_main_exon_rpkm_by_asid(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        # Iterate samples
        for sample_name in sorted(data["samples"].keys(), key=natural_sort_key):
            sample_data = data["samples"][sample_name]

            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]

            # Setup colors
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE

            # If event is not reported in sample, draw monochrome canvas
            if not is_reported:
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY

            # Highlight background for sample of interest
            canvas_background = COLOR_WHITE

            # Initialize canvas and grid it to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")

            # Keep track of canvases used
            self.canvases.append(row_canvas)

            # Setup tagging buttons
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)

            # Set even weight for every row in the center frame
            self.exon_frame.rowconfigure(row_number, weight=1)

            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, heigh=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ######################
            # Draw upstream exon #
            ######################
            # Get RPKM values
            #upstream_exon_rpkm = flanking_exons_data[sample_name][prev_exon_id]["rpkm"]
            #upstream_exon_max_rpkm = flanking_exons_data[sample_name][prev_exon_id]["max_rpkm"]
            upstream_exon_rpkm = sample_data["prev_exon_rpkm"]
            upstream_exon_max_rpkm = sample_data["prev_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(upstream_exon_rpkm) / float(upstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            upstream_exon_start_x = (donor_site_width - upstream_exon_width) / 2

            # Draw exon background
            row_canvas.create_rectangle(upstream_exon_start_x, exon_start_y, upstream_exon_start_x + upstream_exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(upstream_exon_start_x, fill_start_y, upstream_exon_start_x + upstream_exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)

            # Draw coverage text
            upstream_exon_text_start_x = upstream_exon_start_x + (upstream_exon_width / 2)
            upstream_exon_text_start_y = exon_height - 10
            row_canvas.create_text(upstream_exon_text_start_x + 1, upstream_exon_text_start_y + 1, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(upstream_exon_text_start_x, upstream_exon_text_start_y, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ##################
            # Draw main exon #
            ##################
            # Get RPKM values
            #sample_rpkm_data = main_exon_rpkm_data[sample_name]
            #combined_rpkm = sample_rpkm_data["combined_rpkm"]
            #max_combined_rpkm = sample_rpkm_data["max_combined_rpkm"]
            avg_rpkm = sample_data["avg_rpkm"]
            max_avg_rpkm = sample_data["max_avg_rpkm"]
            try:
                percent_of_max_rpkm = (float(avg_rpkm) / float(max_avg_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            # Dimensions
            main_exon_start_x = upstream_exon_start_x + upstream_exon_width
            border_width = 1

            # Draw exon background
            row_canvas.create_rectangle(main_exon_start_x, exon_start_y, main_exon_start_x + main_exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor, width=border_width)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            row_canvas.create_rectangle(main_exon_start_x, fill_start_y, main_exon_start_x + main_exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor, width=border_width)

            # Draw rpkm text
            main_exon_text_start_x = main_exon_start_x + (main_exon_width / 2)
            row_canvas.create_text(main_exon_text_start_x + 1, upstream_exon_text_start_y + 1, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(main_exon_text_start_x, upstream_exon_text_start_y, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            # Get PSI values
            #if not main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = upstream_exon_text_start_y - 30
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ########################
            # Draw downstream exon #
            ########################
            # Get RPKM values
            #downstream_exon_rpkm = flanking_exons_data[sample_name][next_exon_id]["rpkm"]
            #downstream_exon_max_rpkm = flanking_exons_data[sample_name][next_exon_id]["max_rpkm"]
            downstream_exon_rpkm = sample_data["next_exon_rpkm"]
            downstream_exon_max_rpkm = sample_data["next_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(downstream_exon_rpkm) / float(downstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            downstream_exon_start_x = donor_site_width  # Not sure if this is right

            # Draw exon background
            row_canvas.create_rectangle(downstream_exon_start_x, exon_start_y, downstream_exon_start_x + downstream_exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            row_canvas.create_rectangle(downstream_exon_start_x, fill_start_y, downstream_exon_start_x + downstream_exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)

            # Draw rpkm text
            downstream_exon_text_start_x = downstream_exon_start_x + (downstream_exon_width / 2)
            row_canvas.create_text(downstream_exon_text_start_x + 1, upstream_exon_text_start_y + 1, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(downstream_exon_text_start_x, upstream_exon_text_start_y, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            # Prepare for next sample
            row_number += 1

        # Draw the names of the exons above the drawings
        exon_name_canvas.create_text(upstream_exon_text_start_x, top_canvas_height / 2, text=upstream_exon_name, font=self.canvas_font, fill=COLOR_EXON_NAME)
        exon_name_canvas.create_text(main_exon_text_start_x, top_canvas_height / 2, text=exon_name, font=self.canvas_font, fill=COLOR_EXON_NAME)
        exon_name_canvas.create_text(downstream_exon_text_start_x, top_canvas_height / 2, text=downstream_exon_name, font=self.canvas_font, fill=COLOR_EXON_NAME)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def draw_alternative_acceptor_events(self, data):
        """
        Draw alternative acceptor events
        """
        # Fill in samples information column
        self.populate_samples_frame(data)

        # Keep track of current row
        row_number = 0

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        donor_site_width = (canvas_width / 3) * 2
        upstream_exon_width = 60
        main_exon_width = 30
        downstream_exon_width = 60
        exon_height = 80
        exon_start_y = (canvas_height - exon_height) / 2
        canvas_background = COLOR_WHITE
        top_canvas_height = 20

        ###########################################
        # Draw a top canvas containing exon names #
        ###########################################
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Get and trim upstream exon name
        upstream_exon_name = data["prev_exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."
        # Get and trim exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."
        # Get and trim downstream exon name
        downstream_exon_name = data["next_exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + ".."

        ########################
        ###### Draw exons ######
        ########################
        # Get data
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        prev_exon_id = data["prev_exon_id"]
        next_exon_id = data["next_exon_id"]
        as_id = data["as_id"]

        # Get expression values for flanking exons and main exon
        #flanking_exons_data = self.data_processor.get_flanking_exons_rpkm_by_exon_ids(sample_names_sorted, prev_exon_id, next_exon_id)
        #main_exon_rpkm_data = self.data_processor.get_main_exon_rpkm_by_asid(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        # Iterate and draw samples
        for sample_name in sample_names_sorted:
            #########################
            ##### Handle sample #####
            #########################
            sample_data = data["samples"][sample_name]
            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]
            # Setup colors
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE
            canvas_background = COLOR_WHITE
            # If event is not reported in sample, draw monochrome canvas
            if not is_reported:
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY
            # Initialize canvas and grid it to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            # Keep track of canvases used
            self.canvases.append(row_canvas)
            # Setup tagging buttons
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)
            # Set event weight for every row in the center frame
            self.exon_frame.rowconfigure(row_number, weight=1)
            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ######################
            # Draw upstream exon #
            ######################
            # Get RPKM values
            #upstream_exon_rpkm = flanking_exons_data[sample_name][prev_exon_id]["rpkm"]
            #upstream_exon_max_rpkm = flanking_exons_data[sample_name][prev_exon_id]["max_rpkm"]
            upstream_exon_rpkm = sample_data["prev_exon_rpkm"]
            upstream_exon_max_rpkm = sample_data["prev_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(upstream_exon_rpkm) / float(upstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0
            # Exon drawing coordinates
            upstream_exon_start_x = (canvas_width / 3) - (upstream_exon_width / 2)
            # Draw exon background
            row_canvas.create_rectangle(upstream_exon_start_x, exon_start_y, upstream_exon_start_x + upstream_exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm / 100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(upstream_exon_start_x, fill_start_y, upstream_exon_start_x + upstream_exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw coverage text
            upstream_exon_text_start_x = upstream_exon_start_x + (upstream_exon_width / 2)
            text_start_y = exon_height - 10
            row_canvas.create_text(upstream_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(upstream_exon_text_start_x, text_start_y, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ##################
            # Draw main exon #
            ##################
            # Get RPKM values
            #sample_rpkm_data = main_exon_rpkm_data[sample_name]
            #combined_rpkm = sample_rpkm_data["combined_rpkm"]
            #max_combined_rpkm = sample_rpkm_data["max_combined_rpkm"]
            avg_rpkm = sample_data["avg_rpkm"]
            max_avg_rpkm = sample_data["max_avg_rpkm"]
            try:
                percent_of_max_rpkm = (float(avg_rpkm) / float(max_avg_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0
            # Dimensions
            main_exon_start_x = canvas_width / 2
            # Draw exon background
            row_canvas.create_rectangle(main_exon_start_x, exon_start_y, main_exon_start_x + main_exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm / 100) * exon_height))
            row_canvas.create_rectangle(main_exon_start_x, fill_start_y, main_exon_start_x + main_exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw RPKM text
            main_exon_text_start_x = main_exon_start_x + (main_exon_width / 2)
            row_canvas.create_text(main_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(main_exon_text_start_x, text_start_y, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)
            # Get and draw PSI text
            #if main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = text_start_y - 30
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ########################
            # Draw downstream exon #
            ########################
            #downstream_exon_rpkm = flanking_exons_data[sample_name][next_exon_id]["rpkm"]
            #downstream_exon_max_rpkm = flanking_exons_data[sample_name][next_exon_id]["max_rpkm"]
            downstream_exon_rpkm = sample_data["next_exon_rpkm"]
            downstream_exon_max_rpkm = sample_data["next_exon_max_rpkm"]
            try:
                percent_of_max_rpkm = (float(downstream_exon_rpkm) / float(downstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0
            downstream_exon_start_x = main_exon_start_x + main_exon_width
            # Draw exon background
            row_canvas.create_rectangle(downstream_exon_start_x, exon_start_y, downstream_exon_start_x + downstream_exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_rpkm/100) * exon_height))
            row_canvas.create_rectangle(downstream_exon_start_x, fill_start_y, downstream_exon_start_x + downstream_exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw rpkm text
            downstream_exon_text_start_x = downstream_exon_start_x + (downstream_exon_width / 2)
            row_canvas.create_text(downstream_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(downstream_exon_text_start_x, text_start_y, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            # Prepare for next sample
            row_number += 1

        # Draw the names of the exons above the drawings
        exon_name_canvas.create_text(upstream_exon_text_start_x, top_canvas_height / 2, text="%s" % upstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(main_exon_text_start_x, top_canvas_height / 2, text=exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(downstream_exon_text_start_x, top_canvas_height / 2, text=downstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def draw_alternative_terminator_event(self, data):
        """
        Draw alternative terminator events.
        """
        # Fill samples information column
        self.populate_samples_frame(data)

        # Keep track of current row
        row_number = 0

        # Dimensions variables
        canvas_width = 300
        canvas_height = 100
        exon_height = 80
        exon_width = 60
        tail_offset_y = 10
        tail_width = 20
        canvas_background = COLOR_WHITE
        top_canvas_height = 20
        """
        The following variables are coordinates for points in a polygon comprising the exon structure. Exon is drawn
        by connecting lines between 8 points, drawn counter-clockwise beginning in the top-left corner:
        
        1 # # # # # # # # # # # # # # 8
        #                             #
        #                             7 # # # # # # # 6
        #                                             #
        #                                             #
        #                                             #
        #                                             #
        #                             4 # # # # # # # 5
        #                             #
        2 # # # # # # # # # # # # # # 3
        """
        one_x, one_y = (canvas_width / 2) - (exon_width + tail_width) / 2, 0
        two_x, two_y = one_x, one_y + exon_height
        three_x, three_y = two_x + exon_width, two_y
        four_x, four_y = three_x, three_y - tail_offset_y
        five_x, five_y = four_x + tail_width, four_y
        six_x, six_y = five_x, one_y + tail_offset_y
        seven_x, seven_y = six_x - tail_width, six_y
        eight_x, eight_y = seven_x, one_y

        ###########################################
        # Draw a top canvas containing exon names #
        ###########################################
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Get and trim exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."

        # Get data
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        as_id = data["as_id"]
        #main_exon_rpkm_data = self.data_processor.get_main_exon_rpkm_by_asid(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        # Iterate and draw samples
        for sample_name in sample_names_sorted:
            #################
            # Handle sample #
            #################
            sample_data = data["samples"][sample_name]
            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]
            # Setup colors
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE
            canvas_background = COLOR_WHITE
            # If event is not reported in sample, draw monochrome canvas
            if not is_reported:
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY
            # Initialize canvas and grid it to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            # Keep track of canvases used
            self.canvases.append(row_canvas)
            # Setup tagging buttons
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)
            # Set event weight for every row in the center frame
            self.exon_frame.rowconfigure(row_number, weight=1)
            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ##################
            # Draw main exon #
            ##################
            row_canvas.create_polygon(
                [
                    one_x, one_y,
                    two_x, two_y,
                    three_x, three_y,
                    four_x, four_y,
                    five_x, five_y,
                    six_x, six_y,
                    seven_x, seven_y,
                    eight_x, eight_y,
                    one_x, one_y
                ],
                fill=canvas_background,
                outline=exon_bordercolor
            )

            # Get RPKM values
            #sample_rpkm_data = main_exon_rpkm_data[sample_name]
            #combined_rpkm = sample_rpkm_data["combined_rpkm"]
            #max_combined_rpkm = sample_rpkm_data["max_combined_rpkm"]
            avg_rpkm = sample_data["avg_rpkm"]
            max_avg_rpkm = sample_data["max_avg_rpkm"]
            try:
                percent_of_max_rpkm = (float(avg_rpkm) / float(max_avg_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            # Calc fill Y position
            fill_start_y = (one_y + two_y) - (int((percent_of_max_rpkm / 100) * exon_height))

            # Draw fill
            fill_polygons = [
                one_x, fill_start_y,
                two_x, two_y,
                three_x, three_y
            ]
            if (exon_height - fill_start_y) < (exon_height - four_y):
                fill_polygons += [four_x, fill_start_y]
            elif (exon_height - fill_start_y) < (exon_height - six_y):
                fill_polygons += [four_x, four_y, five_x, five_y, six_x, fill_start_y]
            else:
                fill_polygons += [four_x, four_y, five_x, five_y, six_x, six_y, seven_x, seven_y, eight_x, fill_start_y]
            # Lastly, add a line back to the first point
            fill_polygons += [one_x, fill_start_y]
            row_canvas.create_polygon(fill_polygons, fill=exon_color, outline=exon_bordercolor)

            # Draw RPKM text
            main_exon_text_start_x = one_x + (exon_width / 2)
            main_exon_text_start_y = two_y - 10
            row_canvas.create_text(main_exon_text_start_x + 1, main_exon_text_start_y + 1, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(main_exon_text_start_x, main_exon_text_start_y, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            # Get PSI values
            #if main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = main_exon_text_start_y - 30
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(main_exon_text_start_x + 1, psi_text_start_y + 1, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
                row_canvas.create_text(main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

        # Draw the name of the exon in the top canvas
        if len(sample_names_sorted) > 0:
            # TODO: Do this check in every draw-function
            exon_name_canvas.create_text(main_exon_text_start_x, top_canvas_height / 2, text=exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def draw_alternative_promotor_event(self, data):
        """
        Draw alternative promotor event.
        """
        # Fill samples information column
        self.populate_samples_frame(data)

        # Keep track of current row
        row_number = 0

        # Dimensions variables
        canvas_width = 300
        canvas_height = 100
        exon_height = 80
        exon_width = 60
        tail_offset_y = 10
        tail_width = 20
        canvas_background = COLOR_WHITE
        top_canvas_height = 20
        """
        The following variables are coordinates for points (x, y) in a polygon comprising the exon structure. Exon is
        drawn by connecting lines between 8 points, drawn clockwise beginning in the top-right corner:
        
                             8 # # # # # # # # # # # # # # # # # 1
                             #                                   #
             6 # # # # # # # 7                                   #
             #                                                   #
             #                                                   #
             #                                                   #
             #                                                   #
             5 # # # # # # # 4                                   #
                             #                                   #
                             3 # # # # # # # # # # # # # # # # # 2                                  
        """
        one_x, one_y = (canvas_width / 2) + ((exon_width + tail_width) / 2), 0
        two_x, two_y = one_x, one_y + exon_height
        three_x, three_y = two_x - exon_width, two_y
        four_x, four_y = three_x, three_y - tail_offset_y
        five_x, five_y = four_x - tail_width, four_y
        six_x, six_y = five_x, one_y + tail_offset_y
        seven_x, seven_y = six_x + tail_width, six_y
        eight_x, eight_y = seven_x, one_y

        ###########################################
        # Draw a top canvas containing exon names #
        ###########################################
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Get and trim exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."

        # Get data
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        as_id = data["as_id"]
        #main_exon_rpkm_data = self.data_processor.get_main_exon_rpkm_by_asid(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        ##################################
        # Iterate samples and draw exons #
        ##################################
        for sample_name in sample_names_sorted:
            sample_data = data["samples"][sample_name]
            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]
            # Setup colors
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE
            canvas_background = COLOR_WHITE
            # If event is not reported in sample, draw monochrome canvas
            if not is_reported:
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY
            # Initialize canvas and grid it to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            # Keep track of canvases used
            self.canvases.append(row_canvas)
            # Setup tagging buttons
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)
            # Set event weight for every row in the center frame
            self.exon_frame.rowconfigure(row_number, weight=1)
            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ##################
            # Draw main exon #
            ##################
            row_canvas.create_polygon(
                [
                    one_x, one_y,
                    two_x, two_y,
                    three_x, three_y,
                    four_x, four_y,
                    five_x, five_y,
                    six_x, six_y,
                    seven_x, seven_y,
                    eight_x, eight_y,
                    one_x, one_y
                ],
                fill=canvas_background,
                outline=exon_bordercolor
            )

            # Get RPKM values
            #sample_rpkm_data = main_exon_rpkm_data[sample_name]
            #combined_rpkm = sample_rpkm_data["combined_rpkm"]
            #max_combined_rpkm = sample_rpkm_data["max_combined_rpkm"]
            avg_rpkm = sample_data["avg_rpkm"]
            max_avg_rpkm = sample_data["max_avg_rpkm"]
            try:
                percent_of_max_rpkm = (float(avg_rpkm) / float(max_avg_rpkm)) * 100
            except ZeroDivisionError:
                percent_of_max_rpkm = 0

            # Calc fill Y position
            fill_start_y = (one_y + two_y) - (int((percent_of_max_rpkm / 100) * exon_height))

            # Draw fill
            fill_polygons = [
                one_x, fill_start_y,
                two_x, two_y,
                three_x, three_y
            ]
            if (exon_height - fill_start_y) < (exon_height - four_y):
                fill_polygons += [four_x, fill_start_y]
            elif (exon_height - fill_start_y) < (exon_height - six_y):
                fill_polygons += [four_x, four_y, five_x, five_y, six_x, fill_start_y]
            else:
                fill_polygons += [four_x, four_y, five_x, five_y, six_x, six_y, seven_x, seven_y, eight_x, fill_start_y]
            # Lastly, add a line back to the first point
            fill_polygons += [one_x, fill_start_y]
            row_canvas.create_polygon(fill_polygons, fill=exon_color, outline=exon_bordercolor)

            # Draw RPKM text
            main_exon_text_start_x = eight_x + (exon_width / 2)
            main_exon_text_start_y = two_y - 10
            row_canvas.create_text(main_exon_text_start_x + 1, main_exon_text_start_y + 1, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(main_exon_text_start_x, main_exon_text_start_y, text="%.1f" % avg_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            # Get PSI values
            #if main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = main_exon_text_start_y - 30
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(main_exon_text_start_x + 1, psi_text_start_y + 1, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
                row_canvas.create_text(main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

        # Draw the name of the exon in the top canvas
        if len(sample_names_sorted) > 0:
            # TODO: Do this check in every draw-function
            exon_name_canvas.create_text(main_exon_text_start_x, top_canvas_height / 2, text=exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def draw_mutually_exclusive_exons_event(self, data):
        """
        Draws mutually exclusive exons events.
        """
        # Fill samples information column
        self.populate_samples_frame(data)

        # Keep track of current row
        row_number = 0

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        exon_height = 60
        exon_width = 60
        canvas_background = COLOR_WHITE
        top_canvas_height = 20
        exon_start_y = (canvas_height - exon_height) / 2
        canvas_per_exon = canvas_width / 4
        exon_padding = (canvas_per_exon - exon_width) / 2
        text_start_y = exon_height - 1

        #########################################
        # Draw top canvas containing exon names #
        #########################################
        exon_name_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=top_canvas_height)
        exon_name_canvas.grid(column=0, row=row_number, sticky="NEWS")
        row_number += 1

        # Get and trim exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."
        # Get and trim upstream exon name
        upstream_exon_name = data["prev_exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."
        # Get and trim downstream exon name
        downstream_exon_name = data["next_exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + [".."]

        # Get data
        sample_names_sorted = sorted(data["samples"].keys(), key=natural_sort_key)
        as_id = data["as_id"]
        prev_exon_id = data["prev_exon_id"]
        next_exon_id = data["next_exon_id"]
        # Get expression values for flanking exons and main exons
        #flanking_exons_data = self.data_processor.get_flanking_exons_rpkm_by_exon_ids(sample_names_sorted, prev_exon_id, next_exon_id)
        main_exon_rpkm_data = self.data_processor.get_rpkm_for_mutually_exclusive_exons(sample_names_sorted, as_id)
        #main_exon_psi_data = self.data_processor.get_main_exon_psi_by_asid(sample_names_sorted, as_id)

        ##################################
        # Iterate samples and draw exons #
        ##################################
        for sample_name in sample_names_sorted:
            sample_data = data["samples"][sample_name]
            # Sample-specific variables
            is_reported = sample_data["is_reported"]
            sample_tag = sample_data["event_tag"]
            # Setup colors
            exon_color = COLOR_BLUE
            exon_bordercolor = COLOR_DARKBLUE
            canvas_background = COLOR_WHITE
            # If event is not reported in sample, draw monochrome canvas
            if not is_reported:
                exon_color = COLOR_DARKGRAY
                exon_bordercolor = COLOR_DARKGRAY
            # Initialize canvas and grid it to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            # Keep track of canvases used
            self.canvases.append(row_canvas)
            # Setup tagging buttons
            self.add_tagging_buttons(row_number, sample_name, is_reported, sample_tag, as_id)
            # Set event weight for every row in the center frame
            self.exon_frame.rowconfigure(row_number, weight=1)
            # Add separator
            row_number += 1
            sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            sep.grid(row=row_number, column=0, sticky="NEWS")
            button_sep = tk.Frame(self.exon_frame, bg=COLOR_DARKWHITE, height=2)
            button_sep.grid(row=row_number, column=1, sticky="NEWS")

            ######################
            # Draw upstream exon #
            ######################
            upstream_exon_start_x = exon_padding
            #upstream_exon_rpkm = flanking_exons_data[sample_name][prev_exon_id]["rpkm"]
            #upstream_exon_max_rpkm = flanking_exons_data[sample_name][prev_exon_id]["max_rpkm"]
            upstream_exon_rpkm = sample_data["prev_exon_rpkm"]
            upstream_exon_max_rpkm = sample_data["prev_exon_max_rpkm"]
            try:
                upstream_percent_of_max_rpkm = (float(upstream_exon_rpkm) / float(upstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                upstream_percent_of_max_rpkm = 0
            # Draw exon background
            row_canvas.create_rectangle(upstream_exon_start_x, exon_start_y, upstream_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            fill_start_y = (exon_start_y + exon_height) - (int((upstream_percent_of_max_rpkm / 100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            # Draw exon fill
            row_canvas.create_rectangle(upstream_exon_start_x, fill_start_y, upstream_exon_start_x + exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw exon RPKM values
            upstream_exon_text_start_x = upstream_exon_start_x + (exon_width / 2)
            row_canvas.create_text(upstream_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(upstream_exon_text_start_x, text_start_y, text="%.1f" % upstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ########################
            # Draw first main exon #
            ########################
            # Get combined RPKM for first exon(s)
            first_main_exon_rpkm = main_exon_rpkm_data[sample_name]["first_exon"]["combined_rpkm"]
            first_main_exon_max_rpkm = main_exon_rpkm_data[sample_name]["first_exon"]["max_combined_rpkm"]
            try:
                first_main_exon_percent_of_max_rpkm = (float(first_main_exon_rpkm) / float(first_main_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                first_main_exon_percent_of_max_rpkm = 0
            first_main_exon_name = main_exon_rpkm_data[sample_name]["first_exon"]["joined_exon_name"]
            # Draw exon background
            first_main_exon_start_x = upstream_exon_start_x + exon_width + (exon_padding * 2)
            first_main_exon_end_x = first_main_exon_start_x + exon_width
            row_canvas.create_rectangle(first_main_exon_start_x, exon_start_y, first_main_exon_end_x, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            # Draw exon fill
            first_main_exon_fill_start_y = (exon_start_y + exon_height) - (int((first_main_exon_percent_of_max_rpkm / 100) * exon_height))
            first_main_exon_fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(first_main_exon_start_x, first_main_exon_fill_start_y, first_main_exon_end_x, first_main_exon_fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw exon RPKM values
            first_main_exon_text_start_x = first_main_exon_start_x + (exon_width / 2)
            row_canvas.create_text(first_main_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % first_main_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(first_main_exon_text_start_x, text_start_y, text="%.1f" % first_main_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)
            # Draw PSI values
            #if main_exon_psi_data[sample_name]["is_reported"]:
            if sample_data["is_reported"]:
                #sample_psi_data = main_exon_psi_data[sample_name]
                #sample_psi = sample_psi_data["psi"]
                sample_psi = sample_data["psi"]
                #sample_included_counts = sample_psi_data["included_counts"]
                sample_included_counts = sample_data["included_counts"]
                #sample_excluded_counts = sample_psi_data["excluded_counts"]
                sample_excluded_counts = sample_data["excluded_counts"]
                psi_text_start_y = exon_start_y - 10
                psi_text = "PSI: %.2f (%d/%d)" % (sample_psi, sample_included_counts, sample_excluded_counts)
                row_canvas.create_text(first_main_exon_text_start_x, psi_text_start_y, text=psi_text, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW)

            #########################
            # Draw second main exon #
            #########################
            # Get combined RPKM for second exon(s)
            second_main_exon_rpkm = main_exon_rpkm_data[sample_name]["second_exon"]["combined_rpkm"]
            second_main_exon_max_rpkm = main_exon_rpkm_data[sample_name]["second_exon"]["max_combined_rpkm"]
            try:
                second_main_exon_percent_of_max_rpkm = (float(second_main_exon_rpkm) / float(second_main_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                second_main_exon_percent_of_max_rpkm = 0
            second_main_exon_name = main_exon_rpkm_data[sample_name]["second_exon"]["joined_exon_name"]
            # Draw exon background
            second_main_exon_start_x = first_main_exon_start_x + exon_width + (exon_padding * 2)
            second_main_exon_end_x = second_main_exon_start_x + exon_width
            row_canvas.create_rectangle(second_main_exon_start_x, exon_start_y, second_main_exon_end_x, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            # Draw exon fill
            second_main_exon_fill_start_y = (exon_start_y + exon_height) - (int((second_main_exon_percent_of_max_rpkm / 100) * exon_height))
            second_main_exon_fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(second_main_exon_start_x, second_main_exon_fill_start_y, second_main_exon_end_x, second_main_exon_fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw exon RPKM values
            second_main_exon_text_start_x = second_main_exon_start_x + (exon_width / 2)
            row_canvas.create_text(second_main_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % second_main_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(second_main_exon_text_start_x, text_start_y, text="%.1f" % second_main_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

            ########################
            # Draw downstream exon #
            ########################
            downstream_exon_start_x = second_main_exon_end_x + (exon_padding * 2)
            #downstream_exon_rpkm = flanking_exons_data[sample_name][next_exon_id]["rpkm"]
            #downstream_exon_max_rpkm = flanking_exons_data[sample_name][next_exon_id]["max_rpkm"]
            downstream_exon_rpkm = sample_data["next_exon_rpkm"]
            downstream_exon_max_rpkm = sample_data["next_exon_max_rpkm"]
            try:
                downstream_percent_of_max_rpkm = (float(downstream_exon_rpkm) / float(downstream_exon_max_rpkm)) * 100
            except ZeroDivisionError:
                downstream_percent_of_max_rpkm = 0
            # Draw exon background
            row_canvas.create_rectangle(downstream_exon_start_x, exon_start_y, downstream_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, outline=exon_bordercolor)
            fill_start_y = (exon_start_y + exon_height) - (int((downstream_percent_of_max_rpkm / 100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            # Draw exon fill
            row_canvas.create_rectangle(downstream_exon_start_x, fill_start_y, downstream_exon_start_x + exon_width, fill_end_y, fill=exon_color, outline=exon_bordercolor)
            # Draw exon RPKM values
            downstream_exon_text_start_x = downstream_exon_start_x + (exon_width / 2)
            row_canvas.create_text(downstream_exon_text_start_x + 1, text_start_y + 1, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT_SHADOW, tags=TEXTTAG_SHADOW)
            row_canvas.create_text(downstream_exon_text_start_x, text_start_y, text="%.1f" % downstream_exon_rpkm, font=self.canvas_font, fill=COLOR_CANVAS_TEXT, tags=TEXTTAG_COVERAGE)

        # Draw the names of the exons above the drawings
        exon_name_canvas.create_text(upstream_exon_text_start_x, top_canvas_height / 2, text=upstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(first_main_exon_text_start_x, top_canvas_height / 2, text=first_main_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(second_main_exon_text_start_x, top_canvas_height / 2, text=second_main_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)
        exon_name_canvas.create_text(downstream_exon_text_start_x, top_canvas_height / 2, text=downstream_exon_name, fill=COLOR_EXON_NAME, font=self.canvas_font)

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

if __name__ == "__main__":

    # Whether or not to draw a diagram of the callstack. Slows down execution, but useful for debugging
    enable_callgraph = False

    if enable_callgraph:
        from pycallgraph import PyCallGraph
        from pycallgraph.output import GraphvizOutput
        from pycallgraph import Config
        from pycallgraph import GlobbingFilter

        # Exclude Pandas-functions from output
        config = Config()
        config.trace_filter = GlobbingFilter(exclude=[
            "pycallgraph.*",
            "pandas.*",
            "numpy.*"
        ])

        with PyCallGraph(output=GraphvizOutput(), config=config):
            # Create and run app
            app = TINTagger()
            app.wm_title("TIN-Tagger")
            screen_width = app.winfo_screenwidth()
            screen_height = app.winfo_screenheight()
            # Fill 80% of screen
            new_width = int((float(screen_width) / 100.00) * 80)
            new_height = int((float(screen_height) / 100.00) * 80)
            app.geometry("%dx%d" % (new_width, new_height))
            #app.withdraw()
            app.center_window(app)
            #app.update()
            #app.deiconify()
            #app.geometry("1024x576")

            app.mainloop()
    else:
        # Create and run app
        app = TINTagger()
        app.wm_title("TIN-Tagger")
        screen_width = app.winfo_screenwidth()
        screen_height = app.winfo_screenheight()
        # Fill 80% of screen
        new_width = int((float(screen_width) / 100.00) * 80)
        new_height = int((float(screen_height) / 100.00) * 80)
        app.geometry("%dx%d" % (new_width, new_height))
        app.center_window(app)
        app.mainloop()



