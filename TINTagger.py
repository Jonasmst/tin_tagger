import Tkinter as tk
import ttk
import sys
import re
from tkFileDialog import askopenfilename, asksaveasfilename
import pandas as pd
import os
import copy
import random
import subprocess
import tkMessageBox
from TINDataProcessor import TINDataProcessor

# TODO: Cache rows so we don't need to run samtools etc when pressing previous-button
# TODO: Find PSI, included counts, excluded counts for other exons and display if it's available (gonna be an sql-call).

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
COLOR_DARKGRAY = "#7f8c8d"
COLOR_ORANGE = "#e67e22"
COLOR_CANVAS_TEXT = COLOR_WHITE
COLOR_CANVAS_TEXT_SHADOW = "black"
COLOR_SAMPLE_HIGHLIGHT = COLOR_LIGHTRED

COLOR_INTERESTING = COLOR_GREEN
COLOR_NOT_INTERESTING = COLOR_RED
COLOR_UNCERTAIN = COLOR_PURPLE
COLOR_NO_TAG = COLOR_DARKGRAY

TAG_INTERESTING = 0
TAG_NOT_INTERESTING = 1
TAG_UNCERTAIN = 2
TAG_NO_TAG = -1

STYLE_BUTTON_INTERESTING = "Interesting.TButton"
STYLE_BUTTON_NOT_INTERESTING = "Not_interesting.TButton"
STYLE_BUTTON_UNCERTAIN = "Uncertain.TButton"


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
        self.delete("text_shadow")
        for item in self.find_all():
            if self.type(item) == "text":
                # Get current font
                font = self.itemcget(item, "font")

                # Get the current text
                current_text = self.itemcget(item, "text")

                # Get coordinates of the text object
                text_x, text_y = self.coords(item)

                # Draw shadow first, then a new instance of the text
                self.create_text(text_x + 1, text_y + 1, text=current_text, tags="text_shadow", font=font, fill=COLOR_CANVAS_TEXT_SHADOW)
                self.create_text(text_x, text_y, text=current_text, font=font, fill=COLOR_CANVAS_TEXT)

                # Finally, delete this item
                self.delete(item)


class TINTagger(tk.Tk):
    def __init__(self, *args, **kwargs):
        # Super init
        tk.Tk.__init__(self, *args, **kwargs)

        # TEST: Different styles (frame backgrounds)
        #s = ttk.Style()
        #s.theme_use("classic")
        #s.configure("TFrame", background=COLOR_DARKBLUE)
        #s.configure("TLabel", background=COLOR_DARKBLUE, foreground=COLOR_WHITE)
        # END TEST

        # Set testing or not
        self.testing = True

        if self.testing:
            # Mac ghetto-fix for bringing GUI to front when running.
            os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
            pass

        #####################
        # Various variables #
        #####################
        # Keep track of all of the canvases in the center window
        self.canvases = []

        # A processor to handle I/O and system calls
        self.data_processor = TINDataProcessor()

        # No dataset loaded by default
        self.dataset = None

        # Start on row 0 by default
        self.current_row_index = 0

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
        self.bind("<Left>", self.left_arrow_clicked)
        self.bind("<Right>", self.right_arrow_clicked)

        #################################
        # Draw graphical user interface #
        #################################
        # Create a menu
        self.create_menu()

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
        self.load_frame = self.create_load_frame()

        # Create statusbar at the bottom
        self.statusbar = self.create_statusbar()

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
        exon_frame = ttk.Frame(self.left_frame, padding=0, borderwidth=1, relief=tk.SOLID)
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

        # Sample
        sample_label = ttk.Label(sidebar_information, text="Sample name:", font=label_font)
        sample_label.grid(column=0, row=current_row, sticky="W")
        self.sample_text = ttk.Label(sidebar_information, text="<Sample name>", font=text_font)
        self.sample_text.grid(column=1, row=current_row, sticky="W", columnspan=2)
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

        # PSI
        psi_label = ttk.Label(sidebar_information, text="PSI:", font=label_font)
        psi_label.grid(column=0, row=current_row, sticky=tk.W)
        self.psi_text = ttk.Label(sidebar_information, text="<PSI>", font=text_font)
        self.psi_text.grid(column=1, row=current_row, sticky=tk.W, columnspan=2)
        current_row += 1

        # Included counts
        included_counts_label = ttk.Label(sidebar_information, text="Included counts:", font=label_font)
        included_counts_label.grid(column=0, row=current_row, sticky="W")
        self.included_counts_text = ttk.Label(sidebar_information, text="<Included counts>", font=text_font)
        self.included_counts_text.grid(column=1, row=current_row, columnspan=2, sticky="W")
        current_row += 1

        # Excluded counts
        excluded_counts_label = ttk.Label(sidebar_information, text="Excluded counts:", font=label_font)
        excluded_counts_label.grid(column=0, row=current_row, sticky="W")
        self.excluded_counts_text = ttk.Label(sidebar_information, text="<Excluded counts>", font=text_font)
        self.excluded_counts_text.grid(column=1, row=current_row, columnspan=2, sticky="W")
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
        # Display prompt to load file
        load_frame = ttk.Frame(self.main_frame)
        load_frame.grid(column=0, row=0, sticky="NEWS")
        load_label = ttk.Label(load_frame, text="No dataset found, please load one: ")
        load_label.grid(column=0, row=0, sticky="E")
        load_button = ttk.Button(load_frame, text="Open file..", command=self.open_file)
        load_button.grid(column=1, row=0, sticky="W")

        load_frame.columnconfigure(0, weight=1)
        load_frame.columnconfigure(1, weight=1)
        load_frame.rowconfigure(0, weight=1)

        return load_frame

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
        sub_menu.add_command(label="Quit", accelerator="Escape", command=lambda: self.close_window(event=None))

        # Add dataset menu
        # TODO: Add an "information" dataset entry that shows the number of events, number of samples, and column names
        dataset_menu = tk.Menu(main_menu)
        main_menu.add_cascade(label="Dataset", menu=dataset_menu)
        dataset_menu.add_command(label="Filter dataset", command=self.show_filter_dataset_window, accelerator="Ctrl+F")
        dataset_menu.add_separator()
        dataset_menu.add_command(label="Open filters..", command=self.read_dataset_filters)
        dataset_menu.add_command(label="Save current filters", command=self.save_dataset_filters)

    def show_filter_dataset_window(self):

        # Alter a copy of the current settings and only apply them if instructed by the user
        native_filters = self.convert_filters_to_native_types()  # copy.deepcopy() doesn't understand tk.IntVars
        native_filters_deepcopy = copy.deepcopy(native_filters)
        filters_copy = self.convert_filters_to_tk_specific_datatypes(native_filters_deepcopy)

        ############################################
        # Create a window to display everything in #
        ############################################
        filter_window = tk.Toplevel()
        filter_window.bind("<Escape>", lambda event=None: filter_window.destroy())
        filter_window.wm_title("Filter dataset")
        filter_window.geometry("400x500")

        # A frame to hold the filtering options
        filter_frame = ttk.Frame(filter_window, padding=(0, 0, 0, 100))
        filter_frame.grid(column=0, row=0, sticky="NEWS")
        filter_frame.focus_force()  # Force focus to a widget in this window so that binds work

        # A frame on the bottom to hold the buttons
        button_frame = ttk.Frame(filter_window)
        button_frame.grid(column=0, row=1, sticky="NEWS")

        # Set weights for the rows so that filter frame expands, but not button frame.
        filter_window.rowconfigure(0, weight=1)
        filter_window.columnconfigure(0, weight=1)
        filter_frame.columnconfigure(0, weight=1)

        # Add a statusbar to the filters window
        filter_statusbar = ttk.Frame(filter_window, borderwidth=1, relief=tk.SUNKEN)
        filter_statusbar.grid(column=0, row=2, sticky="NEWS")
        self.filter_statusbar_label = ttk.Label(filter_statusbar, text="All filters are valid", font="tkDefaultFont")
        self.filter_statusbar_label.grid(column=0, row=0, sticky="NEWS")

        # Validation functions (see http://infohost.nmt.edu/tcc/help/pubs/tkinter/web/entry-validation.html)
        validate_positive_integer = (self.register(self.validate_positive_integer), "%P", "%W", "%d", "%V")
        validate_positive_float = (self.register(self.validate_positive_float), "%P", "%W", "%d", "%V")

        # Keep track of the rows we've assigned stuff to in the filter frame
        current_filter_row = 0

        #################
        # Setup filters #
        #################

        # Filter entry for included counts
        inc_counts_frame = ttk.Frame(filter_frame)
        inc_counts_frame.grid(column=0, row=current_filter_row, sticky="NEWS")
        inc_counts_label = ttk.Label(inc_counts_frame, text="Min. included counts:")
        inc_counts_label.grid(column=0, row=0, sticky="NEWS")
        inc_counts_entry = ttk.Entry(
            inc_counts_frame,
            textvariable=filters_copy["included_counts"],
            exportselection=0,
            justify=tk.LEFT,
            validatecommand=validate_positive_integer,
            validate="all",
            name="included_counts_filter"
        )
        inc_counts_entry.grid(column=1, row=0, sticky="NEWS")
        inc_counts_frame.columnconfigure(1, weight=1)

        current_filter_row += 1

        # Filter entry for excluded counts
        excl_counts_frame = ttk.Frame(filter_frame)
        excl_counts_frame.grid(column=0, row=current_filter_row, sticky="NEWS")
        excl_counts_label = ttk.Label(excl_counts_frame, text="Min. excluded counts:")
        excl_counts_label.grid(column=0, row=0, sticky="NEWS")
        excl_counts_entry = ttk.Entry(
            excl_counts_frame,
            textvariable=filters_copy["excluded_counts"],
            exportselection=0,
            justify=tk.LEFT,
            validatecommand=validate_positive_integer,
            validate="all",
            name="excluded_counts_filter"
        )
        excl_counts_entry.grid(column=1, row=0, sticky="NEWS")
        excl_counts_frame.columnconfigure(1, weight=1)

        current_filter_row += 1

        # Filter entry for PSI
        psi_frame = ttk.Frame(filter_frame)
        psi_frame.grid(column=0, row=current_filter_row, sticky="NEWS")
        psi_label = ttk.Label(psi_frame, text="Min. PSI:")
        psi_label.grid(column=0, row=0)
        psi_entry = ttk.Entry(
            psi_frame,
            textvariable=filters_copy["psi"],
            exportselection=0,
            justify=tk.LEFT,
            validatecommand=validate_positive_float,
            validate="all",
            name="psi_filter"
        )
        psi_entry.grid(column=1, row=0, sticky="NEWS")
        psi_frame.columnconfigure(1, weight=1)

        current_filter_row += 1

        # Filter entry for gene RPKM
        rpkm_frame = ttk.Frame(filter_frame)
        rpkm_frame.grid(column=0, row=current_filter_row, sticky="NEWS")
        rpkm_label = ttk.Label(rpkm_frame, text="Min. gene RPKM:")
        rpkm_label.grid(column=0, row=0)
        rpkm_entry = ttk.Entry(
            rpkm_frame,
            textvariable=filters_copy["rpkm"],
            exportselection=0,
            justify=tk.LEFT,
            validatecommand=validate_positive_float,
            validate="all",
            name="rpkm_filter"
        )
        rpkm_entry.grid(column=1, row=0, sticky="NEWS")
        rpkm_frame.columnconfigure(1, weight=1)

        current_filter_row += 1

        # Filter choices for splice types
        splice_type_frame = ttk.Frame(filter_frame)
        splice_type_frame.grid(column=0, row=current_filter_row, sticky="NEWS")
        splice_type_label = ttk.Label(splice_type_frame, text="Splice types:")
        # Keep track of splice type row
        current_splice_type_row = 0
        splice_type_label.grid(column=0, row=current_splice_type_row, sticky="NEWS")
        current_splice_type_row += 1

        # Setup checkboxes
        for splice_type in sorted(filters_copy["splice_type"].keys()):
            # Set on or off
            filters_copy["splice_type"][splice_type]["enabled_var"].set(filters_copy["splice_type"][splice_type]["enabled"])
            check_button = ttk.Checkbutton(
                splice_type_frame,
                text=filters_copy["splice_type"][splice_type]["description"],
                variable=filters_copy["splice_type"][splice_type]["enabled_var"],
                onvalue=1,
                offvalue=0
            )
            check_button.grid(column=0, row=current_splice_type_row, sticky="W")
            current_splice_type_row += 1

        current_filter_row += 1

        ################################
        # Add Apply and Cancel buttons #
        ################################
        cancel_button = ttk.Button(button_frame, text="Cancel", command=filter_window.destroy)
        cancel_button.grid(column=0, row=0, sticky="E")
        #self.apply_button = ttk.Button(button_frame, text="Apply", command=self.apply_filters)
        self.apply_button = ttk.Button(button_frame, text="Apply", command=lambda: self.apply_filters(filters_copy))
        self.apply_button.grid(column=1, row=0, sticky="E")
        button_frame.columnconfigure(0, weight=1)
        #button_frame.columnconfigure(1, weight=1)

    def validate_positive_integer(self, string_value, widget_name, edit_type, edit_reason):
        """
        Validates entry fields. Returns true if the text represents a positive integer, false if not.
        """
        style = ttk.Style()
        style.configure("Red.TEntry", foreground=COLOR_RED)
        style.configure("Green.TEntry", foreground=COLOR_GREEN)
        widget = self.nametowidget(widget_name)
        short_name = widget_name.split(".")[-1]
        self.valid_filters[short_name] = 1

        # Check if the new value will be a positive integer
        if string_value.isdigit():
            # Valid (note that isdigit() returns false for negative values
            if edit_reason == "focusout":
                widget.configure(style="TEntry")
            else:
                widget.configure(style="Green.TEntry")
        else:
            # Invalid
            widget.configure(style="Red.TEntry")
            self.filter_statusbar_label["text"] = "Invalid filter: Must be positive integer."
            self.valid_filters[short_name] = 0

        # If all filters are valid, enable the apply-button
        if all(valid == 1 for valid in self.valid_filters.values()):
            self.filter_statusbar_label["text"] = "All filters are valid"
            self.apply_button["state"] = "enabled"
        else:
            self.apply_button["state"] = "disabled"
            if edit_reason == "focusin":
                self.filter_statusbar_label["text"] = "There are invalid filters."

        return True

    def validate_positive_float(self, string_value, widget_name, edit_type, edit_reason):
        """
        Validates entry fields. Returns true if the text represents a positive float, false if not.
        """
        style = ttk.Style()
        style.configure("Red.TEntry", foreground=COLOR_RED)
        style.configure("Green.TEntry", foreground=COLOR_GREEN)
        widget = self.nametowidget(widget_name)
        short_name = widget_name.split(".")[-1]
        self.valid_filters[short_name] = 1

        try:
            float_val = float(string_value)
            if float_val >= 0.0:
                # Valid
                if edit_reason == "focusout":
                    widget.configure(style="TEntry")
                else:
                    widget.configure(style="Green.TEntry")
            else:
                # Invalid
                self.valid_filters[short_name] = 0
                widget.configure(style="Red.TEntry")
                self.filter_statusbar_label["text"] = "Invalid filter: Negative values not allowed."

        except ValueError:
            # Invalid
            self.valid_filters[short_name] = 0
            self.filter_statusbar_label["text"] = "Invalid filter: Non-negative floats only (e.g. 4.65)"
            widget.configure(style="Red.TEntry")

        # Enable or disable apply-button
        if all(valid == 1 for valid in self.valid_filters.values()):
            self.filter_statusbar_label["text"] = "All filters are valid"
            self.apply_button["state"] = "enabled"
        else:
            self.apply_button["state"] = "disabled"
            if edit_reason == "focusin":
                self.filter_statusbar_label["text"] = "There are invalid filters."
        return True

    def get_default_filters(self):
        """
        Returns a default filter: Include all splice types, all minimum values are 0.
        """
        dataset_filter = {
            "included_counts": tk.StringVar(),
            "excluded_counts": tk.StringVar(),
            "psi": tk.StringVar(),
            "rpkm": tk.StringVar(),
            "splice_type": {
                "AP": {
                    "description": "Alternative Promoter (AP)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                },
                "AD": {
                    "description": "Alternative Donor Site (AD)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                },
                "AA": {
                    "description": "Alternative Acceptor Site (AA)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                },
                "ES": {
                    "description": "Exon Skipping (ES)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                },
                "AT": {
                    "description": "Alternative Terminator (AT)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                },
                "RI": {
                    "description": "Retained Intron (RI)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                },
                "ME": {
                    "description": "Mutually Exclusive Exons (ME)",
                    "enabled": 1,
                    "enabled_var": tk.IntVar()
                }
            }
        }

        # Set default values
        dataset_filter["included_counts"].set("0")
        dataset_filter["excluded_counts"].set("0")
        dataset_filter["psi"].set("0.00")
        dataset_filter["rpkm"].set("0.00")

        return dataset_filter

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

    def apply_filters(self, filters=None):
        """
        Applies the currently active filter. Calls the TINDataProcessor to return a filtered dataset and updates
        self.dataset to reflect the changes. Also resets the self.current_row_index to 0. Finally, calls
        self.update_information() to update the UI.
        """

        # If filters are passed in, assign them to self.filters first
        if filters:
            print "Apply() got filters."
            self.filters = filters

        # Convert current filters to native datatypes
        native_filters = self.convert_filters_to_native_types()

        # Get filtered data from data processor
        filtered_dataset = self.data_processor.filter_dataset(self.dataset, native_filters)

        # Check that it isn't None
        try:
            if not filtered_dataset.empty:
                self.dataset = filtered_dataset
                self.current_row_index = 0
                self.update_information()
        except AttributeError as e:
            print "ERROR: Empty dataset after filtering: %s" % e.message

            if not filtered_dataset:
                tkMessageBox.showerror("Error", "Selected filters returned empty dataset. Filters have been reset to defaults.")
                self.filters = self.get_default_filters()
            else:
                print "BUG: Reached unreachable code (TINTagger.apply_filters())"

    def convert_filters_to_native_types(self):
        """
        Converts Tkinter-specific filetypes present in the current filter (e.g. tk.IntVar(), tk.StringVar()), into
        native python datatypes (e.g. int, str) and returns the converted filters.
        """
        native = {
            "included_counts": int(self.filters["included_counts"].get()),
            "excluded_counts": int(self.filters["excluded_counts"].get()),
            "psi": float(self.filters["psi"].get()),
            "rpkm": float(self.filters["rpkm"].get())
        }

        splice_types = {}
        for name, info in self.filters["splice_type"].items():
            splice_types[name] = {
                "description": self.filters["splice_type"][name]["description"],
                "enabled": self.filters["splice_type"][name]["enabled"],
                "enabled_var": int(self.filters["splice_type"][name]["enabled_var"].get())
            }

        native["splice_type"] = splice_types

        return native

    def convert_filters_to_tk_specific_datatypes(self, filters):
        """
        Converts certain dataset filter entries to tk-specific datatypes (e.g. tk.IntVar(), tk.StringVar()) from their
        native python datatype representation. Returns the converted filters.
        """

        # Initialize variables as tk-specific datatypes
        tk_specific = {
            "included_counts": tk.StringVar(),
            "excluded_counts": tk.StringVar(),
            "psi": tk.StringVar(),
            "rpkm": tk.StringVar()
        }

        # Assign values to the variables, based on the values from the native filters
        tk_specific["included_counts"].set(str(filters["included_counts"]))
        tk_specific["excluded_counts"].set(str(filters["excluded_counts"]))
        tk_specific["psi"].set(str(filters["psi"]))
        tk_specific["rpkm"].set(str(filters["rpkm"]))

        splice_types = {}
        for name, info in filters["splice_type"].items():
            splice_types[name] = {
                "description": filters["splice_type"][name]["description"],
                "enabled": filters["splice_type"][name]["enabled"],
                "enabled_var": tk.IntVar()
            }
            splice_types[name]["enabled_var"].set(splice_types[name]["enabled"])

        tk_specific["splice_type"] = splice_types

        return tk_specific

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

    def open_file(self):
        """
        Presents a filedialog for the user to choose a file. If a file is selected, the dataset is read
        into a Pandas.DataFrame
        """
        # Open file chooser dialog
        filename = askopenfilename()

        # Check if file was chosen
        if len(filename) > 0:
            self.read_dataset(filename)

    def read_dataset(self, filepath):
        """
        Calls the data processor to read file by path.
        """
        # Feedback: Show busy cursor while dataset is loaded
        self.config(cursor="wait")

        # Get dataset
        self.dataset = self.data_processor.load_dataset(filepath)

        # Default to the first line of the file
        self.current_row_index = 0

        # Handle data
        self.update_information()

        # Show normal cursor again
        self.config(cursor="")

    def update_information(self):
        """
        Handles top-level stuff for each event (row):
            - Collects row data
            - Orchestrates UI update
        """

        # If no dataset is loaded, prompt user to load one
        try:
            if len(self.dataset) > 0:
                self.load_frame.destroy()

                # Show right sidebar
                self.right_sidebar.grid()
                self.left_frame.grid()

        except TypeError:
            # Hide right sidebar
            self.right_sidebar.grid_remove()
            self.left_frame.grid_remove()

            self.set_statusbar_text("No dataset loaded.")
            self.load_frame.grid()

            # Dodge the open-file-dialog when testing/debugging
            if self.testing:
                self.read_dataset("/Users/jonas/Dropbox/phd/code/tin_tagger/datasets/mikes_query_with_tags.tsv")

            return

        # Make sure the dataset is not empty
        if len(self.dataset) == 0:
            # Hide right sidebar
            self.right_sidebar.grid_remove()
            # Hide sample frame
            self.left_frame.grid_remove()
            # Update statusbar text
            self.set_statusbar_text("Loaded dataset was empty, please select another.")
            self.load_frame.grid()
            return

        # Display which event we're at in the statusbar
        self.set_statusbar_text("Splicing event %d/%d" % (self.current_row_index + 1, len(self.dataset)))

        # Update tagging progress information in statusbar
        self.update_tag_information()

        # Get data for this row
        data = self.get_row_data()

        # Populate the sidebar with general information
        self.location_text["text"] = data["location"]
        self.splice_type_text["text"] = "%s (%s)" % (self.splice_type_map[data["splice_type"]], data["splice_type"])
        self.exons_text["text"] = data["exons"]
        self.sample_text["text"] = data["sample_of_interest"]
        self.gene_text["text"] = data["gene_symbol"]
        self.strand_text["text"] = data["strand"]
        self.psi_text["text"] = "%.2f" % data["exon_psi"]
        self.included_counts_text["text"] = str(data["included_counts"])
        self.excluded_counts_text["text"] = str(data["excluded_counts"])

        # Clear all canvases before drawing new ones
        self.clear_all_canvases()

        # Draw events
        if data["splice_type"] == "AT":
            self.draw_alternative_terminator_event(data)
        elif data["splice_type"] == "ES":
            #self.draw_exon_skipping_events(data)
            self.exon_skipping_new(data)
        elif data["splice_type"] == "AD":
            self.draw_alternative_donor_events(data)

    def update_tag_information(self):
        """
        Reads the number of tagged events in total and populates the statusbar with text about how many
        positive, negative, and neutral tags are in the dataset.
        """

        # Find events that are tagged as interesting
        positive_tags = len(self.dataset.loc[self.dataset["event_tag"] == TAG_INTERESTING])
        self.statusbar_text_interesting["text"] = "%d" % positive_tags

        # Find events that are tagged as not interesting
        negative_tags = len(self.dataset.loc[self.dataset["event_tag"] == TAG_NOT_INTERESTING])
        self.statusbar_text_not_interesting["text"] = "%d" % negative_tags

        # Find events that are tagged as uncertain
        neutral_tags = len(self.dataset.loc[self.dataset["event_tag"] == TAG_UNCERTAIN])
        self.statusbar_text_uncertain["text"] = "%d" % neutral_tags

        # Find how many events are tagged, in total
        total_events = len(self.dataset)
        total_tags = positive_tags + negative_tags + neutral_tags
        self.statusbar_text_progress["text"] = "%d/%d" % (total_tags, total_events)

    def get_row_data(self):
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
        # TODO: Sanitize everything
        # TODO: Find included_counts and excluded_counts for other exons than the one in question

        # Set waiting cursor
        self.config(cursor="wait")

        # Get all sample names
        sample_names = list(self.dataset["name"].unique())

        # Get information
        splice_type = self.dataset.iloc[self.current_row_index]["splice_type"]
        sample_name = self.dataset.iloc[self.current_row_index]["name"]
        as_id = self.dataset.iloc[self.current_row_index]["as_id"]
        psi = self.dataset.iloc[self.current_row_index]["psi"]
        gene_symbol = self.dataset.iloc[self.current_row_index]["symbol"]
        strand = self.dataset.iloc[self.current_row_index]["strand"]
        exons = self.dataset.iloc[self.current_row_index]["exons"]
        chrom = self.dataset.iloc[self.current_row_index]["chr"]
        splice_start = self.dataset.iloc[self.current_row_index]["first_exon_start"]
        splice_stop = self.dataset.iloc[self.current_row_index]["last_exon_stop"]
        prev_exon_start = self.dataset.iloc[self.current_row_index]["prev_exon_start"]  # NaN for AT/AP
        prev_exon_stop = self.dataset.iloc[self.current_row_index]["prev_exon_stop"]  # NaN for AT/AP
        next_exon_start = self.dataset.iloc[self.current_row_index]["next_exon_start"]  # NaN for AT/AP
        next_exon_stop = self.dataset.iloc[self.current_row_index]["next_exon_stop"]  # NaN for AT/AP
        # Handle negative strand start- and stop- coordinates
        if strand == "-":
            splice_start = self.dataset.iloc[self.current_row_index]["last_exon_stop"]  # NaN for AT/AP
            splice_stop = self.dataset.iloc[self.current_row_index]["first_exon_start"]  # NaN for AT/AP
            prev_exon_start = self.dataset.iloc[self.current_row_index]["prev_exon_stop"]  # NaN for AT/AP
            prev_exon_stop = self.dataset.iloc[self.current_row_index]["prev_exon_start"]  # NaN for AT/AP
            next_exon_start = self.dataset.iloc[self.current_row_index]["next_exon_stop"]  # NaN for AT/AP
            next_exon_stop = self.dataset.iloc[self.current_row_index]["next_exon_start"]  # NaN for AT/AP
        included_counts = self.dataset.iloc[self.current_row_index]["included_counts"]
        excluded_counts = self.dataset.iloc[self.current_row_index]["excluded_counts"]
        prev_exon_name = self.dataset.iloc[self.current_row_index]["exon1"]  # NaN for AT/AP
        next_exon_name = self.dataset.iloc[self.current_row_index]["exon2"]  # NaN for AT/AP
        # prev_exon_id = self.dataset.iloc[self.current_row_index]["start_ex"]  # NaN for AT/AP
        # next_exon_id = self.dataset.iloc[self.current_row_index]["end_ex"]  # NaN for AT/AP

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
            is_reported = self.data_processor.is_event_reported_in_sample(sample_name, as_id, self.dataset)
            samples_data["is_reported"] = is_reported
            # Default to sample not being tagged
            sample_tag = TAG_NO_TAG
            # Default to RPKM being 0 (in case it's not reported)
            gene_rpkm = 0

            if is_reported:
                gene_rpkm = self.data_processor.get_gene_rpkm_by_sample_name(s_name, gene_symbol, self.dataset)
                sample_tag = self.data_processor.get_sample_tag_by_as_id(s_name, as_id, self.dataset)

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
                    "coverage": self.data_processor.get_coverage_by_coordinates(upstream_exon_coords, self.bam_paths[s_name], self.testing)
                }
                sample_exons["upstream_exon"] = upstream_exon

                # Handle downstream exon
                downstream_exon_coords = str(chrom) + ":" + str(int(next_exon_start)) + "-" + str(int(next_exon_stop))
                downstream_exon = {
                    "exon_name": next_exon_name,
                    "coverage": self.data_processor.get_coverage_by_coordinates(downstream_exon_coords, self.bam_paths[s_name], self.testing)
                }
                sample_exons["downstream_exon"] = downstream_exon

                # Keep track of coverage values
                upstream_exon_coverages.append(upstream_exon["coverage"])
                downstream_exon_coverages.append(downstream_exon["coverage"])

            # Handle the main exon
            main_exon_coords = str(chrom) + ":" + str(int(splice_start)) + "-" + str(int(splice_stop))
            exon_of_interest = {
                "exon_name": exons,
                "coverage": self.data_processor.get_coverage_by_coordinates(main_exon_coords, self.bam_paths[s_name], self.testing),
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

        # Reset cursor to normal again
        self.config(cursor="")

        return row_data

    def save_file(self):
        print "Bleep, blop, saving file."

    def draw_alternative_terminator_event(self, data):
        pass

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

    def next_button_clicked(self):
        """
        Handles next-button presses: Update row index and initiate reading of new row.
        """

        self.current_row_index += 1
        if self.current_row_index > len(self.dataset):
            print "Wops, no more rows (reached end of dataset)"
            self.current_row_index -= 1
            self.set_statusbar_text("Reached end of dataset.")
        else:
            self.update_information()

    def previous_button_clicked(self):
        """
        Handles previous-button presses: Update row index and initiate reading of new row.
        """

        self.current_row_index -= 1
        if self.current_row_index < 0:
            print "Wops, no more rows (reached start of dataset)"
            self.current_row_index += 1
            self.set_statusbar_text("Reached beginning of dataset.")
        else:
            self.update_information()

    def random_button_clicked(self):
        """
        Handles random-button presses: Pick and display a random event.
        """

        self.current_row_index = random.randint(0, len(self.dataset))
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
            gene_rpkm_percent_of_max = (float(gene_rpkm) / float(data["max_gene_rpkm"])) * 100

            # Create frame for this sample
            frame_for_sample = ttk.Frame(self.sample_frame, padding=5)
            frame_for_sample.grid(column=0, row=row_number, sticky="NEWS")
            self.sample_frame.rowconfigure(row_number, weight=1)

            # Display RPKM info in statusbar on hover
            display_text = "%s gene RPKM: %s/%s (%.0f%%)" % (
            sample_name, str(gene_rpkm), str(data["max_gene_rpkm"]), gene_rpkm_percent_of_max)
            frame_for_sample.bind("<Enter>", store_rpkm_text(display_text))

            # Create text for sample name
            sample_text = ttk.Label(frame_for_sample, text=sample_name, font="TkDefaultFont", anchor=tk.CENTER)
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
                remainder_frame = tk.Frame(rpkm_frame, width=10, background=COLOR_RED)
            remainder_frame.grid(column=0, row=0, sticky="NEWS")
            # Fill for RPKM
            fill_frame = tk.Frame(rpkm_frame, bg=COLOR_RED)
            fill_frame.grid(column=0, row=1, sticky="NEWS")
            rpkm_frame.rowconfigure(0, weight=100 - int(gene_rpkm_percent_of_max))
            rpkm_frame.rowconfigure(1, weight=int(gene_rpkm_percent_of_max))

            row_number += 1

    def draw_exon_skipping_events(self, data):
        """
        Draw exons using canvas instead of frame (one canvas per row, so x, y coordinates are the same for each sample
        """

        # Fill in samples information column
        self.populate_samples_frame(data)

        # Get strand information to know which way to draw exons
        strand = data["strand"]

        # Keep track of row number
        row_number = 0

        # First, draw a top frame containing exon names
        exons = sorted(data["samples"].values()[0]["exons"].keys(), key=natural_sort_key)
        if strand == "-":
            exons.reverse()
        exon_name_frame = ttk.Frame(self.exon_frame)
        exon_name_frame.grid(column=0, row=row_number, sticky="NEW")
        # Get names of exons
        exon_column = 0
        for exon_name in exons:
            # Create label for this exon
            exon_font = "tkDefaultFont 16"
            if exon_name == data["exons"]:
                exon_font = "tkDefaultFont 16 bold"
            exon_text = exon_name
            if len(exon_name) > 15:
                exon_text = exon_name[:14] + ".."
            exon_label = ttk.Label(exon_name_frame, text=exon_text, anchor=tk.CENTER, font=exon_font)
            exon_label.grid(row=0, column=exon_column, sticky="NEWS")
            exon_name_frame.columnconfigure(exon_column, weight=1)
            exon_column += 1

        # TEST: Horizontal separator
        sep = ttk.Separator(self.exon_frame, orient=tk.HORIZONTAL)
        sep.grid(column=0, row=row_number, sticky="EWS")
        # END TEST

        row_number += 1

        # Draw exons for each sample
        for sample_name in sorted(data["samples"].keys(), key=natural_sort_key):

            # Get data for this sample
            sample_data = data["samples"][sample_name]
            # Get sample exons
            sample_exons = sorted(sample_data["exons"].keys(), key=natural_sort_key)
            # Reverse exons if we're on the minus trand
            if strand == "-":
                sample_exons.reverse()

            # Dimension variables
            canvas_width = 300
            canvas_height = 100
            width_per_exon_container = canvas_width / 3
            x_offset = 0  # This will increase with width_per_exon_container for each exon drawn
            exon_width = 60
            exon_height = 80
            exon_start_y = (canvas_height - exon_height) / 2

            # Highlight background for sample of interest
            canvas_background = COLOR_WHITE
            if sample_name == data["sample_of_interest"]:
                canvas_background = COLOR_LIGHTRED

            # Initialize canvas and grid it to this row in the center frame.
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")
            self.exon_frame.columnconfigure(0, weight=1)

            # Keep track of canvases used
            self.canvases.append(row_canvas)

            # Set even weight for every row in the center-frame
            self.exon_frame.rowconfigure(row_number, weight=1)

            # exon_number = 0  # Do I not need this?
            for exon_name in sample_exons:
                exon_data = sample_data["exons"][exon_name]
                coverage = exon_data["coverage"]
                max_coverage = exon_data["max_coverage"]
                percent_of_max_coverage = (float(coverage) / float(max_coverage)) * 100

                # Update exon start position
                exon_start_x = x_offset + ((width_per_exon_container - exon_width) / 2)  # /2 because there's padding at both side

                # TEST: Highlight exon of interest
                rect_width = 1
                rect_dash = None
                outline_color = COLOR_DARKBLUE
                if exon_name == data["exons"] and sample_name == data["sample_of_interest"]:
                    rect_width = 3
                    #rect_dash = (5, 5)
                    outline_color = COLOR_DARKBLUE

                # Draw background rectangle
                row_canvas.create_rectangle(exon_start_x, exon_start_y, exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, outline=outline_color, width=rect_width, dash=rect_dash)

                # Draw fill rectangle based on relative coverage
                fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
                fill_end_y = exon_start_y + exon_height
                row_canvas.create_rectangle(exon_start_x, fill_start_y, exon_start_x + exon_width, fill_end_y, fill=COLOR_BLUE, outline=outline_color, width=rect_width, dash=rect_dash)

                # Draw coverage text (text shadow first)
                text_start_x = exon_start_x + (exon_width / 2)
                text_start_y = exon_height - 10
                text_color = COLOR_WHITE
                #row_canvas.create_text(text_start_x + 1, text_start_y + 1, anchor=tk.CENTER, text=str(coverage), font="tkDefaultFont 16", fill=shadow_color)
                row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(coverage), tags="text_shadow", font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW)

                row_canvas.create_text(text_start_x, text_start_y, anchor=tk.CENTER, text=str(coverage), font="tkDefaultFont 16", fill=text_color)

                # Increment x_offset
                x_offset += width_per_exon_container

            # Increment row number
            row_number += 1

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

        # First, draw a top frame containing exon names
        exon_name_frame = ttk.Frame(self.exon_frame)
        exon_name_frame.grid(column=0, row=row_number, sticky="NEWS")

        row_number += 1

        # Add upstream exon name
        upstream_exon_name = data["samples"].values()[0]["exons"]["upstream_exon"]["exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."

        upstream_label = ttk.Label(exon_name_frame, text=upstream_exon_name, font="tkDefaultFont 16", anchor=tk.E)
        upstream_label.grid(column=0, row=0, sticky="NEWS")

        # Add exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."
        exon_of_interest_label = ttk.Label(exon_name_frame, text=exon_name, font="tkDefaultFont 16 bold", anchor=tk.CENTER)
        exon_of_interest_label.grid(column=1, row=0, sticky="NEWS")

        # Add downstream exon name
        downstream_exon_name = data["samples"].values()[0]["exons"]["downstream_exon"]["exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + ".."

        downstream_label = ttk.Label(exon_name_frame, text=downstream_exon_name, font="tkDefaultFont 16", anchor=tk.W)
        downstream_label.grid(column=2, row=0, sticky="NEWS")

        # Finally, add even weight for all exon names
        exon_name_frame.columnconfigure(0, weight=1)
        exon_name_frame.columnconfigure(1, weight=1)
        exon_name_frame.columnconfigure(2, weight=1)

        ############################################
        ################ Draw exons ################
        ############################################

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        donor_site_width = (canvas_width / 3) * 2
        upstream_exon_width = 60
        main_exon_width = 30
        downstream_exon_width = 60
        exon_height = 80
        exon_start_y = (canvas_height - exon_height) / 2

        # Iterate samples
        for sample_name in sorted(data["samples"].keys(), key=natural_sort_key):
            sample_data = data["samples"][sample_name]

            # Get sample exons
            upstream_exon = sample_data["exons"]["upstream_exon"]
            exon_of_interest = sample_data["exons"]["exon_of_interest"]
            downstream_exon = sample_data["exons"]["downstream_exon"]

            # Highlight background for sample of interest
            canvas_background = COLOR_WHITE
            if sample_name == data["sample_of_interest"]:
                canvas_background = COLOR_SAMPLE_HIGHLIGHT

            # Initialize canvas and grid it to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")

            # Keep track of canvases uses
            self.canvases.append(row_canvas)

            # Set even weight for every row in the center frame
            self.exon_frame.rowconfigure(row_number, weight=1)

            ######################
            # Draw upstream exon #
            ######################
            upstream_coverage = upstream_exon["coverage"]
            max_coverage = data["max_upstream_exon_coverage"]
            percent_of_max_coverage = (float(upstream_coverage) / float(max_coverage)) * 100

            upstream_exon_start_x = (donor_site_width - upstream_exon_width) / 2

            # Draw exon background
            row_canvas.create_rectangle(upstream_exon_start_x, exon_start_y, upstream_exon_start_x + upstream_exon_width, exon_start_y + exon_height, fill=canvas_background)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(upstream_exon_start_x, fill_start_y, upstream_exon_start_x + upstream_exon_width, fill_end_y, fill=COLOR_BLUE, outline=COLOR_DARKBLUE)

            # Draw coverage text
            text_start_x = upstream_exon_start_x + (upstream_exon_width / 2)
            text_start_y = exon_height - 10
            row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(upstream_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW, tags="text_shadow")
            row_canvas.create_text(text_start_x, text_start_y, text=str(upstream_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT)

            ##################
            # Draw main exon #
            ##################
            main_exon_coverage = exon_of_interest["coverage"]
            main_exon_max_coverage = data["max_exon_of_interest_coverage"]
            percent_of_max_coverage = (float(main_exon_coverage) / float(main_exon_max_coverage)) * 100

            main_exon_start_x = upstream_exon_start_x + upstream_exon_width

            border_width = 1
            if sample_name == data["sample_of_interest"]:
                border_width = 3

            # Draw exon background
            row_canvas.create_rectangle(main_exon_start_x, exon_start_y, main_exon_start_x + main_exon_width, exon_start_y + exon_height, fill=canvas_background, width=border_width)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
            row_canvas.create_rectangle(main_exon_start_x, fill_start_y, main_exon_start_x + main_exon_width, fill_end_y, fill=COLOR_BLUE, outline=COLOR_DARKBLUE, width=border_width)

            # Draw coverage text
            text_start_x = main_exon_start_x + (main_exon_width / 2)
            row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(main_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW, tags="text_shadow")
            row_canvas.create_text(text_start_x, text_start_y, text=str(main_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT)

            ########################
            # Draw downstream exon #
            ########################
            downstream_exon_coverage = downstream_exon["coverage"]
            downstream_exon_max_coverage = data["max_downstream_exon_coverage"]
            percent_of_max_coverage = (float(downstream_exon_coverage) / float(downstream_exon_max_coverage)) * 100

            downstream_exon_start_x = donor_site_width  # Not sure if this is right

            # Draw exon background
            row_canvas.create_rectangle(downstream_exon_start_x, exon_start_y, downstream_exon_start_x + downstream_exon_width, exon_start_y + exon_height, fill=canvas_background)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
            row_canvas.create_rectangle(downstream_exon_start_x, fill_start_y, downstream_exon_start_x + downstream_exon_width, fill_end_y, fill=COLOR_BLUE, outline=COLOR_DARKBLUE)

            # Draw coverage text
            text_start_x = downstream_exon_start_x + (downstream_exon_width / 2)
            row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(downstream_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW, tags="text_shadow")
            row_canvas.create_text(text_start_x, text_start_y, text=str(downstream_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT)

            # Prepare for next sample
            row_number += 1

        # Expand exon frame horizontally
        self.exon_frame.columnconfigure(0, weight=1)

    def exon_skipping_new(self, data):
        """
        Draw exon skipping events
        """
        ######################################
        # Fill in samples information column #
        ######################################
        self.populate_samples_frame(data)

        ################################
        # Draw exon names in top frame #
        ################################

        # Keep track of row number
        row_number = 0

        # Draw exon names on top
        exon_name_frame = ttk.Frame(self.exon_frame)
        exon_name_frame.grid(column=0, row=row_number, sticky="NEWS")

        self.exon_frame.columnconfigure(0, weight=1)

        row_number += 1

        # Add upstream exon name
        upstream_exon_name = data["samples"].values()[0]["exons"]["upstream_exon"]["exon_name"]
        if len(upstream_exon_name) > 15:
            upstream_exon_name = upstream_exon_name[:14] + ".."

        upstream_label = ttk.Label(exon_name_frame, text=upstream_exon_name, font="tkDefaultFont 16", anchor=tk.CENTER)
        upstream_label.grid(column=0, row=0, sticky="NEWS")

        # Add exon of interest name
        exon_name = data["exons"]
        if len(exon_name) > 15:
            exon_name = exon_name[:14] + ".."
        exon_of_interest_label = ttk.Label(exon_name_frame, text=exon_name, font="tkDefaultFont 16 bold", anchor=tk.CENTER)
        exon_of_interest_label.grid(column=1, row=0, sticky="NEWS")

        # Add downstream exon name
        downstream_exon_name = data["samples"].values()[0]["exons"]["downstream_exon"]["exon_name"]
        if len(downstream_exon_name) > 15:
            downstream_exon_name = downstream_exon_name[:14] + ".."

        downstream_label = ttk.Label(exon_name_frame, text=downstream_exon_name, font="tkDefaultFont 16", anchor=tk.CENTER)
        downstream_label.grid(column=2, row=0, sticky="NEWS")

        # Finally, add even weight for all exon names
        exon_name_frame.columnconfigure(0, weight=1)
        exon_name_frame.columnconfigure(1, weight=1)
        exon_name_frame.columnconfigure(2, weight=1)

        ##############
        # Draw exons #
        ##############

        # Dimension variables
        canvas_width = 300
        canvas_height = 100
        width_per_exon_container = canvas_width / 3
        exon_width = 60
        exon_height = 80
        exon_start_y = (canvas_height - exon_height) / 2

        # Iterate samples
        for sample_name in sorted(data["samples"].keys(), key=natural_sort_key):
            sample_data = data["samples"][sample_name]

            # Get sample exons
            upstream_exon = sample_data["exons"]["upstream_exon"]
            main_exon = sample_data["exons"]["exon_of_interest"]
            downstream_exon = sample_data["exons"]["downstream_exon"]

            # Highlight background for sample of interest
            canvas_background = "white"
            if sample_name == data["sample_of_interest"]:
                canvas_background = COLOR_SAMPLE_HIGHLIGHT

            # Initialize canvas and grid to this row in the exon frame
            row_canvas = ResizingCanvas(self.exon_frame, bg=canvas_background, highlightthickness=0, width=canvas_width, height=canvas_height)
            row_canvas.grid(row=row_number, column=0, sticky="NEWS")

            # TEST: Create button frame
            # TODO: If this is tagged, highlight the corresponding button
            style = ttk.Style()
            style.configure(STYLE_BUTTON_INTERESTING, foreground=COLOR_DARKBLUE)
            style.configure(STYLE_BUTTON_NOT_INTERESTING, foreground=COLOR_DARKBLUE)
            style.configure(STYLE_BUTTON_UNCERTAIN, foreground=COLOR_DARKBLUE, font="tkDefaultFont 16 bold")
            style.map(
                STYLE_BUTTON_INTERESTING,
                foreground=[("active", COLOR_INTERESTING)],
            )
            style.map(
                STYLE_BUTTON_NOT_INTERESTING,
                foreground=[("active", COLOR_NOT_INTERESTING)]
            )
            style.map(
                STYLE_BUTTON_UNCERTAIN,
                foreground=[("active", COLOR_UNCERTAIN)]
            )
            test = random.randint(0, 2)
            bg_color = COLOR_INTERESTING
            if test == 1:
                bg_color = COLOR_NOT_INTERESTING
            if test == 0:
                bg_color = COLOR_UNCERTAIN
            button_frame = tk.Frame(self.exon_frame, bg=bg_color, padx=5, pady=0, borderwidth=0, relief=tk.SOLID)
            button_frame.grid(row=row_number, column=1, sticky="NEWS")
            up_button = ttk.Button(button_frame, text=u"\u25B2", style=STYLE_BUTTON_INTERESTING)
            up_button.grid(column=0, row=0, sticky="NEWS")
            down_button = ttk.Button(button_frame, text=u"\u25BC", style=STYLE_BUTTON_NOT_INTERESTING)
            down_button.grid(column=0, row=1, sticky="NEWS")
            uncertain_button = ttk.Button(button_frame, text="?", style=STYLE_BUTTON_UNCERTAIN)
            uncertain_button.grid(column=0, row=2, sticky="NEWS")
            button_frame.rowconfigure(0, weight=1)
            button_frame.rowconfigure(1, weight=1)
            button_frame.rowconfigure(2, weight=1)
            # END TEST

            # Keep track of canvases used
            self.canvases.append(row_canvas)

            # Set even weight for every row in the exon frame
            self.exon_frame.rowconfigure(row_number, weight=1)

            ######################
            # Draw upstream exon #
            ######################
            upstream_exon_coverage = upstream_exon["coverage"]
            upstream_exon_max_coverage = data["max_upstream_exon_coverage"]
            percent_of_max_coverage = (float(upstream_exon_coverage) / float(upstream_exon_max_coverage)) * 100

            upstream_exon_start_x = (width_per_exon_container - exon_width) / 2

            # Draw exon background
            row_canvas.create_rectangle(upstream_exon_start_x, exon_start_y, upstream_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
            fill_end_y = exon_start_y + exon_height
            row_canvas.create_rectangle(upstream_exon_start_x, fill_start_y, upstream_exon_start_x + exon_width, fill_end_y, fill=COLOR_BLUE, outline=COLOR_DARKBLUE)

            # Draw coverage text
            text_start_x = upstream_exon_start_x + (exon_width / 2)
            text_start_y = exon_height - 10
            row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(upstream_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW, tags="text_shadow")
            row_canvas.create_text(text_start_x, text_start_y, text=str(upstream_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT)

            ##################
            # Draw main exon #
            ##################
            main_exon_coverage = main_exon["coverage"]
            main_exon_max_coverage = data["max_exon_of_interest_coverage"]
            percent_of_max_coverage = (float(main_exon_coverage) / float(main_exon_max_coverage)) * 100

            main_exon_start_x = upstream_exon_start_x + width_per_exon_container  # Not sure if this is correct

            border_width = 1
            if sample_name == data["sample_of_interest"]:
                border_width = 3

            # Draw exon background
            row_canvas.create_rectangle(main_exon_start_x, exon_start_y, main_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background, width=border_width)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
            row_canvas.create_rectangle(main_exon_start_x, fill_start_y, main_exon_start_x + exon_width, fill_end_y, fill=COLOR_BLUE, outline=COLOR_DARKBLUE, width=border_width)

            # Draw coverage text
            text_start_x = main_exon_start_x + (exon_width / 2)
            row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(main_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW, tags="text_shadow")
            row_canvas.create_text(text_start_x, text_start_y, text=str(main_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT)

            ########################
            # Draw downstream exon #
            ########################
            downstream_exon_coverage = downstream_exon["coverage"]
            downstream_exon_max_coverage = data["max_downstream_exon_coverage"]
            percent_of_max_coverage = (float(downstream_exon_coverage) / float(downstream_exon_max_coverage)) * 100

            downstream_exon_start_x = main_exon_start_x + width_per_exon_container

            # Draw exon background
            row_canvas.create_rectangle(downstream_exon_start_x, exon_start_y, downstream_exon_start_x + exon_width, exon_start_y + exon_height, fill=canvas_background)

            # Draw exon fill
            fill_start_y = (exon_start_y + exon_height) - (int((percent_of_max_coverage/100) * exon_height))
            row_canvas.create_rectangle(downstream_exon_start_x, fill_start_y, downstream_exon_start_x + exon_width, fill_end_y, fill=COLOR_BLUE, outline=COLOR_DARKBLUE)

            # Draw coverage text
            text_start_x = downstream_exon_start_x + (exon_width / 2)
            row_canvas.create_text(text_start_x + 1, text_start_y + 1, text=str(downstream_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT_SHADOW, tags="text_shadow")
            row_canvas.create_text(text_start_x, text_start_y, text=str(downstream_exon_coverage), font="tkDefaultFont 16", fill=COLOR_CANVAS_TEXT)

            # Update row index for next sample
            row_number += 1


if __name__ == "__main__":
    # Create and run app
    app = TINTagger()
    app.wm_title("TIN-Tagger")
    screen_width = app.winfo_screenwidth()
    screen_height = app.winfo_screenheight()
    print "Screen width / height: %d / %d" % (screen_width, screen_height)
    # Fill 80% of screen
    new_width = int((float(screen_width) / 100.00) * 80)
    new_height = int((float(screen_height) / 100.00) * 80)
    print "New width / height: %d / %d" % (new_width, new_height)
    app.geometry("%dx%d" % (new_width, new_height))
    #app.withdraw()
    app.center_window(app)
    #app.update()
    #app.deiconify()
    #app.geometry("1024x576")

    app.mainloop()


