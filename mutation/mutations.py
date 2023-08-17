import string, math

from functools import cache

import Bio
from Bio import Graphics, Phylo
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Mutations:
    """ Get mutation info from an alignment """

    def __init__(self, alignment):
        """ Initialize the Mutations object """

        self.alignment = alignment
        self.mutations: dict[dict[int: list]] = {}
        self.reference: int = 0

    def list_mutations(self, *, reference: int|str=0, apobec: bool=False, g_to_a: bool=False) -> dict[int: list]:
        """ Get mutations from a sequence and a reference sequence """

        reference_str: str = ""
        mutations: list[dict[str: list]] = []

        if isinstance(reference, str):
            reference = self.get_seq_index_by_id(reference)
        elif not isinstance(reference, int):
            raise TypeError(f"Expected reference to be an int or str, got {type(reference)}")
        
        self.reference = reference

        if isinstance(self.alignment[reference], SeqRecord):
            reference_str = str(self.alignment[reference].seq)
        elif isinstance(self.alignment[reference], Seq):
            reference_str = str(self.alignment[reference])
        elif isinstance(self.alignment[reference], str):
            reference_str = self.alignment[reference]

        if reference_str == "":
            raise ValueError(f"Reference sequence {reference} is empty")
        
        for sequence in self.alignment:
            if isinstance(sequence, SeqRecord):
                sequence_str = str(sequence.seq)
            elif isinstance(sequence, Seq):
                sequence_str = str(sequence)
            elif isinstance(sequence, str):
                sequence_str = sequence

            mutations.append(self.get_mutations(sequence=sequence_str, reference=reference_str, apobec=apobec, g_to_a=g_to_a))

        return mutations
        
    def get_seq_index_by_id(self, id: str) -> str:
        """ Get a sequence from the alignment by its id """

        for index, sequence in enumerate(self.alignment):
            if sequence.id == id:
                return index
        
        raise IndexError(f"Could not find sequence with id {id}")

    @staticmethod
    def get_mutations(*, sequence: str|Seq|SeqRecord, reference: str|Seq|SeqRecord, apobec: bool=False, g_to_a: bool=False) -> dict[int: list]:
        """ Get mutations from a a sequence and a reference sequence 
        returns a dictionary of mutations where the key is the position of the mutation and the value is a list of types of mutations """

        if type(sequence) is Seq:
                sequence = str(sequence)
        elif type(sequence) is  SeqRecord:
                sequence = str(sequence.seq)
        elif type(sequence) is not str:
                raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")

        if type(reference) is Seq:
                reference = str(reference)
        elif type(sequence) is  SeqRecord:
                reference = str(reference.seq)
        elif type(reference) is not str:
                raise TypeError(f"Expected reference to be a string, Seq, or SeqRecord, got {type(reference)}")
            
        if len(sequence) != len(reference):
            raise ValueError("Reference and sequence must be the same length")
        
        return Mutations.get_mutations_from_str(sequence=sequence, reference=reference, apobec=apobec, g_to_a=g_to_a)
        
    @cache
    @staticmethod
    def get_mutations_from_str(*, sequence: str, reference: str, apobec: bool, g_to_a: bool) -> dict[int: list]:
        """ Get mutations from a sequence and a reference sequence
        separated out so it can be cached (Seq and SeqRecord are not hashable) """

        mutations: dict = {}

        if sequence == reference:
            return mutations

        for base_index in range(len(sequence)):
            if reference[base_index] != sequence[base_index]:
                mutations[base_index] = []
                
                if sequence[base_index] != "-":
                    mutations[base_index].append(sequence[base_index])
                else:
                    mutations[base_index].append("Gap")

                if reference[base_index] == "G" and sequence[base_index] == "A":
                    if g_to_a:
                        mutations[base_index].append("G->A mutation")

                    if apobec and base_index < len(sequence)-3 and sequence[base_index+1] in "AG" and sequence[base_index+2] != "C":
                        mutations[base_index].append("APOBEC")
        
        return mutations
    
AlignInfo.Mutations = Mutations


from reportlab.lib import colors
from reportlab.lib.colors import Color
from reportlab.lib.units import inch

from reportlab.pdfbase.pdfmetrics import stringWidth

from reportlab.graphics.shapes import Drawing, String, Line, Rect, Circle, PolyLine

from Bio.Graphics import _write
from Bio.Align import AlignInfo


class MutationPlot:
    """ Create and output a mutation plot """

    def __init__(self, alignment, *, tree: str|object=None, output_format: str="svg", seq_name_font: str="Helvetica", seq_name_font_size: int=20, left_margin: float=.25, top_margin: float=.25, botom_margin: float=0, right_margin: float=0, mark_reference: bool=True, title: str=None, title_font="Helvetica", title_font_size: int=30, ruler: bool=True, ruler_font: str="Helvetica", ruler_font_size: int=15, ruler_major_ticks: int=10, ruler_minor_ticks=3):
        """ Initialize the MutationPlot object """

        self.alignment = alignment

        if tree is not None:
            if isinstance(tree, Bio.Phylo.BaseTree.Tree):
                self.tree = tree
            else:
                raise TypeError("tree must be a Bio.Phylo.BaseTree.Tree object (or a derivative)")

        self._mutations = AlignInfo.Mutations(alignment)
        self._seq_count = len(alignment)
        self._seq_length = len(alignment[0])
        self.mark_reference: bool = mark_reference

        self.output_format: str = output_format
        self.seq_name_font: str = seq_name_font
        self.seq_name_font_size: int = seq_name_font_size

        self.top_margin: float = inch * top_margin
        self.bottom_margin: float = inch * botom_margin
        self.left_margin: float = inch * left_margin
        self.right_margin: float = inch * right_margin

        self.title: str = title
        self.title_font: str = title_font
        self.title_font_size: int = title_font_size
        self._title_font_height: float = self._font_height(self.title_font, self.title_font_size)
        self._title_height = 0

        self.ruler: bool = ruler
        self.ruler_font: str = ruler_font
        self.ruler_font_size: int = ruler_font_size
        self.ruler_major_ticks: int = ruler_major_ticks
        self.ruler_minor_ticks: int = ruler_minor_ticks
        self._ruler_font_height: float = self._font_height(self.ruler_font, self.ruler_font_size)
        self._ruler_height = 0 if not ruler else self._ruler_font_height * 3

        self._plot_width: float = 14 * inch
        self._plot_floor: float = self.bottom_margin + self._ruler_height

        self._seq_name_width: float = self._max_seq_name_width
        self._width: float = self.left_margin + self._plot_width + (inch/4) + self._seq_name_width + self.right_margin
        
        self._seq_height: float = self._font_height(self.seq_name_font, self.seq_name_font_size)
        self._seq_gap: float = self._seq_height / 5
        self._height: float = len(self.alignment) * (self._seq_height + self._seq_gap) + self.top_margin + self.bottom_margin + self._title_height + self._ruler_height

        self._plot_colors: dict[str: str] = {"A": "#42FF00", "C": "#41B8EE", "G": "#FFA500", "T": "#EE0B10", "Gap": "#666666"}

    def draw(self, output_file, reference: str|int=0, apobec: bool=False, g_to_a: bool=False, sort: str="similar", narrow_markers: bool=True, min_marker_width: float=1):
        """ Writes out the mutation plot to a file """
        
        drawing = self.drawing = Drawing(self._width, self._height)
        self.mutations_list = self._mutations.list_mutations(reference=reference, apobec=apobec, g_to_a=g_to_a)
        self.narrow_markers: bool = narrow_markers
        self.min_marker_width: float = min_marker_width

        if sort == "similar":
            sorted_keys = self._sort_similar()
        
        elif sort == "tree":
            if self.tree is None:
                raise ValueError("Cannot sort by tree if no tree is provided")
            
            sorted_keys = self._indexes_by_tree_order()

        else: 
            sorted_keys = range(len(self.mutations_list))

        self._draw_ruler()

        # for seq_index, mutations in enumerate(self.mutations_list):
        for plot_index, seq_index in enumerate(sorted_keys):
            mutations = self.mutations_list[seq_index]

            # Add label for sequence
            id = self.alignment[seq_index].id
            if self.mark_reference and seq_index == self._mutations.reference:
                id += " (r)"

            x: float = self.left_margin + self._plot_width + (inch/4)
            y: float = ((self._seq_count-(plot_index + .5)) * (self._seq_height + self._seq_gap)) + self._plot_floor
            sequence_str: String = String(x, y, id, fontName="Helvetica", fontSize=self.seq_name_font_size)
            drawing.add(sequence_str)

            # Add base line for sequence
            x1: float = self.left_margin
            x2: float = self.left_margin + self._plot_width
            y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self._seq_gap) + self._seq_gap + self._plot_floor
            sequence_baseline: Line = Line(x1, y, x2, y, strokeColor=colors.lightgrey)
            drawing.add(sequence_baseline)

            self._draw_mutations(plot_index, mutations)

        return _write(drawing, output_file, self.output_format)
    
    def _draw_mutations(self, plot_index: int, mutations: dict[int: list]) -> None:
        """ Draw mutations for a sequence """

        for base, mutation in mutations.items():
            for code in mutation:
                if code in self._plot_colors:
                    x1: float = self.left_margin + self._base_left(base)
                    x2: float = self.left_margin + self._base_left(base+1)

                    y1: float = ((self._seq_count-plot_index) * (self._seq_height + self._seq_gap)) + (self._seq_gap/2) + self._plot_floor
                    y2: float = ((self._seq_count-(plot_index+1)) * (self._seq_height + self._seq_gap)) + self._seq_gap + self._plot_floor

                    if code != "Gap" and self.narrow_markers and x2-x1 > self.min_marker_width:
                        x1, x2 = (
                            x1+(((x2-x1)-self.min_marker_width)/2), 
                            x2-(((x2-x1)-self.min_marker_width)/2)
                        )

                    base_color: Color = self._hex_to_color(self._plot_colors[code])
                    base_mark: Rect = Rect(x1, y1, x2-x1, y2-y1, fillColor=base_color, strokeColor=base_color, strokeWidth=0.1)
                    self.drawing.add(base_mark)
        
        # APOBEC and G->A go second so they go on top of other elements
        for base, mutation in mutations.items():
            if "APOBEC" in mutation:
                x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
                y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self._seq_gap) + self._seq_gap + self._plot_floor
                
                circle = Circle(x, y, (self._seq_height/3)/2, fillColor=self._hex_to_color("#FF00FF"), strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=0.1)
                
                self.drawing.add(circle)

            elif "G->A mutation" in mutation:
                x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
                y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self._seq_gap) + self._seq_gap + self._plot_floor
                
                diamond = self._g_to_a_diamond(x, y)
                
                self.drawing.add(diamond)

    def _draw_ruler(self) -> None:
        """ Draw the ruler at the bottom of the plot """

        label_width = stringWidth(str(self._seq_length), self.seq_name_font, self.seq_name_font_size)
        marks: list = []

        if self._seq_length <= 20:
            self._ruler_marks(range(self._seq_length))
        else:
            bases = [0]
            last_base: int = self.significant_digits(self._seq_length)-1
            
            dist = last_base / self.ruler_major_ticks

            for base in range(1, self.ruler_major_ticks+1):
                bases.append(int(base * dist))

            bases.append(last_base)

            self._ruler_marks(bases)

        # Draw vertical line
        x1: float = self.left_margin
        x2: float = self.left_margin + self._plot_width
        y: float = self.bottom_margin + (self._ruler_font_height * 3)

        ruler_line: Line = Line(x1, y, x2, y, strokeColor=colors.black, strokeWidth=2)
        self.drawing.add(ruler_line)

    def _ruler_marks(self, marked_bases: list) -> None:
        """ Draw marks on the ruler """

        mark_locations: list = []

        for index, base in enumerate(marked_bases):
                mark_locations.append(self._ruler_label(base))
                self._ruler_heavy_tick(base)

                if index:
                    spacing = (mark_locations[index][0] - mark_locations[index-1][0]) / (self.ruler_minor_ticks+1)
                    left = mark_locations[index-1][0]

                    for tick in range(1, self.ruler_minor_ticks+1):
                        self._ruler_light_tick(left + (tick * spacing))
                    

    def _ruler_label(self, base: int) -> tuple[float]:
        """ Draw a label on the ruler
         returns the coordinates of the label """

        x: float = self.left_margin + self._base_center(base)
        y: float = self.bottom_margin+self._ruler_font_height

        self.drawing.add(String(x, y, str(base + 1), textAnchor="middle", fontName=self.ruler_font, fontSize=self.ruler_font_size))

        return (x, y)

    def _ruler_heavy_tick(self, base: int) -> tuple[float]:
        """ Draw a heavy tick on the ruler
        returns the x, top, and bottom of the tick """

        x: float = self.left_margin + self._base_center(base)
        top: float = self.bottom_margin + (self._ruler_font_height * 3)
        bottom: float = self.bottom_margin + (self._ruler_font_height * 2)

        self.drawing.add(Line(x, top, x, bottom, strokeColor=colors.black, strokeWidth=1))

        return (x, top, bottom)

    def _ruler_light_tick(self, x: float) -> tuple[float]:
        """ Draw a light tick on the ruler 
        x is the raw x, including the _left_margin
        returns the x, top, and bottom of the tick"""

        top: float = self.bottom_margin + (self._ruler_font_height * 3)
        bottom: float = self.bottom_margin + (self._ruler_font_height * 2.5)

        self.drawing.add(Line(x, top, x, bottom, strokeColor=colors.black, strokeWidth=1))

        return (x, top, bottom)

    def _base_left(self, base: int) -> float:
        """ Get the left coordinate of a base """

        if base == self._seq_length:
            return self._plot_width
        
        return (base / self._seq_length) * self._plot_width
    
    def _base_center(self, base: int) -> float:
        """ Get the center coordinate of a base """

        left_x: float = self._base_left(base)
        right_x: float = self._base_left(base + 1)

        return left_x + ((right_x-left_x)/2)
    
    @property
    def _max_seq_name_width(self) -> float:
        """ Get the width of the longest sequence name """

        max_width: float = 0

        for sequence in self.alignment:
            width: float = stringWidth(f"{sequence.id} (r)", self.seq_name_font, self.seq_name_font_size)
            if width > max_width:
                max_width = width
        
        return max_width
    
    def _font_height(self, font, size) -> float:
        """ Get the height of a font """

        _, bottom, _, top = String(0, 0, string.ascii_letters + string.digits + "_", fontName=font, fontSize=size).getBounds()
        return top - bottom
    
    def _hex_to_color(self, hex: str) -> Color:
        """ Convert a hex color to rgb """

        hex = hex.lstrip("#")
        color_list = [int(hex[i:i+2], 16)/256 for i in (0, 2, 4)]
        return Color(color_list[0], color_list[1], color_list[2])
    
    def _sort_similar(self) -> list[int]:
        """ Sort sequences by similarity to the reference sequence 
        returns list of indexes"""

        return sorted(range(len(self.mutations_list)), key=lambda x: len(self.mutations_list[x]))
    
    def _g_to_a_diamond(self, x: float, y: float) -> Rect:
        """ Draw a rectangle for a G->A mutation """

        diamond = PolyLine([x, y-((self._seq_height/3)/2), x-((self._seq_height/3)/2), y, x, y+((self._seq_height/3)/2), x+((self._seq_height/3)/2), y, x, y-((self._seq_height/3)/2), x-((self._seq_height/3)/2), y], strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=2)

        return diamond
    
    def _get_index_by_id(self, id: str) -> int:
        """ Get the index of a sequence by its id """

        for index, sequence in enumerate(self.alignment):
            if sequence.id == id:
                return index
        
        return None
    
    def _indexes_by_tree_order(self) -> list[int]:
        """ Get the indexes of the sequences in the order they appear in the tree """

        indexes: list[int] = []

        for terminal in self.tree.get_terminals():
            indexes.append(self._get_index_by_id(terminal.name))

        return indexes
    
    @staticmethod
    def significant_digits(value: float, digits: int=2) -> float:
        """ Round a number to a certain number of significant digits """

        if value == 0:
            return 0

        if value > 0:
            value_str = str(value)
            value_len = len(value_str)

            if value_len > digits:
                return(int(value_str[:digits] + "0" * (value_len-digits)))
            else:
                return value

        return value

Graphics.MutationPlot = MutationPlot