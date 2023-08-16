import string

from functools import cache

from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Graphics

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

    def __init__(self, alignment, *, output_format="svg", seq_name_font="Helvetica", seq_name_size=20, left_margin: float=.25, top_margin: float=.25, botom_margin: float=0, mark_reference: bool=True):
        """ Initialize the MutationPlot object """

        self.alignment = alignment
        self.mutations = AlignInfo.Mutations(alignment)
        self.seq_count = len(alignment)
        self.seq_length = len(alignment[0])
        self.mark_reference: bool = mark_reference

        self.output_format: str = output_format
        self.seq_name_font: str = seq_name_font
        self.seq_name_size: int = seq_name_size

        self.top_margin: float = inch * top_margin
        self.bottom_margin: float = inch * botom_margin
        self.left_margin: float = inch * left_margin

        self.plot_width: float = 14 * inch
        self.seq_name_width: float = self._max_seq_name_width
        self.width: float = self.left_margin + self.plot_width + (inch/4) + self.seq_name_width
        
        self.seq_height: float = self._font_height
        self.seq_gap: float = self.seq_height / 5
        self.height: float = len(self.alignment) * (self.seq_height + self.seq_gap) + self.top_margin + self.bottom_margin

        self.plot_colors: dict[str: str] = {"A": "#42FF00", "C": "#41B8EE", "G": "#FFA500", "T": "#EE0B10", "Gap": "#666666"}

        # self.apobec_circle: Circle = Circle(0, 0, (self.seq_height/3)/2, fillColor=self._hex_to_color("#FF00FF"), strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=0.1)

    def draw(self, output_file, reference: str|int=0, apobec: bool=False, g_to_a: bool=False, sort: str="similar", narrow_markers: bool=True, min_marker_width: float=1):
        """ Writes out the mutation plot to a file """
        
        drawing = self.drawing = Drawing(self.width, self.height)
        self.mutations_list = self.mutations.list_mutations(reference=reference, apobec=apobec, g_to_a=g_to_a)
        self.narrow_markers: bool = narrow_markers
        self.min_marker_width: float = min_marker_width

        if sort == "similar":
            sorted_keys = self._sort_similar()
        else: 
            sorted_keys = range(len(self.mutations_list))

        # for seq_index, mutations in enumerate(self.mutations_list):
        for plot_index, seq_index in enumerate(sorted_keys):
            mutations = self.mutations_list[seq_index]

            # Add label for sequence
            id = self.alignment[seq_index].id
            if self.mark_reference and seq_index == self.mutations.reference:
                id += " (r)"

            x: float = self.left_margin + self.plot_width + (inch/4)
            y: float = ((self.seq_count-(plot_index + .5)) * (self.seq_height + self.seq_gap)) + self.bottom_margin
            sequence_str: String = String(x, y, id, fontName="Helvetica", fontSize=self.seq_name_size)
            drawing.add(sequence_str)

            # Add base line for sequence
            x1: float = self.left_margin
            x2: float = self.left_margin + self.plot_width
            y: float = (self.seq_count-(plot_index + .5)) * (self.seq_height + self.seq_gap) + self.seq_gap + self.bottom_margin
            sequence_baseline: Line = Line(x1, y, x2, y, strokeColor=colors.lightgrey)
            drawing.add(sequence_baseline)

            self._draw_mutations(plot_index, mutations)

        return _write(drawing, output_file, self.output_format)
    
    def _draw_mutations(self, plot_index: int, mutations: dict[int: list]) -> None:
        """ Draw mutations for a sequence """

        for base, mutation in mutations.items():
            for code in mutation:
                if code in self.plot_colors:
                    x1: float = self.left_margin + self._base_left(base)
                    x2: float = self.left_margin + self._base_left(base+1)

                    y1: float = ((self.seq_count-plot_index) * (self.seq_height + self.seq_gap)) + (self.seq_gap/2) + self.bottom_margin
                    y2: float = ((self.seq_count-(plot_index+1)) * (self.seq_height + self.seq_gap)) + self.seq_gap + self.bottom_margin

                    if code != "Gap" and self.narrow_markers and x2-x1 > self.min_marker_width:
                        x1, x2 = (
                            x1+(((x2-x1)-self.min_marker_width)/2), 
                            x2-(((x2-x1)-self.min_marker_width)/2)
                        )

                    base_color: Color = self._hex_to_color(self.plot_colors[code])
                    base_mark: Rect = Rect(x1, y1, x2-x1, y2-y1, fillColor=base_color, strokeColor=base_color, strokeWidth=0.1)
                    self.drawing.add(base_mark)
        
        # APOBEC and G->A go second so they go on top of other elements
        for base, mutation in mutations.items():
            if "APOBEC" in mutation:
                print(f"Base: {base}, apobec")
                x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
                y: float = (self.seq_count-(plot_index + .5)) * (self.seq_height + self.seq_gap) + self.seq_gap + self.bottom_margin
                
                circle = Circle(x, y, (self.seq_height/3)/2, fillColor=self._hex_to_color("#FF00FF"), strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=0.1)
                
                self.drawing.add(circle)

            elif "G->A mutation" in mutation:
                print(f"Base: {base}, G->A")
                x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
                y: float = (self.seq_count-(plot_index + .5)) * (self.seq_height + self.seq_gap) + self.seq_gap + self.bottom_margin
                
                diamond = self._g_to_a_diamond(x, y)
                
                self.drawing.add(diamond)
    
    def _base_left(self, base: int) -> float:
        """ Get the left coordinate of a base """

        if base == self.seq_length:
            return self.plot_width
        
        return (base / self.seq_length) * self.plot_width
    
    @property
    def _max_seq_name_width(self) -> float:
        """ Get the width of the longest sequence name """

        max_width: float = 0

        for sequence in self.alignment:
            width: float = stringWidth(f"{sequence.id} (r)", self.seq_name_font, self.seq_name_size)
            if width > max_width:
                max_width = width
        
        return max_width
    
    @property
    def _font_height(self) -> float:
        """ Get the height of a font """

        _, left, _, right = String(0, 0, string.ascii_letters + string.digits + "_", fontName=self.seq_name_font, fontSize=self.seq_name_size).getBounds()
        return right - left
    
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

        #diamond = Rect(x-((self.seq_height/3)/2), y-((self.seq_height/3)/2), self.seq_height/3, self.seq_height/3, strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=1)
        diamond = PolyLine([x, y-((self.seq_height/3)/2), x-((self.seq_height/3)/2), y, x, y+((self.seq_height/3)/2), x+((self.seq_height/3)/2), y, x, y-((self.seq_height/3)/2), x-((self.seq_height/3)/2), y], strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=2)

        return diamond

Graphics.MutationPlot = MutationPlot