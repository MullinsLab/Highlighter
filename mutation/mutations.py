from functools import cache
from typing import Union

import Bio
from Bio import Graphics
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Mutations:
    """ Get mutation info from an alignment """

    def __init__(self, alignment, *, seq_type: str=None, codon_offset: int=0):
        """ Initialize the Mutations object """

        self.alignment = alignment
        self.codon_offet: int=codon_offset % 3

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")
        else:
            self.seq_type = seq_type
        
    def get_seq_index_by_id(self, id: str) -> str:
        """ Get a sequence from the alignment by its id """

        for index, sequence in enumerate(self.alignment):
            if sequence.id == id:
                return index
        
        raise IndexError(f"Could not find sequence with id {id}")

    def list_mismatches(self, *, references: Union[int, str]=0, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, codon_offset: int=0) -> list[dict[str: list]]:
        """ Get matches from a sequence and a reference sequence """

        mismatches: list[dict[str: list]] = []

        if isinstance(references, str):
            references = self.get_seq_index_by_id(references)
        elif not isinstance(references, int):
            raise TypeError(f"Expected reference to be an int or str, got {type(references)}")

        self.references = references

        reference_object = self.alignment[references]

        if not reference_object:
            raise ValueError(f"Reference sequence {references} is empty")
        
        for sequence in self.alignment:
            mismatches.append(self.get_mismatches(sequence=sequence, references=reference_object, seq_type=self.seq_type, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=codon_offset))

        return mismatches

    @staticmethod
    def get_mismatches(*, sequence: Union[str, Seq, SeqRecord], references: Union[str, Seq, SeqRecord], seq_type: str=None, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, codon_offset: int=0) -> dict[int: list]:
        """ Get mismatches from a sequence and a reference sequence 
        returns a dictionary of mismatches where the key is the position of the mismatch and the value is a list of types of mismatches """

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")

        if isinstance(sequence, Seq):
                sequence = str(sequence)
        elif isinstance(sequence, SeqRecord):
                sequence = str(sequence.seq)
        elif not isinstance(sequence, str):
                raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")

        if isinstance(references, Seq):
                references = str(references)
        elif isinstance(references, SeqRecord):
                references = str(references.seq)
        elif not isinstance(references, str):
                raise TypeError(f"Expected reference to be a string, Seq, or SeqRecord, got {type(references)}")
            
        if len(sequence) != len(references):
            raise ValueError("Reference and sequence must be the same length")
        
        return Mutations.get_mismatches_from_str(sequence=sequence, references=references, seq_type=seq_type, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=codon_offset)
    
    @staticmethod
    @cache
    def get_mismatches_from_str(*, sequence: str, references: str, seq_type: str, apobec: bool, g_to_a: bool, stop_codons: bool=False, glycosylation: bool, codon_offset: int=0) -> dict[int: list]:
        """ Get mutations from a sequence and a reference sequence
        separated out so it can be cached (Seq and SeqRecord are not hashable) """

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")

        mismatches: dict = {}

        if sequence == references and (seq_type == "NT" and not stop_codons) and (seq_type == "AA" and not glycosylation):
            return mismatches

        for base_index in range(len(sequence)):
            if references[base_index] != sequence[base_index]:
                mismatches[base_index] = []
                
                if sequence[base_index] != "-":
                    mismatches[base_index].append(sequence[base_index])
                else:
                    mismatches[base_index].append("Gap")

                # APOBEC and G->A mutations only apply to NT sequences
                if seq_type == "NT":
                    if references[base_index] == "G" and sequence[base_index] == "A":
                        if g_to_a:
                            mismatches[base_index].append("G->A mutation")

                        if apobec and base_index <= len(sequence)-3 and sequence[base_index+1] in "AG" and sequence[base_index+2] != "C":
                            mismatches[base_index].append("APOBEC")

            # Stop codons only apply to NT sequences
            if seq_type == "NT":
                if stop_codons and sequence[base_index] in "TU" and base_index <= len(sequence)-3 and SeqUtils.codon_position(sequence, base_index, codon_offset=codon_offset) == 0:
                    base_snippet: str = ""
                    snippet_index: int = base_index+1

                    while len(base_snippet) < 2 and snippet_index < len(sequence):
                        base_snippet += sequence[snippet_index] if sequence[snippet_index] != "-" else ""
                        snippet_index += 1

                        if base_snippet in ("AA", "AG", "GA"):
                            if base_index not in mismatches:
                                mismatches[base_index] = []

                            mismatches[base_index].append("Stop codon")

            # Glycosylation only applies to AA sequences
            elif seq_type == "AA":
                if glycosylation and sequence[base_index] == "N" and base_index <= len(sequence)-3:
                    base_snippet: str = ""
                    snippet_index: int = base_index+1

                    while len(base_snippet) < 2 and snippet_index < len(sequence):
                        base_snippet += sequence[snippet_index] if sequence[snippet_index] != "-" else ""
                        snippet_index += 1
                    
                    if base_snippet[0] != "P" and base_snippet[1] in "ST":
                        if base_index not in mismatches:
                            mismatches[base_index] = []

                        mismatches[base_index].append("Glycosylation")
        
        return mismatches
    
    def list_matches(self, *, references: tuple[Union[int, str]]=0) -> list[dict[str: list]]:
        """ Get matches from a sequence and a reference sequence """

        matches: list[dict[str: list]] = []

        if not isinstance(references, list):
            references = [references]

        reference_objects: list = []
        self.references = []

        for reference in references:
            if isinstance(reference, str):
                self.references.append(self.get_seq_index_by_id(reference))

            elif isinstance(reference, int):
                if reference >= len(self.alignment):
                    raise IndexError(f"Reference index {reference} is out of range")
                
                self.references.append(reference)

            else:
                raise TypeError(f"Expected reference to be an int or str, got {type(reference)}")
            
            reference_objects.append(self.alignment[self.references[-1]])

        for sequence in self.alignment:
            matches.append(self.get_matches(sequence=sequence, references=reference_objects, seq_type=self.seq_type))

        return matches
    
    @staticmethod
    def get_matches(*, sequence: Union[str, Seq, SeqRecord], references: Union[list[str, Seq, SeqRecord], str, Seq, SeqRecord], seq_type: str=None) -> dict[int: list]:
        """ Get matches of a sequence to one or more reference sequences
        returns a dictionary of matches where the key is the position of the match and the value is a list of types of matches """
        
        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")
        
        if isinstance(sequence, Seq):
                sequence = str(sequence)
        elif isinstance(sequence, SeqRecord):
                sequence = str(sequence.seq)
        elif not isinstance(sequence, str):
                raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")
        
        if not isinstance(references, list):
            references = [references]

        new_references: list = []

        for reference in references:
            if isinstance(reference, Seq):
                reference = str(reference)
            elif isinstance(reference, SeqRecord):
                reference = str(reference.seq)
            elif not isinstance(reference, str):
                raise TypeError(f"Expected reference to be a string, Seq, or SeqRecord, got {type(reference)}")
            
            if len(sequence) != len(reference):
                raise ValueError("All references and sequence must be the same length")
            
            new_references.append(reference)

        return Mutations.get_matches_from_str(sequence=sequence, references=tuple(new_references), seq_type=seq_type)

    @staticmethod
    @cache
    def get_matches_from_str(*, sequence: str, references: tuple[str], seq_type: str) -> dict[int: list]:
        """ Get matches from a sequence and a reference sequence
        separated out so it can be cached (Seq and SeqRecord are not hashable) """

        if seq_type not in ("NT", "AA"):
            raise ValueError("type must be provided (either 'NT' or 'AA')")

        matches: dict = {}
        
        if sequence in references:
            return matches

        for base_index in range(len(sequence)):
            matches[base_index] = []
            for reference_index, reference in enumerate(references):
                if reference[base_index] == sequence[base_index]:
                    matches[base_index].append(reference_index)

            if not matches[base_index]:
                matches[base_index].append("Unique")

            elif len(matches[base_index]) == len(references):
                del matches[base_index]

        return matches
    
    def export_matches(self, output_file, *, references: tuple[Union[int, str]]=0) -> None:
        """ Export matches to a .txt file """

        output: str = ""
        matches = self.list_matches(references=references)
        
        for sequence_index, sequence in enumerate(self.alignment):
            matched: dict = {}

            for base, codes in matches[sequence_index].items():
                if "Unique" in codes:
                    if "Unique" not in matched:
                        matched["Unique"] = []
                    
                    matched["Unique"].append(base+1)
                
                elif len(codes) > 1:
                    if "Multiple" not in matched:
                        matched["Multiple"] = []
                    
                    matched["Multiple"].append(base+1)

                else:
                    if codes[0] not in matched:
                        matched[codes[0]] = []
                    
                    matched[codes[0]].append(base+1)

            if sequence.id == "696": 
                print(matched)
                print(matches[sequence_index])
                print(references)

            if sequence_index in self.references:
                output += f"{sequence.id} (R{self.references.index(sequence_index)+1})\n"
            else:
                output += f"{sequence.id}\n"

            for code in list(range(10)) + ["Unique", "Multiple"]:
                if code not in matched:
                    continue

                if code == "Unique":
                    output += f"Unique in query "

                elif code == "Multiple":
                    output += f"Multiple matches "
                    
                else:
                    output += f"R{code+1} "

                output += f"[{' '.join([str(thing) for thing in matched[code]])}]\n"

            output += "\n"

        with open(output_file, mode="wt") as file:
            file.write(output)
    
AlignInfo.Mutations = Mutations


from reportlab.lib import colors
from reportlab.lib.colors import Color
from reportlab.lib.units import inch

from reportlab.pdfbase.pdfmetrics import stringWidth

from reportlab.graphics.shapes import Drawing, String, Line, Rect, Circle, PolyLine, Polygon

from Bio import SeqUtils
from Bio.Graphics import _write
from Bio.Align import AlignInfo


class MutationPlot:
    """ Create and output a mutation plot """

    def __init__(self, alignment, *, type: str=None, tree: Union[str, object]=None, plot_width: int = 4*inch, seq_name_font: str="Helvetica", seq_name_font_size: int=8, seq_gap: int=None, left_margin: float=.25*inch, top_margin: float=.25*inch, bottom_margin: float=0, right_margin: float=0, mark_reference: bool=True, title_font="Helvetica", title_font_size: int=12, ruler: bool=True, ruler_font: str="Helvetica", ruler_font_size: int=6, ruler_major_ticks: int=10, ruler_minor_ticks=3, codon_offset: int=0):
        """ Initialize the MutationPlot object """

        self.alignment = alignment
        self.codon_offset = codon_offset % 3

        if type not in ("NT", "AA"):
            self.type = self.guess_alignment_type(alignment)
        else:
            self.type = type

        if tree is not None:
            if isinstance(tree, Bio.Phylo.BaseTree.Tree):
                self.tree = tree
            else:
                raise TypeError("Tree must be a Bio.Phylo.BaseTree.Tree object (or a derivative)")

        self._seq_count = len(alignment)
        self._seq_length = len(alignment[0])
        self.mark_reference: bool = mark_reference

        self.seq_name_font: str = seq_name_font
        self.seq_name_font_size: int = seq_name_font_size

        self.top_margin: float = top_margin
        self.bottom_margin: float = bottom_margin
        self.left_margin: float = left_margin
        self.right_margin: float = right_margin

        self.title_font: str = title_font
        self.title_font_size: int = title_font_size

        self.ruler: bool = ruler
        self.ruler_font: str = ruler_font
        self.ruler_font_size: int = ruler_font_size
        self.ruler_major_ticks: int = ruler_major_ticks
        self.ruler_minor_ticks: int = ruler_minor_ticks
        self._ruler_font_height: float = self.ruler_font_size

        self.plot_width: float = plot_width

        self._seq_height: float = self.seq_name_font_size
        self.seq_gap: float = self._seq_height / 5 if seq_gap is None else seq_gap

        self.match_plot_colors: dict[str: dict] = {
            "ML": {
                "references": ["#FF0000", "#537EFF", "#00CB85", "#000000", "#FFA500"],
                "unique": "#EFE645",
                "multiple": "#808080",
            },
            "LANL": {
                "references": ["#ED1C24", "#235192", "#FFC20E", "#00A651", "#8DC73F"],
                "unique": "#000000",
                "multiple": "#666666",
            }
        }
        self.match_plot_unique_color: str = "#EFE645"
        self.match_plot_multiple_color: str = "#808080"

        self.mismatch_plot_colors: dict[str: [dict[str: str]]] = {
            "NT": {
                "LANL": {
                    "A": "#42FF00", 
                    "C": "#41B8EE", 
                    "G": "#FFA500", 
                    "T": "#EE0B10", 
                    "Gap": "#666666",
                },
                "ML": {
                    "A": "#36b809", 
                    "C": "#1282b5", 
                    "G": "#FFA500", 
                    "T": "#c48002", 
                    "Gap": "#666666",
                }
            },
            "AA": {
                "LANL": {
                    "H": "#FF0000",
                    "D": "#302ECD",
                    "E": "#302ECD",
                    "K": "#23659B",
                    "N": "#23659B",
                    "Q": "#23659B",
                    "R": "#23659B",
                    "M": "#2F9A2F",
                    "I": "#42FF00",
                    "L": "#42FF00",
                    "V": "#42FF00",
                    "F": "#F900FF",
                    "W": "#F900FF",
                    "Y": "#F900FF",
                    "C": "#CD2F2E",
                    "A": "#F9CE2E", 
                    "G": "#F9CE2E",
                    "S": "#F9CE2E",
                    "T": "#F9CE2E",
                    "P": "#FBFF00",
                    "Other": "#000000",
                    "Gap": "#bebebe",
                }
            }
        }

    def setup_drawing(self, *, plot_type: str, output_format: str="svg", title: str=None, sort: str="similar", mark_width: float=1, scale: float=1, scheme: str="LANL"):
        """ Setus up the drawing """
        
        self.output_format: str = output_format

        self.title: str = title
        self._title_font_height: float = self.title_font_size
        self._title_height = 0 if not title else self._title_font_height*2

        self._ruler_height = 0 if not self.ruler else self._ruler_font_height * 3

        self._plot_floor: float = self.bottom_margin + self._ruler_height

        self._seq_name_width: float = self._max_seq_name_width
        self._width: float = self.left_margin + self.plot_width + (inch/4) + self._seq_name_width + self.right_margin
        
        self._height: float = len(self.alignment) * (self._seq_height + self.seq_gap) + self.top_margin + self.bottom_margin + self._title_height + self._ruler_height

        self.mark_width: float = mark_width

        self.scale = scale

        self.drawing = Drawing(self._width, self._height)

        if sort == "similar":
            self.sorted_keys = self._sort_similar()
        
        elif sort == "tree":
            if self.tree is None:
                raise ValueError("Cannot sort by tree if no tree is provided")
            
            self.sorted_keys = self._indexes_by_tree_order()

        else: 
            self.sorted_keys = range(len(self.matches_list))

        self._draw_title()
        if self.ruler:
            self._draw_ruler()

        for plot_index, seq_index in enumerate(self.sorted_keys):

            # Add label for sequence
            id = self.alignment[seq_index].id
            if self.mark_reference:
                if isinstance(self._mutations.references, int):
                    if seq_index == self._mutations.references:
                        id += " (r)"
                elif seq_index in self._mutations.references:
                    id += f" (r{self._mutations.references.index(seq_index)+1})"
                

            x: float = self.left_margin + self.plot_width + (inch/4)
            y: float = ((self._seq_count-(plot_index + .75)) * (self._seq_height + self.seq_gap))  + self.seq_gap + self._plot_floor
            sequence_str: String = String(x, y, id, fontName="Helvetica", fontSize=self.seq_name_font_size)
            self.drawing.add(sequence_str, id)

            # Add base line for sequence
            if plot_type == "match" and seq_index in self._mutations.references:
                color: Color = self._hex_to_color(self._current_scheme[self._mutations.references.index(seq_index)])
            else:
                color: Color = colors.lightgrey

            x1: float = self.left_margin
            x2: float = self.left_margin + self.plot_width
            y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self.seq_gap) + self.seq_gap + self._plot_floor
            sequence_baseline: Line = Line(x1, y, x2, y, strokeColor=color)
            self.drawing.add(sequence_baseline)
    
    def draw_mismatches(self, output_file, *, output_format: str="svg", title: str=None, reference: Union[str, int]=0, apobec: bool=False, g_to_a: bool=False, stop_codons: bool=False, glycosylation: bool=False, sort: str="similar", mark_width: float=1, scheme: str="LANL", scale: float=1):
        """ Draw mismatches compared to a reference sequence """

        self._mutations = AlignInfo.Mutations(self.alignment, seq_type=self.type)
        
        self.matches_list = self._mutations.list_mismatches(references=reference, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons, glycosylation=glycosylation, codon_offset=self.codon_offset)
        self.references = self._mutations.references

        self.setup_drawing(output_format=output_format, title=title, sort=sort, mark_width=mark_width, scale=scale, plot_type="mismatch")

        self.scheme: str = scheme
        self._current_scheme: dict = self.mismatch_plot_colors[self.type][self.scheme] if self.scheme in self.mismatch_plot_colors[self.type] else self.mismatch_plot_colors[self.type]["LANL"]

        self._apobec: bool = apobec
        self._g_to_a: bool = g_to_a
        self._glycosylation: bool = glycosylation
        self._stop_codons: bool = stop_codons

        for plot_index, seq_index in enumerate(self.sorted_keys):
            matches = self.matches_list[seq_index]
            self._draw_marks_mismatch(plot_index, matches, is_reference=(seq_index == self.references))

        return _write(self.drawing, output_file, self.output_format, dpi=288*self.scale)

    def _draw_marks_mismatch(self, plot_index: int, matches: dict[int: list], is_reference: bool=False) -> None:
        """ Draw marks for a mismatch sequence """

        for base, match_item in matches.items():
            for code in match_item:
                if code in self._current_scheme:
                    color: Color = self._hex_to_color(self._current_scheme[code])
                    self.drawing.add(self._base_mark(plot_index, base, color))
        
        # Symboloic markers need to be drawn second so they are on top of the rectangles
        for base, match_item in matches.items():
            x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
            y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self.seq_gap) + self.seq_gap + self._plot_floor
                
            if "APOBEC" in match_item:
                self.draw_circle(x, y)
                
            elif "G->A mutation" in match_item:
                self.draw_diamond(x, y)

            elif "Glycosylation" in match_item:
                if is_reference:
                    self.draw_circle(x, y)
                else:
                    if "Glycosylation" not in self.matches_list[self.references].get(base, {}):
                        self.draw_diamond(x, y, filled=True)

            elif "Stop codon" in match_item:
                self.draw_diamond(x, y, color="#0000FF")

        if self.type == "AA" and self._glycosylation:
            for base, match_item in self.matches_list[self.references].items():
                if "Glycosylation" in match_item and "Glycosylation" not in matches.get(base, {}):
                    x: float = self.left_margin + self._base_left(base) + ((self._base_left(base+1)-self._base_left(base))/2)
                    y: float = (self._seq_count-(plot_index + .5)) * (self._seq_height + self.seq_gap) + self.seq_gap + self._plot_floor

                    self.draw_diamond(x, y, color="#0000FF")

    def draw_matches(self, output_file, *, output_format: str="svg", title: str=None, references: list[Union[str, int]]=0, sort: str="similar", mark_width: float=1, scheme: str="LANL", scale: float=1):
        """ Draw mismatches compared to a reference sequence """

        self._mutations = AlignInfo.Mutations(self.alignment, seq_type=self.type)
        
        self.matches_list = self._mutations.list_matches(references=references)
        self.references = self._mutations.references

        self.scheme: str = scheme
        self._current_scheme: list = self.match_plot_colors[self.scheme]["references"]
        self._current_unique_color: str = self.match_plot_colors[self.scheme]["unique"]
        self._current_multiple_color: str = self.match_plot_colors[self.scheme]["multiple"]

        self.setup_drawing(output_format=output_format, title=title, sort=sort, mark_width=mark_width, scale=scale, plot_type="match", scheme=scheme)

        for plot_index, seq_index in enumerate(self.sorted_keys):
            matches = self.matches_list[seq_index]
            self._draw_marks_match(plot_index, matches, is_reference=(seq_index in self.references))

        return _write(self.drawing, output_file, self.output_format, dpi=288*self.scale)

    def _draw_marks_match(self, plot_index: int, matches: dict[int, list], is_reference: bool) -> None:
        """ Draw the marks for a match sequence """

        for base, match_item in matches.items():
            if "Unique" in match_item:
                color: Color = self._hex_to_color(self._current_unique_color)
                self.drawing.add(self._base_mark(plot_index, base, color))

            elif len(match_item) > 1:
                color: Color = self._hex_to_color(self._current_multiple_color)
                self.drawing.add(self._base_mark(plot_index, base, color))

            else:
                for code in match_item:
                    color: Color = self._hex_to_color(self._current_scheme[code])
                    self.drawing.add(self._base_mark(plot_index, base, color))

    def _base_mark(self, plot_index, base, color):
        """ Returns a mark for a particular base """
        x1: float = self.left_margin + self._base_left(base)
        x2: float = self.left_margin + self._base_left(base+self.mark_width)

        y1: float = ((self._seq_count-plot_index) * (self._seq_height + self.seq_gap)) + (self.seq_gap/2) + self._plot_floor
        y2: float = ((self._seq_count-(plot_index+1)) * (self._seq_height + self.seq_gap)) + self.seq_gap + self._plot_floor

        return Rect(x1, y1, x2-x1, y2-y1, fillColor=color, strokeColor=color, strokeWidth=0.1)

    def _draw_title(self) -> None:
        """ Draw the title at the top of the plot """
            
        if self.title:
            x: float = self.left_margin + (self.plot_width/2)
            y: float = self._height - self.top_margin - (self._title_font_height/2)

            self.drawing.add(String(x, y, self.title, textAnchor="middle", fontName=self.title_font, fontSize=self.title_font_size))

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
        x2: float = self.left_margin + self.plot_width
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

        # if base == self._seq_length:
        #     return self.plot_width
        
        return (base / self._seq_length) * self.plot_width
    
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
            width: float = stringWidth(f"{sequence.id} (r10)", self.seq_name_font, self.seq_name_font_size)
            if width > max_width:
                max_width = width
        
        return max_width
    
    def _hex_to_color(self, hex: str) -> Color:
        """ Convert a hex color to rgb """

        hex = hex.lstrip("#")
        color_list = [int(hex[i:i+2], 16)/256 for i in (0, 2, 4)]
        return Color(color_list[0], color_list[1], color_list[2])
    
    def _sort_similar(self) -> list[int]:
        """ Sort sequences by similarity to the reference sequence 
        returns list of indexes"""

        return sorted(range(len(self.matches_list)), key=lambda x: len(self.matches_list[x]))
    
    def draw_diamond(self, x: float, y: float, color: str="#FF00FF", filled: bool=False) -> None:
        """ Draw a rectangle on the plot """
        
        fill_color = self._hex_to_color(color) if filled else None

        diamond = Polygon([x, y-((self._seq_height/3)/2), x-((self._seq_height/3)/2), y, x, y+((self._seq_height/3)/2), x+((self._seq_height/3)/2), y], strokeColor=self._hex_to_color(color), strokeWidth=2, fillColor=fill_color)
        
        self.drawing.add(diamond)

    def draw_circle(self, x: float, y: float, color: str="#FF00FF", filled: bool=True) -> None:
        """ Draw a circle on the plot"""
        
        circle = Circle(x, y, (self._seq_height/3)/2, fillColor=self._hex_to_color(color), strokeColor=self._hex_to_color("#FF00FF"), strokeWidth=0.1)
        self.drawing.add(circle)
    
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
    
    @staticmethod
    def guess_alignment_type(alignment) -> str:
        """ Returns either 'NT' (nucleotide) or 'AA' (amino acid) """

        nt_codes: str = "ACGTUiRYKMSWBDHVN-" # IUPAC nucleotide codes

        for sequence in alignment:
            if isinstance(sequence, SeqRecord):
                sequence_str = str(sequence.seq)
            elif isinstance(sequence, Seq):
                sequence_str = str(sequence)
            elif isinstance(sequence, str):
                sequence_str = sequence
                
            for symbol in sequence_str.strip():
                if symbol not in nt_codes:
                    return "AA"
                    
        return "NT"

Graphics.MutationPlot = MutationPlot


from Bio import SeqUtils

@cache
def codon_position(sequence: Union[str, Seq, SeqRecord], base: int, *, codon_offset: int=0) -> int:
    """ Get the codon position of a base in a sequence
    returns 0 for the first base of a codon, 1 for the second, or 2 for the third) """

    if isinstance(sequence, Seq):
            sequence = str(sequence)
    elif isinstance(sequence, SeqRecord):
            sequence = str(sequence.seq)
    elif not isinstance(sequence, str):
            raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")
    
    if not isinstance(base, int):
        raise TypeError(f"Expected position to be an int, got {type(base)}")

    if base > len(sequence):
        raise ValueError(f"Position {base} is greater than the length of the sequence ({len(sequence)})")
    
    if sequence[base] == "-":
        raise ValueError(f"Position {base} is a gap")
    
    adjusted_base: int =  base-sequence[:base+1].count("-")
    return ((adjusted_base + codon_offset) % 3)

SeqUtils.codon_position = codon_position