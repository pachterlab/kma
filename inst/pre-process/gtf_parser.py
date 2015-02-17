# A Python GTF Parser
# Copyright (C) 2015 Riyaz Faizullabhoy
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# Riyaz Faizullabhoy
# 7/9/2012
# GTF --> GTF Objects

import argparse
import bisect
import copy
import sys
import warnings

from bisect import *


def run():
    """Runs gtf_parser.py and creates transcript objects for the input GTF
    file.

    Usage: python gtf_parser.py <gtf-input-file-name>
    The GTF file must be in the current directory

    """
    parser = argparse.ArgumentParser()

    parser.add_argument('gi',
        metavar='gtf_input',
        help="The GTF annotation input file name")

    try:
        args = parser.parse_args()
    except ValueError as IOError:
        print >> sys.stderr, "Parsing Error, please use the following ' \
                'command-line format:\npython gtf_parser.py <GTF_FILE_NAME>"

    return gtf_parse(args.gi)


class Transcript(object):

    """An object for a given transcript id in the GTF file.

    Will be in the output.

    """

    def __init__(
            self,
            transcript_id,
            refname,
            strand,
            frame,
            gene_id_attributes,
            gene_id,
            score = None,
            source = None):
        self.transcript_id = transcript_id
        self.refname = refname
        self.strand = strand
        if strand == "+":
            self._is_reverse = False
        else:
            self._is_reverse = True

        self.frame = frame
        self.gene_id_attributes = gene_id_attributes
        self.gene_id = gene_id
        self.exons = []
        self._furthest_added_exon = None

        if score == '.':
            self.score = None
        else:
            self.score = score

        self.source = source

    def add_exon(self, exon):
        """Adds an exon to the transcript object.

        Throws an exception if the exon is out of order in the GTF file,
        or when it overlaps another exon.  If an exception is thrown for
        a given transcript, it will not be included for the output.

        """

        if exon[0] > exon[1]:
            raise Exception(
                "Invalid exon start/stop in transcript: " + \
                        str(self.transcript_id))

        if exon[0] == self._furthest_added_exon:
            raise Exception( 'Non-sensical exons. ' + \
                    'One begins right after the other ends.' +  \
                    ' This should all be one exon: {0}'.format( self ) )

        # Add an exon to the end (i.e. exons are in order)
        if self._furthest_added_exon < exon[0] or \
                self._furthest_added_exon is None:
            self.exons += [exon]
            self._furthest_added_exon = exon[1]
            return

        for e in self.exons:
            if (e[1] >= exon[0]) and (exon[1] >= e[0]):
                raise Exception(
                    "Overlapping exon in transcript: " + str(self.transcript_id))
        # Add an exon elsewhere, to keep exons sorted
        index = bisect_left( self.exons, exon )
        self.exons.insert(index, exon)

    @property
    def front_coordinate(self):
        """Returns the left-most coordinate of the left-most exon of the
        transcript object."""

        return self.exons[0][0]

    @property
    def end_coordinate(self):
        """Returns the right-most coordinate of the last (right-most) exon of
        the transcript object."""

        return self.exons[-1][1]

    def compatible(self, pos):
        """Returns the exon index if 'pos' is inside one of the exons.

        Else, returns None

        """
        for idx, exon in enumerate(self.exons):
            if pos >= exon[0] and pos < exon[1]:
                return idx
        return None

    def to_gtf(self):
        if len(self.exons) == 0:
            raise Exception( 'Cannot print to GTF if no exons' )
        trans_str = []
        trans_str.append(self.refname)
        trans_str.append(self.source if self.source is not None else 'NA')
        trans_str.append('exon')

        all_exons = []
        for ex in self.exons:
            exon_str = copy.deepcopy(trans_str)
            exon_str.append(str(ex[0] + 1))
            exon_str.append(str(ex[1]))
            if self.score is None:
                exon_str.append('.')
            else:
                exon_str.append(self.score)
            exon_str.append('-' if self.is_reverse else '+')
            exon_str.append(str(self.frame))
            exon_str.append('gene_id "' + self.gene_id + '";' +
                            ' transcript_id "' + self.transcript_id + '";')
            # FIXME: after gene_id_attributes is parsers correctly, add this
            # line.  (and a space)
            # self.gene_id_attributes)
            all_exons.append('\t'.join(exon_str))

        return '\n'.join(all_exons)

    @property
    def is_reverse(self):
        return self._is_reverse

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '{0}:{1}:{2}-{3}'.format(
            self.transcript_id,
            self.refname,
            self.front_coordinate,
            self.end_coordinate)


def gtf_parse(input_gtf):
    """Parses the input GTF file by line.

    Creates a list of transcript objects that include each transcript's id,
    refname, strand, frame, gene id, and its exons.  Uses pythonic 0-based
    right-exclusive coordinates.  For example, start: 1, end: 2 would be an
    exon of length one (at 2) in genomic one-based coordinates.

    """
    transcript_dictionary = {}
    bad_transcripts = []
    current_transcript = None
    gtf_file = open(input_gtf, 'r')
    for line in gtf_file:

        if line.startswith('#'):
            continue

        gtf_line = (line.split("\t"))

        if len(gtf_line) < 9 or gtf_line[2] != "exon":
            continue

        try:
            transcript_id = (
                gtf_line[8].split(";")[1].split(" ")[2].replace(
                    "\"",
                    ""))
            if transcript_id in bad_transcripts:
                continue

        except IndexError as i:
            print >> sys.stderr, "GTF File Input missing 'transcript_id' field"

        if current_transcript is None:
            try:
                # FIXME: gene_id_attributes isn't being parsed correctly
                gene_id = gtf_line[8].split(";")[0].split(
                    " ")[1].replace("\"", "")
                current_transcript = Transcript(
                    transcript_id,
                    gtf_line[0],
                    gtf_line[6],
                    gtf_line[7],
                    gtf_line[8].split(" ")[3],
                    gene_id,
                    gtf_line[5],
                    gtf_line[1])
                transcript_dictionary[transcript_id] = current_transcript
                transcript_dictionary[transcript_id].add_exon(
                    ( int(gtf_line[3]) - 1, int(gtf_line[4]) ) )
            except IndexError as i:
                print >> sys.stderr, "GTF File Input missing fields"
            except Exception as e:
                print >> sys.stderr, e
                bad_transcripts += current_transcript.transcript_id
                current_transcript = None  # throw out the bad transcript
                transcript_dictionary[transcript_id] = None

        elif transcript_id != current_transcript.transcript_id:
            try:
                transcript_dictionary[transcript_id].add_exon(
                    ( int(gtf_line[3]) - 1,
                        int(gtf_line[4]) ) )
                current_transcript = transcript_dictionary[transcript_id]
            except KeyError as k:
                gene_id = gtf_line[8].split(";")[0].split(
                    " ")[1].replace("\"", "")
                # FIXME: gene_id_attributes isn't being parsed correctly
                current_transcript = Transcript(
                    transcript_id,
                    gtf_line[0],
                    gtf_line[6],
                    gtf_line[7],
                    gtf_line[8].split(" ")[
                        4:],
                    gene_id,
                    gtf_line[5],
                    gtf_line[1])
                transcript_dictionary[transcript_id] = current_transcript
                transcript_dictionary[transcript_id].add_exon(
                    (int(gtf_line[3]) - 1,
                        int(gtf_line[4])))

            except IndexError as i:
                print >> sys.stderr, "GTF File Input missing fields"
            except Exception as e:
                print >> sys.stderr, (e)
                bad_transcripts += current_transcript.transcript_id
                current_transcript = None
                transcript_dictionary[transcript_id] = None

        else:
            try:
                current_transcript.add_exon(
                    (int(gtf_line[3]) - 1,
                        int(gtf_line[4])))
            except IndexError as i:
                print >> sys.stderr, "GTF File Input missing fields"
            except Exception as e:
                print >> sys.stderr, e
                bad_transcripts += current_transcript.transcript_id
                current_transcript = None
                transcript_dictionary[transcript_id] = None

    gtf_file.close()

    return transcript_dictionary


def gtf_write(all_trans, out_handle):
    all_trans = sorted( all_trans.values(), key = lambda t: (t.refname,
        t.front_coordinate) )
    print all_trans
    for trans in all_trans:
        print >> out_handle, trans.to_gtf()

if (__name__ == "__main__"):
    run()
