# intron ops: tools for intron operations
# Copyright (C) 2015 Harold Pimentel
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import os
import sys

import pysam

from gtf_parser import Transcript

def print_trans_names(trans):
    if len(trans) == 0:
        return

    print trans[0].gene_id
    for t in trans:
        print '\t', t.transcript_id

def reduce_to_gene(trans_list, process_by_gene = None):
    """ Performs 'process_by_gene' on a list of transcripts which all have the
    same gene names. returns the gene -> transcript mapping and returns the gene
    -> result mapping in a tuple. *IMPORTANT* Assumes that transcripts are
    sorted by position"""

    cur_gene = None

    ret_dict = {}
    gene_to_trans = {}

    # makes no assumption about overlaps
    for trans in trans_list:
        cur_gene = trans.gene_id
        gene_list = gene_to_trans.get(cur_gene, [])
        # XXX: need to decide what best design is here
        # if len(gene_list) > 0:
            # if gene_list[0].refname != trans.refname:
            #     print >> sys.stderr, 'ERROR: Gene showed up several times, discarding mismatch'
            #     continue
        gene_list.append(trans)
        gene_to_trans[cur_gene] = gene_list

    if process_by_gene is not None:
        for gene in gene_to_trans:
            ret_dict[gene] = process_by_gene( gene_to_trans[gene] )
        return (gene_to_trans, ret_dict)

    return gene_to_trans

class Intron:
    def __init__(self, ref, start, stop, extend_start = 0, extend_stop = 0):
        self.refname = ref
        self.coords = (start, stop)
        self.extension = (extend_start, extend_stop)

    def __getitem__(self, i):
        if i == 0:
            return self.coords[i] - self.extension[i]
        elif i == 1:
            return self.coords[i] + self.extension[i]

        # should raise an out of bounds exception
        return self.coords[i]

    @property
    def start(self):
        return self.coords[0]

    @property
    def stop(self):
        return self.coords[1]

    def __eq__(self, other):
        return self.refname == other.refname and self.coords == other.coords

    def __repr__(self):
        return '{0}:{1}-{2}'.format(self.refname, self[0],
                                    self[1])

    def to_string_noext(self):
        return '{0}:{1}-{2}'.format(self.refname, self.coords[0],
                                    self.coords[1])

    @staticmethod
    def from_string(intron_string, extend_start = 0, extend_stop = 0):
        ref, rest = intron_string.split(":")
        start, stop = rest.split("-")

        return Intron(ref, int(start) + extend_start, int(stop) - extend_stop,
                      extend_start, extend_stop)

def get_introns(trans, extend = 0):
    """ Given a transcript, return a set of introns with (start, stop)
    locations. If extend > 0, then extends the intron on the left and right
    side.  """

    if extend < 0:
        raise Exception("Non-sensical value for extend (must be >= 0)")

    introns = []
    for i in xrange(len(trans.exons) - 1):
        intron = Intron(trans.refname, trans.exons[i][1],
                        trans.exons[i + 1][0], extend, extend)
        introns.append( intron )

    return introns

def intron_all_trans(trans_list):
    """ Given a list of transcripts, return a sorted list of introns that every
    transcript shares. """

    if len(trans_list) == 0:
        return []
    all_introns = map(get_introns, trans_list)
    print all_introns
    # TODO: fixme
    all_introns = [introns.coords for introns in all_introns]
    all_introns = [set(introns) for introns in all_introns]
    all_introns = list(reduce(set.intersection, all_introns))
    all_introns.sort()

    return all_introns

def intron_intersection(i1, i2):
    """ Given two introns, finds their maximal overlap """
    left = max(i1[0], i2[0])
    right = min(i1[1], i2[1])
    if right <= left:
        return None
    return (left, right)

# @profile
def transcript_union(trans_list):
    """ Given a list of transcripts, return a transcript that is the 'union' of
    them. That is, take the union overlap region of every exon."""

    all_exons = [exon for trans in trans_list for exon in trans.exons]
    all_exons = sorted(list(set(all_exons)))

    exon_union = []
    candidate = all_exons[0]
    for it in xrange(1, len(all_exons)):
        cur_exon = all_exons[it]
        if candidate[1] < cur_exon[1] and cur_exon[0] < candidate[1]:
            candidate = (candidate[0], cur_exon[1])
        elif candidate[1] <= cur_exon[0]:
            if candidate[1] == cur_exon[0]:
                candidate = (candidate[0], cur_exon[1])
            else:
                exon_union.append( candidate )
                candidate = cur_exon
    exon_union.append(candidate)

    t0 = trans_list[0]
    strand = '-' if t0.is_reverse else '+'
    t = Transcript(t0.gene_id,
                   t0.refname,
                   strand,
                   t0.frame,
                   t0.gene_id_attributes,
                   t0.gene_id,
                   t0.score,
                   t0.source)
    t.exons = exon_union

    return t

def intron_trans_compat(intron_list, trans_list):
    """ Under the assumption that the intron_list is derived from the
    union gene from the trans_list, returns a dictionary where intron =>
    list of transcript ids that it is compatible with """

    intron_to_trans = {}
    for intron in intron_list:
        key = str(intron)
        for trans in trans_list:
            if trans.front_coordinate < intron[0] and \
                    intron[1] < trans.end_coordinate and \
                    intron.refname == trans.refname:
                matches = intron_to_trans.get(key, [])
                matches.append(trans.transcript_id)
                intron_to_trans[key] = matches

    return intron_to_trans

def get_overlapping_transcripts(transcripts):
    """ Iterate through a list of transcripts. Every yield returns a list of
    transcripts that overlap in genomic coordinates. """

    transcripts = sorted(transcripts, key = lambda x: (x.refname,
                                                       x.front_coordinate,
                                                       x.end_coordinate))
    cur_end = None
    cur_ref = None
    trans_list = []
    for t in transcripts:
        if cur_end is not None and t.front_coordinate < cur_end and \
                t.refname == cur_ref:
            cur_end = max(cur_end, t.end_coordinate)
        else:
            if cur_end is not None:
                yield trans_list
            cur_end = t.end_coordinate
            cur_ref = t.refname
            trans_list[:] = []
        trans_list.append(t)

    yield trans_list

def discard_overlapping_introns(transcripts, extend = 0):
    """ Given a dictionary which maps from gene_to_union (reduce_to_gene(trans,
    transcript_union)), and gene_to_intron, remove introns that overlap with
    other genes. Side effects: gene_to_introns has an updated set of introns """

    gene_to_introns = {}

    for overlap_trans in get_overlapping_transcripts( transcripts ):
        # get the gene union then the introns from this union
        g2t, g2u = reduce_to_gene( overlap_trans, transcript_union )
        g2i = { gene: get_introns(g_union, extend) for (gene, g_union) in \
               g2u.iteritems() }

        t_unions = sorted(g2u.values(), key = lambda x: (x.front_coordinate,
                              x.end_coordinate))
        intersection = []
        for j in xrange(len(t_unions)):
            for k in xrange(j + 1, len(t_unions)):
                if t_unions[k].front_coordinate < t_unions[j].end_coordinate:
                    l_isect = t_unions[k].front_coordinate
                    r_isect = min(t_unions[j].end_coordinate,
                                  t_unions[k].end_coordinate)
                    intersection.append( (l_isect, r_isect) )
                else:
                    break

        # given a set of candidate intersections, again, find the union of
        # intersection
        # print 'the intersection', intersection
        intersection = unionize_regions(intersection)
        # print 'the intersection', intersection

        # we have a list of intersection regions, now need to find introns
        # that cross intersection
        for cur_gene in t_unions:
            valid_introns = []
            for intron in g2i[cur_gene.gene_id]:
                valid_intron = True
                for invalid_region in intersection:
                    # print 'int: ', intron, ' ir ', invalid_region
                    if (intron[1] <= invalid_region[1] and \
                            invalid_region[0] < intron[1]) or \
                            (invalid_region[1] <= intron[1] and \
                             intron[0] < invalid_region[1]):
                        valid_intron = False
                        break
                if valid_intron:
                    valid_introns.append( intron )

            # print valid_introns
            introns = gene_to_introns.get(cur_gene.gene_id, [])
            introns.extend( valid_introns )
            gene_to_introns[cur_gene.gene_id] = introns

    return gene_to_introns

def unionize_regions(regions):
    """ Given a set of (sorted) regions, take the union of the overlaps """

    cand = None
    unionized = []
    for reg in regions:
        if cand is None:
            cand = reg
        elif reg[0] <= cand[1]:
            cand = (cand[0], max(reg[1], cand[1]))
        else:
            unionized.append(cand)
            cand = reg
    if cand is not None:
        unionized.append(cand)

    return unionized

def intron_all_junction_left(trans_list):
    """ Given a list of transcripts, return a sorted list of intronic regions
    that every transcript shares. """

    if len(trans_list) == 0:
        return []

    trans_list = sorted(trans_list, key = lambda x: x.front_coordinate)


    all_introns = [get_introns(trans) for trans in trans_list]

    if len(all_introns) == 1:
        return all_introns

    # TODO: get the union of all introns and see if it still works -- should be
    # more accurate

    candidates_l = all_introns[0]
    mark_for_removal = []

    # iterate each transcript
    for i in xrange(1, len(all_introns)):
        # XXX: might need to move this into the next loop
        c_start = 0
        for intron in all_introns[i]:
            for j in xrange(c_start, len(candidates_l)):
                cand = candidates_l[j]
                intersection = intron_intersection(intron, cand)
                if intron[0] == cand[0]:
                    candidates_l[j] = intersection
                    c_start += 1
                elif intersection is not None:
                    if trans_list[i].front_coordinate > cand[0]:
                        candidates_l[j] = (cand[0],
                                           trans_list[i].front_coordinate)
                    else:
                        mark_for_removal.append(j)
                elif intron[0] < cand[0] and intron[1] > cand[0]:
                    # there is an overlap from an intron further on the left
                    # side. shouldn't include it from the left side but from the
                    # right
                    mark_for_removal.append(j)
                elif cand[0] > intron[1]:
                    break

        candidates_l = [ c for idx, c in enumerate(candidates_l)
                        if idx not in mark_for_removal]
        mark_for_removal[:] = []

    # cleanup the introns that overlap an exonic region
    for trans in trans_list:
        for j, cand in enumerate(candidates_l):
            if trans.compatible(cand[0]) is not None or \
                trans.compatible(cand[1] - 1) is not None:
                mark_for_removal.append(j)

    candidates_l = [ c for idx, c in enumerate(candidates_l)
                    if idx not in mark_for_removal]

    return candidates_l

# mapping: dictionary where keys are gene names and values is a list of
# transcripts intersection: dictionary where keys are gene names and values are
# intronic regions that are common amonst all isoforms
# def intron_retained_transcripts(mapping, intersection):


class IntronCoverage:
    def __init__(self, ref = None, coords = None, coverage = 0,
            support = (0, 0)):
        self.ref = ref
        self.coords = coords
        self.coverage = coverage
        self.support = support

def junction_support( ps_handle, ref, intron, read_len ):
    # TODO: test me

    left_start = (intron[0] - 1) - read_len + 1
    left_end = (intron[0] - 1) + read_len - 2

    try:
        left_count = sum(read.rlen == read.overlap(left_start, left_end)
                for read in ps_handle.fetch(ref, left_start, left_end))
    except ValueError:
        return (0, 0)

    right_start = (intron[1] - 1) - read_len + 2
    right_end = (intron[1] - 1) + read_len - 1

    right_count = sum(read.rlen == read.overlap(right_start, right_end)
            for read in ps_handle.fetch(ref, right_start, right_end))

    return (left_count, right_count)



def compute_coverage( ps_handle, ref, intron, read_len ):
    # TODO: test me
    left_start = intron[0] - read_len + 1
    right_end = (intron[1] - 1) + read_len - 1

    try:
        count = sum(read.rlen == read.overlap(left_start, right_end) and
                    (read.aend - read.pos) == read.rlen for read in
                    ps_handle.fetch(ref, left_start, right_end))
        cov = float(count) / (right_end - left_start + 1)
    except ValueError:
        return 0.0

    return cov

# bam_fname: a BAM file name
#
def bam_to_measurable(bam_fname, gene_to_trans, gene_to_introns):
    # TODO: test me

    bam_handle = pysam.Samfile(bam_fname, 'rb')
    tmp = bam_handle.next()
    read_len = tmp.rlen

    gene_to_max_introns = {}
    gene_to_measurable_introns = {}

    for gene_name, all_introns in gene_to_introns.iteritems():
        if len(all_introns) == 0:
            continue
        max_intron = None
        max_cov = 0.0

        measurable_introns = []

        ref = gene_to_trans[gene_name][0].refname
        for intron in all_introns:
            cur_cov = None
            junc_supp = junction_support(bam_handle, ref, intron, read_len)
            cov = compute_coverage(bam_handle, ref, intron, read_len)
            cur_cov = IntronCoverage(ref, intron, cov, junc_supp)

            if cov > max_cov:
                max_intron = cur_cov
                max_cov = cov
            measurable_introns.append( cur_cov )

            # For now compute for every intron
            # if sum(junc_supp) > 0:
            #     cov = compute_coverage(ref, intron)
            #     cur_cov = IntronCoverage(ref, intron, cov, True)
            # else:
            #     cur_cov = IntronCoverage(ref, intron)

        gene_to_max_introns[ gene_name ] = max_intron
        gene_to_measurable_introns[ gene_name ] = measurable_introns

    bam_handle.close()

    return (gene_to_max_introns, gene_to_measurable_introns)

def print_measurable(gene2trans, gene2intersect, gene2max, gene2measurable, handle):
    print >> handle, 'gene\treference\tstart\tend\tcoverage\tsupport_start\tsupport_end\ttranscripts'
    for gene in gene2max:
        for measurable in gene2measurable[gene]:
            cur_line = []
            cur_line.append(gene)
            cur_line.append( measurable.ref )
            cur_line.append( str(measurable.coords[0]) )
            cur_line.append( str(measurable.coords[1]) )
            cur_line.append( str(measurable.coverage) )
            cur_line.append( str(measurable.support[0]) )
            cur_line.append( str(measurable.support[1]) )
            cur_line.append( ','.join(trans.transcript_id
                for trans in gene2trans[gene]) )
            print >> handle, '\t'.join( cur_line )
