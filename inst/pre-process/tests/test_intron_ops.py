# intron ops
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


import unittest

from gtf_parser import *
from intron_ops import *

class TestIntronOps(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_reduce_to_gene(self):
        gtf_dict = gtf_parse('/Users/hjp/lmcb/unittests/3gene.gtf')
        gtf_list = sorted(gtf_dict.values(), key = lambda x: x.front_coordinate)
        g2t = reduce_to_gene(gtf_list)

        self.assertEqual(len(g2t), 3)
        self.assertEqual(len(g2t['ENSG00000124193']), 2)
        self.assertEqual(len(g2t['ENSG00000078399']), 5)
        self.assertEqual(len(g2t['ENSG00000151303']), 1)

        g2t, gene_len = reduce_to_gene(gtf_list, lambda x: len(x))
        self.assertEqual(gene_len['ENSG00000124193'], 2)
        self.assertEqual(gene_len['ENSG00000078399'], 5)
        self.assertEqual(gene_len['ENSG00000151303'], 1)

    def test_get_introns(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 10))
        t1.add_exon((15, 20))

        introns = get_introns(t1)
        self.assertEqual( len(introns), 1)
        self.assertEqual( introns[0], Intron('chr', 10, 15) )

        t1.add_exon((21, 30) )
        introns = get_introns(t1)

        self.assertEqual(introns, [Intron('chr', 10, 15),
                                   Intron('chr', 20, 21)])

    # FIXME: after class Intron was introduced, broke this... not too difficult
    # to fix..
    # def test_intron_all_trans(self):
    #     t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
    #     t1.add_exon((3, 10))
    #     t1.add_exon((15, 20))

    #     t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
    #     t2.add_exon((3, 10))
    #     t2.add_exon((15, 20))
    #     t2.add_exon((21, 30) )

    #     introns = intron_all_trans([t1, t2])

    #     self.assertEqual(introns, [Intron('chr', 10, 15)])

    #     self.assertEqual(intron_all_trans([]), [])

    def test_intron_intersection(self):
        i1 = (3, 10)
        i2 = (2, 8)
        self.assertEqual(intron_intersection(i1, i2), (3,8))

        i2 = (10, 12)
        self.assertEqual(intron_intersection(i1, i2), None)

    # def test_intron_all_junction_left(self):
    #     t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
    #     t1.add_exon((3, 10))
    #     t1.add_exon((15, 20))

    #     t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
    #     t2.add_exon((3, 10))
    #     t2.add_exon((18, 25))

    #     intronic_regions = intron_all_junction_left([t1, t2])
    #     self.assertEqual(intronic_regions, [(10, 15)])

    #     t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
    #     t3.add_exon((3, 10))
    #     t3.add_exon((16, 23))

    #     intronic_regions = intron_all_junction_left([t1, t2, t3])
    #     self.assertEqual(intronic_regions, [(10, 15)])

    #     # this intron shouldn't match on the left side
    #     t1.add_exon((40, 44))
    #     t2.add_exon((40, 44))
    #     t3.add_exon((40, 44))

    #     intronic_regions = intron_all_junction_left([t1, t2, t3])
    #     self.assertEqual(intronic_regions, [(10, 15)])

    #     # this intron should match on the left side
    #     t1.add_exon((49, 57))
    #     t2.add_exon((47, 53))
    #     t3.add_exon((47, 54))
    #     intronic_regions = intron_all_junction_left([t1, t2, t3])
    #     self.assertEqual(intronic_regions, [(10, 15),
    #                                         (44, 47)])

    # def test_intron_all_junction_left_exon_overlap(self):
    #     # test case where exon overlaps introns
    #     t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
    #     t1.add_exon((3, 10))
    #     t1.add_exon((15, 20))

    #     t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
    #     t2.add_exon((3, 10))
    #     t2.add_exon((18, 25))

    #     t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
    #     t3.add_exon((3, 25))

    #     intronic_regions = intron_all_junction_left([t1, t2, t3])
    #     self.assertEqual(intronic_regions, [])

    # def test_intron_all_junction_left_no_overlap(self):
    #     t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
    #     t1.add_exon((3, 10))
    #     t1.add_exon((15, 20))
    #     t1.add_exon((30, 40))

    #     t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
    #     t2.add_exon((18, 25))
    #     t2.add_exon((30, 40))

    #     t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
    #     t3.add_exon((18, 25))
    #     t3.add_exon((30, 40))

    #     intronic_regions = intron_all_junction_left([t1, t2, t3])
    #     self.assertEqual(intronic_regions, [(10, 15)])

    #     intronic_regions = intron_all_junction_left([t2, t1, t3])
    #     self.assertEqual(intronic_regions, [(10, 15)])

    # def test_intron_all_junction_left_gtf(self):
    #     trans = gtf_parse('tests/inputs/refGene_07.23.2014_CHL1.gtf')
    #     trans = sorted(trans.values(), key = lambda x: x.front_coordinate)

    #     gene_to_trans, gene_to_introns = \
    #         reduce_to_gene(trans, intron_all_junction_left)
    #     self.assertEqual(gene_to_introns['CHL1'][0], (238746, 239325))
    #     self.assertEqual(gene_to_introns['CHL1'][1], (286375, 288250))
    #     self.assertEqual(gene_to_introns['CHL1'][2], (361550, 367641))
    #     self.assertEqual(gene_to_introns['CHL1'][3], (367747, 369849))
    #     self.assertEqual(gene_to_introns['CHL1'][4], (370037, 382476))
    #     self.assertEqual(gene_to_introns['CHL1'][5], (382599, 383594))
    #     self.assertEqual(gene_to_introns['CHL1'][6], (383765, 384666))
    #     self.assertEqual(gene_to_introns['CHL1'][7], (386392, 391041))
    #     self.assertEqual(gene_to_introns['CHL1'][8], (391226, 396322))
    #     self.assertEqual(gene_to_introns['CHL1'][9], (396454, 401966))
    #     self.assertEqual(gene_to_introns['CHL1'][10], (402107, 403381))
    #     self.assertEqual(gene_to_introns['CHL1'][11], (403493, 404899))
    #     self.assertEqual(gene_to_introns['CHL1'][12], (405066, 407632))
    #     self.assertEqual(gene_to_introns['CHL1'][13], (407798, 419500))
    #     self.assertEqual(gene_to_introns['CHL1'][14], (419625, 423861))
    #     self.assertEqual(gene_to_introns['CHL1'][15], (423963, 424156))
    #     self.assertEqual(gene_to_introns['CHL1'][16], (424354, 425498))
    #     self.assertEqual(gene_to_introns['CHL1'][17], (425569, 430934))
    #     self.assertEqual(gene_to_introns['CHL1'][18], (431157, 432383))
    #     self.assertEqual(gene_to_introns['CHL1'][19], (432499, 432637))
    #     self.assertEqual(gene_to_introns['CHL1'][20], (432842, 433357))
    #     self.assertEqual(gene_to_introns['CHL1'][21], (433480, 436375))
    #     self.assertEqual(gene_to_introns['CHL1'][22], (436555, 439909))
    #     self.assertEqual(gene_to_introns['CHL1'][23], (440831, 443308))
    #     self.assertEqual(gene_to_introns['CHL1'][24], (443381, 447177))

    def test_get_introns_gtf(self):
        trans = gtf_parse('tests/inputs/refGene_07.23.2014_CHL1.gtf')
        introns = get_introns(trans['NM_001253388'])
        self.assertEqual(len(introns), 24)

    def test_transcript_union_simple2(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 10))
        t1.add_exon((15, 20))
        t1.add_exon((30, 40))

        t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
        t2.add_exon((18, 25))
        t2.add_exon((30, 40))

        tu = transcript_union([t1, t2]).exons
        self.assertEqual(tu, [(3, 10), (15, 25), (30, 40)])

        tu = transcript_union([t2, t1]).exons
        self.assertEqual(tu, [(3, 10), (15, 25), (30, 40)])

        t2.add_exon((45, 50))
        tu = transcript_union([t1, t2]).exons
        self.assertEqual(tu, [(3, 10), (15, 25), (30, 40), (45, 50)])

        t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
        t3.add_exon((18, 25))
        t3.add_exon((30, 40))

        tu = transcript_union([t1, t2, t3]).exons
        self.assertEqual(tu, [(3, 10), (15, 25), (30, 40), (45, 50)])

        t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
        t3.add_exon((12, 25))
        t3.add_exon((30, 40))

        tu = transcript_union([t1, t2, t3]).exons
        self.assertEqual(tu, [(3, 10), (12, 25), (30, 40), (45, 50)])

        t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
        t3.add_exon((9, 25))
        t3.add_exon((30, 40))

        tu = transcript_union([t1, t2, t3]).exons
        self.assertEqual(tu, [(3, 25), (30, 40), (45, 50)])

    def test_transcript_union_no_intron(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 10))
        t1.add_exon((15, 20))

        t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
        t2.add_exon((10,15))

        tu = transcript_union([t1, t2]).exons
        self.assertEqual(tu, [(3,20)])

    def test_transcript_union_internal_exon(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 10))
        t1.add_exon((15, 20))

        t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
        t2.add_exon((11,15))

        tu = transcript_union([t1, t2]).exons
        self.assertEqual(tu, [(3,10), (11, 20)])

    def test_transcript_union_disjoint(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 10))
        t1.add_exon((15, 20))

        t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
        t2.add_exon((25, 30))
        t2.add_exon((40, 45))

        tu = transcript_union([t1, t2]).exons
        self.assertEqual(tu, [(3,10), (15, 20), (25, 30), (40, 45)])

    def test_transcript_union_gtf(self):
        trans = gtf_parse('tests/inputs/refGene_07.23.2014_CHL1.gtf')
        trans = sorted(trans.values(), key = lambda x: x.front_coordinate)

        gene_to_trans, gene_to_union = \
            reduce_to_gene(trans, transcript_union)
        chl1 = gene_to_union['CHL1'].exons
        self.assertEqual(chl1[0], (238278, 238746))
        self.assertEqual(chl1[1], (239325, 239775))
        self.assertEqual(chl1[2], (286295, 286375))
        self.assertEqual(chl1[3], (288250, 290282))
        self.assertEqual(chl1[4], (361365, 361550))
        self.assertEqual(chl1[5], (367641, 367747))
        self.assertEqual(chl1[6], (369849, 370037))
        self.assertEqual(chl1[7], (382476, 382599))
        self.assertEqual(chl1[8], (383594, 383765))
        self.assertEqual(chl1[9], (384666, 384714))
        self.assertEqual(chl1[10], (386271, 386392))
        self.assertEqual(chl1[11], (391041, 391226))
        self.assertEqual(chl1[12], (396322, 396454))
        self.assertEqual(chl1[13], (401966, 402107))
        self.assertEqual(chl1[14], (403381, 403493))
        self.assertEqual(chl1[15], (404899, 405066))
        self.assertEqual(chl1[16], (407632, 407798))
        self.assertEqual(chl1[17], (419500, 419625))
        self.assertEqual(chl1[18], (423861, 423963))
        self.assertEqual(chl1[19], (424156, 424354))
        self.assertEqual(chl1[20], (425498, 425569))
        self.assertEqual(chl1[21], (430934, 431157))
        self.assertEqual(chl1[22], (432383, 432499))
        self.assertEqual(chl1[23], (432637, 432842))
        self.assertEqual(chl1[24], (433357, 433480))
        self.assertEqual(chl1[25], (436375, 436555))
        self.assertEqual(chl1[26], (439909, 440068))
        self.assertEqual(chl1[27], (440699, 440831))
        self.assertEqual(chl1[28], (443308, 443381))
        self.assertEqual(chl1[29], (447177, 451097))

    def test_transcript_union_start(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 10))
        t1.add_exon((20, 25))
        t1.add_exon((30, 35))

        t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
        t2.add_exon((12, 17))
        t2.add_exon((20, 25))

        tu = transcript_union([t1, t2]).exons
        self.assertEqual(tu, [(3, 10), (12, 17), (20, 25), (30, 35)])

    def test_transcript_union_clk1(self):
        trans = gtf_parse('tests/inputs/refGene_07.23.2014_CLK1.gtf')
        trans = sorted(trans.values(), key = lambda x: x.front_coordinate)

        gene_to_trans, gene_to_union = \
            reduce_to_gene(trans, transcript_union)
        clk1 = gene_to_union['CLK1'].exons

        self.assertEqual(clk1[8], (201724402, 201726189))
        self.assertEqual(clk1[9], (201726424,201726585))
        self.assertEqual(clk1[-1], (201729286, 201729467))
        self.assertEqual(clk1[-2], (201728823, 201729284))
        self.assertEqual(len(clk1), 12)

    def test_intron_from_string(self):
        i1_str = "chr1:30-40"
        i1 = Intron.from_string(i1_str)

        self.assertEqual(i1[0], 30)
        self.assertEqual(i1[1], 40)

        i2 = Intron("chrX", 10, 100, 3, 2)
        i2_from_str = Intron.from_string(str(i2), 3, 2)

        self.assertEqual(str(i2), str(i2_from_str))
        self.assertEqual(i2.to_string_noext(), i2_from_str.to_string_noext())

class TestTransOps(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_iterate_overlapping_transcripts(self):
        t1 = Transcript('t1', 'chr', '+', '.', None, 'g1')
        t1.add_exon((3, 35))

        t2 = Transcript('t2', 'chr', '+', '.', None, 'g1')
        t2.add_exon((12, 17))
        t2.add_exon((20, 25))

        t3 = Transcript('t3', 'chr', '+', '.', None, 'g1')
        t3.add_exon((100, 120))

        it = get_overlapping_transcripts([t3, t2, t1])
        self.assertEqual(it.next(), [t1, t2])
        self.assertEqual(it.next(), [t3])
        self.assertRaises(StopIteration, it.next)

class TestIntronTransCompat(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_intron_trans_compat_simple(self):
        t1 = Transcript('t1', 'chr1', '+', '.', None, 'g1')
        t1.add_exon((3, 35))
        t1.add_exon((39, 50))

        t2 = Transcript('t2', 'chr2', '+', '.', None, 'g1')
        t2.add_exon((3, 35))
        t2.add_exon((39, 50))

        i1 = Intron('chr1', 35, 39)

        compat = intron_trans_compat([i1], [t2, t1])
        self.assertEqual(compat, {str(i1) : ['t1']})

    def test_intron_trans_compat(self):
        # TODO: more testing required !
        trans = gtf_parse('tests/inputs/refGene_07.23.2014_CHL1.gtf')
        trans = sorted(trans.values(), key = lambda x: x.front_coordinate)

        gene_to_trans, gene_to_union = \
            reduce_to_gene(trans, transcript_union)

        introns = get_introns(gene_to_union['CHL1'])
        i2t = intron_trans_compat(introns, gene_to_trans['CHL1'])
        self.assertEqual(len(i2t), len(introns))

    def test_discard_overlapping_introns_simple(self):
        g1 = Transcript('g1', 'chr1', '+', '.', None, 'g1')
        g1.add_exon( (0, 3) )
        g1.add_exon( (7, 10) )
        g1.add_exon( (12, 35) )
        g1.add_exon( (40, 100) )

        g2 = Transcript('g2', 'chr1', '+', '.', None, 'g2')
        g2.add_exon( (0, 4) )
        g2.add_exon( (6, 9) )
        g2.add_exon( (15, 30) )

        g3 = Transcript('g3', 'chr1', '+', '.', None, 'g3')
        g3.add_exon( (80, 82) )
        g3.add_exon( (87, 90) )

        trans = [g2, g1, g3]
        g2i = discard_overlapping_introns(trans)
        self.assertEqual(g2i['g1'], [Intron('chr1',35, 40)])
        self.assertEqual(len(g2i['g2']), 0)
        self.assertEqual(len(g2i['g3']), 0)

        g1.exons[3] = (40, 60)

        g2i = discard_overlapping_introns(trans)
        self.assertEqual(g2i['g1'], [Intron('chr1', 35, 40)])
        self.assertEqual(len(g2i['g2']), 0)
        self.assertEqual(g2i['g3'], [Intron('chr1', 82, 87)])

        # TODO: test with empty introns

    def test_discard_overlapping_introns_partial_overlap(self):
        g1 = Transcript('g1', 'chr1', '+', '.', None, 'g1')
        g1.add_exon( (0, 3) )
        g1.add_exon( (7, 10) )
        g1.add_exon( (12, 20) )

        g2 = Transcript('g2', 'chr1', '+', '.', None, 'g2')
        g2.add_exon( (70, 72) )
        g2.add_exon( (74, 81) )
        g2.add_exon( (83, 85) )

        g3 = Transcript('g3', 'chr1', '+', '.', None, 'g3')
        g3.add_exon( (80, 82) )
        g3.add_exon( (87, 90) )

        trans = [g2, g1, g3]
        g2i = discard_overlapping_introns(trans)
        self.assertEqual(g2i['g1'], [Intron('chr1', 3, 7),
                                     Intron('chr1', 10, 12)])
        self.assertEqual(g2i['g2'], [Intron('chr1', 72, 74)])
        self.assertEqual(g2i['g3'], [])

    def test_discard_overlapping_introns_same_coords_diff_gene(self):
        g1 = Transcript('g1', 'chr1', '+', '.', None, 'g1')
        g1.add_exon( (0, 30) )
        g1.add_exon( (40, 70) )

        g2 = Transcript('g2', 'chr1', '+', '.', None, 'g2')
        g2.add_exon( (0, 30) )
        g2.add_exon( (40, 70) )

        trans = [g2, g1]

        g2i = discard_overlapping_introns(trans)

        self.assertEqual(g2i['g1'], [])
        self.assertEqual(g2i['g2'], [])

    def test_discard_overlapping_introns_gtfird2(self):
        # In this test there are two genes that have the same exact front and
        # end coordinates for the gene...
        gtf_in = 'tests/inputs/GTF2IRD2.gtf'
        gtf_dict = gtf_parse(gtf_in)
        gtf_list = gtf_dict.values()

        # XXX: seems like GTFIRD2B isn't getting added to this overlap... weird
        g2i = discard_overlapping_introns(gtf_list)
        self.assertEqual(g2i['GTF2IRD2B'], [])

    def test_unionize_regions(self):
        r = [(0, 5), (10, 20)]
        self.assertEqual( unionize_regions(r), r )

        r = [(0, 5), (3, 7), (3, 20)]
        self.assertEqual( unionize_regions(r), [(0, 20)] )

        r = [(0, 5), (3, 7), (10, 20)]
        self.assertEqual( unionize_regions(r), [(0, 7), (10, 20)] )

        r = [(0, 30), (3, 7), (10, 20)]
        self.assertEqual( unionize_regions(r), [(0, 30)] )

        r = [(0, 5), (5, 7)]
        self.assertEqual( unionize_regions(r), [(0, 7)] )

class TestIntronClass(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_intron(self):
        int1 = Intron('sup', 3, 7)
        self.assertEqual(int1[0], 3)
        self.assertEqual(int1[1], 7)

        # self.assertRaises(Exception, int1.__getitem__, 2)

    def test_repr(self):
        int1 = Intron('hi', 2, 4)
        self.assertEqual(str(int1), 'hi:2-4')
        self.assertEqual(int1.to_string_noext(), 'hi:2-4')

        int1 = Intron('hi', 2, 4, 1, 1)
        self.assertEqual(str(int1), 'hi:1-5')
        self.assertEqual(int1.to_string_noext(), 'hi:2-4')

    def test_extension(self):
        i1 = Intron('hi', 2, 4, 1, 2)
        self.assertEqual(i1[0], 1)
        self.assertEqual(i1[1], 6)

