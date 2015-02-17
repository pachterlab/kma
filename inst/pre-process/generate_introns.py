# Pre-process step for "KeepMeAround"
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

import argparse
import logging
import os
import sys

import gtf_parser
import intron_ops

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyfasta import Fasta

def bed_to_introns(bed_in, fasta_in, fasta_out):
    logging.info("Opening FASTA: {0}".format(fasta_in))
    logging.info("Note: will take a while the first time it is opened.")
    fasta = Fasta(fasta_in)

    bed_h = open(bed_in, 'r')

    all_keys = {}
    output_seq = []
    count = 0
    for line in bed_h:
        if count % 10000 == 0 and count > 0:
            logging.info("On intron: {0}".format(count))
        ref, start, stop, genename = line.split()
        key = ref + ':' + start + '-' + stop
        if key in all_keys:
            logging.warning("ERROR: {0} appears once already".format(key))
            continue
        all_keys[key] = True
        if int(start) > int(stop):
            logging.warning("Intron coords greater than reference: {0}:{1}-{2}".format(ref, start, stop)
)
            logging.warning("Reference length: {0}".format(len(fasta[ref])))
            continue
        seq = fasta[ref][int(start):int(stop)]
        if len(seq) == 0:
            logging.warning("Intron length is 0? {0}:{1}-{2}".format(ref, start, stop))
            continue
        record = SeqRecord(Seq(seq), key, '', '')
        output_seq.append(record)
        count += 1

    bed_h.close()

    with open(fasta_out, 'w') as outf:
        logging.info("Writing intron sequences out to {0}".format(fasta_out))
        SeqIO.write(output_seq, outf, 'fasta')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--extend', type = int, default = 0)
    parser.add_argument('--gtf', required = True,
                        help = 'The transcript file')
    parser.add_argument('--out', required = True,
                        help = 'The directory to output to')
    parser.add_argument('--genome', required = True,
                        help = 'The Multi-FASTA file for the genome')

    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    # check it output directory exists. if not, make it
    if not os.path.exists( args.out ):
        logging.info("Directory '{0}' does not exist -- making it.".format(args.out))
        os.makedirs( args.out )

    if not os.path.isfile( args.genome ):
        logging.error("Genome FASTA '{0}' does not seem to exist".format(args.genome))
        sys.exit(1)

    logging.info('Reading in GTF: {0}'.format(args.gtf))
    gtf_dict = gtf_parser.gtf_parse(args.gtf)
    gtf_list = gtf_dict.values()

    logging.info('Grouping transcripts by gene')
    g2t = intron_ops.reduce_to_gene(gtf_list)

    g2i = intron_ops.discard_overlapping_introns(gtf_list, args.extend)

    bed_out = args.out + os.sep + "introns.bed"

    # TODO: eventually print out GTF union... maybe. might be helpful for
    # debugging
    # print 'printing out gtf'
    # with open(gtf_out, 'w') as outh:
    #     for gene in g2u:
    #         trans = g2u[gene]
    #         outh.write(trans.to_gtf() + '\n')

    i2g = {}
    logging.info("Writing intron BED file: {0}".format(bed_out))
    with open(bed_out, 'w') as outh:
        for gene in g2i:
            for intron in g2i[gene]:
                # for now, output the gene id as the name.. need to come up with
                # better way to do this later for better housekeeping
                intron_list = [intron.refname, str(intron[0]), str(intron[1]),
                               intron.to_string_noext() + ";" + gene]
                i2g[str(intron)] = gene
                outh.write('\t'.join(intron_list) + '\n')

    logging.info('Computing intron-to-transcript compatability')
    i2t = {}
    for gene, introns in g2i.iteritems():
        trans = g2t[gene]
        intron_compat = intron_ops.intron_trans_compat(introns, trans)
        for intron, tlist in intron_compat.iteritems():
            i2t[intron] = tlist

    intron_trans_out = args.out + os.sep + "intron_to_transcripts.txt"
    with open(intron_trans_out, 'w') as outh:
        print >> outh, "intron\ttranscript\tgene\tintron_extension"
        for intron, tlist in i2t.iteritems():
            intron_obj = intron_ops.Intron.from_string(intron, args.extend, args.extend)
            for trans in tlist:
                print >> outh, intron_obj.to_string_noext(), '\t', trans, '\t', i2g[str(intron)], \
                    '\t', intron

    introns_out = args.out + os.sep + "introns.fa"
    # read BED file and FASTA file, write out FASTA file
    bed_to_introns(bed_out, args.genome, introns_out)

if __name__ == '__main__':
    main()
