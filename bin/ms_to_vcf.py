from __future__ import division
import sys, itertools
import argparse
from read_ms import process_ms_block_to_genotypes


import time
start_time = time.time()

import locale
locale.setlocale(locale.LC_ALL, 'en_US')



parser = argparse.ArgumentParser(description='Convert an ms file to a VCF file - in a very basic way.')
parser.add_argument('msfile', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--output-file', type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-winlen', '--window-length', type=int, default=50000)
parser.add_argument('-winstep', '--window-step', type=int, default=20000)
parser.add_argument('-ninds', '--num_inds', type=int, required=True)
parser.add_argument('-arc-num-chrs', '--archaic-num-chrs', type=int, nargs='+', required=False, default=[])
parser.add_argument('-oarc', '--archaic-output-files', nargs='+', type=argparse.FileType('w'), default=None)

opts = parser.parse_args()

if opts.archaic_output_files != None and len(opts.archaic_output_files) != len(opts.archaic_num_chrs):
    print "Require the same number of archaic output files as archaic populations"
    sys.exit()
    pass

if opts.archaic_output_files == None:
    opts.archaic_output_files = [sys.stdout for _ in opts.archaic_num_chrs]
    pass


# ms 42 1 -s 10 
# 44126 40565 42561

# //
# segsites: 10
# positions: 0.1717 0.2230 0.2277 0.4523 0.4598 0.5201 0.7094 0.8533 0.9100 0.9894 
# 0010000001
# 0110000001
# 1000000000
# 1000000000
# 0010000001

# header
header = opts.msfile.readline().strip().split()

# parse params from header
(_, nsam, howmany) = header[:3]
nsam = int(nsam)
# num_inds = nsam // 2
#setattr(opts, 'num_inds', int(nsam) // 2)

if nsam != 2 * opts.num_inds + sum(opts.archaic_num_chrs):
    # sys.stderr.write('Not an even number of chromosomes - dropping last haplotype.\n')
    print "incorrect number of modern human samples / archaic chromosomes specified", nsam, opts.num_inds, opts.archaic_num_chrs
    sys.exit(-1)
    pass

## print header
opts.output_file.write( '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + \
                                      ['i%d' % d for d in xrange(1,opts.num_inds+1)]) + '\n' )

# seed params
seed_params = opts.msfile.readline()

line = opts.msfile.readline()
while not line.startswith('//'):
    if opts.debug: print line
    line = opts.msfile.readline()
    pass

ms_block_results = process_ms_block_to_genotypes(opts.msfile, opts.num_inds, opts.window_length, 'ms1', opts.archaic_num_chrs)

if ms_block_results != None:
    snps, archaic_vcfs = ms_block_results
    
    if opts.debug: 
        snps = list( snps )
        print snps
        pass
    
    blocknum = 1
    for pos, gts in snps:
        opts.output_file.write('\t'.join(str(s) for s in ['ms%d' % blocknum, pos, '.', 'A', 'G', '.', '.', '.', 'GT'] + \
                                             ['%d|%d' % (gt[0], gt[1]) for gt in gts]) + '\n')
        pass

    for vcf_index, archaic_vcf in enumerate(archaic_vcfs):

        opts.archaic_output_files[vcf_index].write( '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 
                                                               'QUAL', 'FILTER', 'INFO', 'FORMAT', 'arc%d' % (vcf_index+1)]) + '\n' )

        # there should only be one chromosome per ms block..
        for chrom in sorted(archaic_vcf.iterkeys()):
            for pos in sorted(archaic_vcf[chrom].iterkeys()):
                gt = archaic_vcf[chrom][pos][1]
                opts.archaic_output_files[vcf_index].write('\t'.join(str(s) for s in [chrom, pos, '.', 'A', 'G', '.', '.', '.', 'GT', gt]) + '\n')
                pass
            pass
        pass

    pass

