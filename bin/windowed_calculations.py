from __future__ import division
import sys
import random
import argparse
from read_vcf import vcf_to_genotypes_windowed
from read_ms import ms_to_genotypes_windowed
import read_archaic_vcf
import gzip
from custom_argparse import (BinarySeqFileAction, BinaryBedFileAction,
                             IntersectBinaryBedFilesAction,
                             MergeBinaryBedFilesAction,
                             munge_regions)
from vcf_readers import ancestral_vcf
import pandas
import tables
import locale

import s_star_fns
import arc_match_pval_tables
import d_stat_fns
import test_fns
import report_mappable_fns
import random_region_pvals

locale.setlocale(locale.LC_ALL, 'en_US')

parser = argparse.ArgumentParser(description='Calculate S*')

parser.add_argument('-vcf', '--vcf-file', required=False,
                    type=argparse.FileType('r'), default=None)
parser.add_argument('-vcfz', '--gzip-vcf-file', required=False, default=None)
parser.add_argument('-indf', '--ind-pop-file', required=False,
                    type=argparse.FileType('r'), default=None)
parser.add_argument('-ptable', '--match-pval-table', default=None)
parser.add_argument('-ptable-precomputed', '--match-pval-table-precomputed',
                    default=None, type=argparse.FileType('r'))
parser.add_argument('-ptmode', '--table-pval-mode',
                    default='pytables', help='could be pandas')
parser.add_argument('-pninds', '--ptables-ninds-for-sfs',
                    type=int, default=108)
parser.add_argument('-padjsfs', '--ptables-adjust-sfs-from-target-to-ref',
                    type=float, default=1)
parser.add_argument('-pfile', '--pickle-file-for-match-pvals', default=None)
parser.add_argument('-winlen', '--window-length', type=int, default=50000)
parser.add_argument('-winstep', '--window-step', type=int, default=20000)
parser.add_argument('-winfile', '--window-file',
                    type=argparse.FileType('r'), default=None)
parser.add_argument('-winchrom', '--process-chromosome',
                    type=int, default=None)
parser.add_argument('-l', '--limit-wins', type=int, default=0)
parser.add_argument('-p', '--progress', type=int, default=0)
parser.add_argument('-range', '--window-range', type=int,
                    nargs=2, default=None,
                    help='Only consider windows within this position range.',
                    metavar=('START', 'STOP'))
parser.add_argument('-ref-pops', '--reference-populations',
                    default=[], nargs='+')
parser.add_argument('-ref-inds', '--reference-individuals',
                    default=[], nargs='+')
parser.add_argument('-target-pops', '--target-populations',
                    default=[], nargs='+')
parser.add_argument('-target-inds', '--target-individuals',
                    default=[], nargs='+')
parser.add_argument('-exclude-pops', '--exclude-populations',
                    default=[], nargs='+', help='Exclude snps where these '
                    'individuals have a nonref allele (or derived, '
                    'if an ancestral genome is given)')
parser.add_argument('-exclude-inds', '--exclude-individuals',
                    default=[], nargs='+', help='Exclude snps where these '
                    'individuals have a nonref allele (or derived, if an '
                    'ancestral genome is given)')
parser.add_argument('-tag-ids', '--tag-ids', default=[], nargs='+')
parser.add_argument('-tags', '--tags', default=[], nargs='+')
parser.add_argument('-no-pvalues', '--no-pvalues', action='store_true',
                    help='Do not calculate archaic match pvalues for '
                    's-star mode.')
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-dpv', '--debug-pval-lookup', action='store_true')
parser.add_argument('-ms', '--vcf-is-ms-file', action='store_true')
parser.add_argument('-mspops', '--ms-pop-sizes', default=None, nargs='+',
                    type=int, help='This is identical to the -I argument for '
                    'ms. WRT target and reference populations, numbering '
                    'starts from 0.')
parser.add_argument('-msinds', '--ms-num-diploid-inds', default=None, type=int,
                    help='The number of diploid individuals considered. This '
                    'is important because we sometimes simulate a single '
                    'archaic chromosome.')
parser.add_argument('-msarc', '--ms-archaic-populations', default=[],
                    nargs='+', type=int,
                    help='The archaic populations, if simulated.')
parser.add_argument('-msarcjt', '--ms-archaic-populations-join-times',
                    default=[], nargs='+', type=float,
                    help='The archaic population join times with the '
                    '*modern human ancestors*. These are used to identify '
                    'introgressed haplotypes.')
parser.add_argument('-msarc-to-process', '--ms-archaic-population-to-process',
                    default=None, type=int, help='The archaic population to '
                    'process, if more than one is given.  Default is the '
                    'first archaic population!')
parser.add_argument('-msintrbed', '--report-intr-bed', action='store_true',
                    help='Report introgressed haplotypes in bed format, in '
                    'addition to other output.')
parser.add_argument('-mssimlen', '--ms-simulated-region-length', default=None,
                    type=int, help='The number of bases simulated in ms '
                    '(i.e., the second argument to -r).')
parser.add_argument('-illumina-chrom', '--vcf-has-illumina-chrnums',
                    action='store_true')
parser.add_argument('-archaic-vcf', '--archaic-vcf', required=False, nargs='+',
                    help='VCF file listing archaic sites', default=None)
parser.add_argument('-ancbsg', '--ancestral-bsg', action=BinarySeqFileAction,
                    required=False, help='BSG file listing ancestral sites '
                    '(CAnc or just chimp)')
parser.add_argument('-ancvcf', '--ancestral-vcf', required=False,
                    help='VCF file listing ancestral sites (CAnc or just '
                    'chimp).  This is cumbersome, and should only be used '
                    'for testing!')

parser.add_argument('-r', '--regions', action=BinaryBedFileAction,
                    required=False, default=None, help='A bbg file that '
                    'specifies which regions to consider.  Only snps in this '
                    'region are loaded.')
parser.add_argument('-ir', '--intersect-region',
                    action=IntersectBinaryBedFilesAction, nargs='+',
                    required=False, default=None, help='A bbg file that is '
                    'intersected with the --regions file to produce a new '
                    'set of regions to consider.')
parser.add_argument('-x', '--exclude-region', nargs='+',
                    action=MergeBinaryBedFilesAction, required=False,
                    default=None, help='bbg file(s) that specify which '
                    'regions should be excluded from the analysis.  If more '
                    'than one file is given, the files are merged.')
parser.add_argument('-regions-min', '--regions-min-mapped',
                    type=int, default=0)

parser.add_argument('-table-query', '--table-query', nargs='+')
parser.add_argument('-len-eps', '--len-eps', type=int, default=1000)
parser.add_argument('-mapped-eps', '--mapped-eps', type=int, default=1000)

parser.add_argument('-o', '--output-file', type=argparse.FileType('w'),
                    required=False, default=sys.stdout, help='output file')

# s* params
parser.add_argument('-s-star-match-bonus', '--s-star-match-bonus',
                    type=int, default=5000)
parser.add_argument('-s-star-max-mismatch', '--s-star-max-mismatch',
                    type=int, default=5)
parser.add_argument('-s-star-mismatch-penalty', '--s-star-mismatch-penalty',
                    type=int, default=-10000)

parser.add_argument('-random-seed', '--random-seed', type=int, default=None)


analysis_group = parser.add_mutually_exclusive_group(required=True)
analysis_group.add_argument('-s-star', '--s-star',
                            action='store_const', dest='analysis',
                            const='s-star', default='s-star')
analysis_group.add_argument('-test-fns', '--test-fns',
                            action='store_const', dest='analysis',
                            const='test-fns')
analysis_group.add_argument('-report-mappable', '--report-mappable',
                            action='store_const', dest='analysis',
                            const='report-mappable')
analysis_group.add_argument('-random-region-pvals', '--random-region-pvals',
                            action='store_const', dest='analysis',
                            const='random-region-pvals')
analysis_group.add_argument('-match-table', '--make-arc-match-pval-tables',
                            action='store_const', dest='analysis',
                            const='match-table')
analysis_group.add_argument('-d-stats', '--d-statistics',
                            action='store_const', dest='analysis',
                            const='d-stats')


opts = parser.parse_args()
setattr(opts, 'first_line', True)

# set stdout to the -o file, if applicable
sys.stdout = opts.output_file

# set the random seed
random.seed(opts.random_seed)

# select the appropriate set of functions
supported_analyses = {
    's-star': s_star_fns,
    'test-fns': test_fns,
    'report-mappable': report_mappable_fns,
    'random-region-pvals': random_region_pvals,
    'match-table': arc_match_pval_tables,
    'd-stats': d_stat_fns,
}

if opts.analysis not in supported_analyses:
    raise ValueError('Unsupported analysis selected!')

module = supported_analyses[opts.analysis]

# process ancestral vcf if given
if opts.ancestral_vcf is not None and opts.ancestral_bsg is None:
    opts.ancestral_bsg = ancestral_vcf(opts.ancestral_vcf)

# THIS IS SUCH A HACK - ONLY USING ONE ARCHAIC VCF AT THIS POINT,
# SO IF THERE'S MORE THAN ONE...... JUST REMOVE THEM
if opts.archaic_vcf is not None:
    if len(opts.archaic_vcf) > 1:
        print("REMOVING ALL BUT ONE ARCHAIC VCF!!!")
    opts.archaic_vcf = read_archaic_vcf.vcf_class(
        opts.archaic_vcf[0], opts.ancestral_bsg,
        opts.vcf_has_illumina_chrnums, opts.regions)

munge_regions(opts)


if opts.vcf_file is None and opts.gzip_vcf_file is None:
    raise ValueError("Require at least one vcf file option.")

elif opts.vcf_file is not None and opts.gzip_vcf_file is not None:
    raise ValueError("Require exactly one vcf file option.")

if not opts.vcf_is_ms_file and opts.ind_pop_file is None:
    raise ValueError("VCF file requires --ind-pop-file.")

elif opts.vcf_is_ms_file and opts.ms_pop_sizes is None:
    raise ValueError("ms file requires --ms-pop-sizes.")

elif opts.vcf_is_ms_file and opts.ms_num_diploid_inds is None:
    raise ValueError("ms file requires --ms-num-diploid-inds.")

elif opts.vcf_is_ms_file and opts.ms_simulated_region_length is None:
    raise ValueError("ms file requires --ms-simulated-region-length.")

elif len(opts.ms_archaic_populations) != \
        len(opts.ms_archaic_populations_join_times):
    raise ValueError(
        "-msarc and -msarcjt must both be given if one is:\n"
        "-msarc:" + opts.ms_archaic_populations + '\n'
        "-msarcjt:" + opts.ms_archaic_populations_join_times)


if opts.gzip_vcf_file is not None:
    opts.vcf_file = gzip.open(opts.gzip_vcf_file)

if opts.vcf_is_ms_file:
    read_genotype_fn = ms_to_genotypes_windowed
else:
    read_genotype_fn = vcf_to_genotypes_windowed

if opts.match_pval_table is not None:

    if opts.table_query is None:
        raise ValueError(
            "must specify what to match when calculating table pvalues..\n"
            'len, mapped, mh, sfs\n'
            '--table-query')

    query_fail = False
    for q in opts.table_query:
        if q not in ['len', 'mapped', 'mh', 'sfs']:
            print('ERROR: invalid match parameter: %s' % q)
            print('must be in: len, mapped, mh, sfs')
            query_fail = True

    if query_fail:
        raise ValueError("query failure")

    if opts.table_pval_mode == 'pytables':
        sys.stderr.write('opening pytables file %s..\n' %
                         opts.match_pval_table)
        opts.match_pval_table = tables.open_file(opts.match_pval_table,
                                                 'r', 'haplotypes')
        sys.stderr.write('...finished opening pytables file\n')

    else:
        sys.stderr.write('loading pval table: %s\n' %
                         opts.match_pval_table)
        dt = pandas.read_table(
            open(opts.match_pval_table, 'rb'), sep=' ',
            header=None, compression='gzip',
            names=['count', 'mapped_bases_bin', 'len',
                   'mh_sites', 'tot_sites', 'sfs', 'std_dev', 'match'])
        dt = dt.query('len > 10000 & mapped_bases_bin > 5000 & tot_sites >= 8')
        opts.match_pval_table = dt
        sys.stderr.write('...finished loading pval tables\n')

    if opts.match_pval_table_precomputed is not None:
        sys.stderr.write('reading precomputed pval file %s..\n' %
                         opts.match_pval_table_precomputed)
        mydict = {}
        c = 1
        header = opts.match_pval_table_precomputed.readline().strip().split()

        for line in opts.match_pval_table_precomputed:
            c += 1

            [count,
             hap_1_window_pval_table,
             hap_1_window_match_pct_table,
             hap_1_window_match_N_table,
             hap_1_window_match_len_table,
             hap_1_window_match_mapped_table,
             hap_1_window_match_mh_table] = line.strip().split()

            # if count == '1': continue

            mydict[(hap_1_window_match_pct_table,
                    hap_1_window_match_N_table,
                    hap_1_window_match_len_table,
                    hap_1_window_match_mapped_table,
                    hap_1_window_match_mh_table)] = hap_1_window_pval_table

            if c % 1000000 == 0:
                print('read..', len(mydict))

        opts.match_pval_table_precomputed = mydict
        sys.stderr.write('...finished reading precomputed pval file\n')


module.initialize_analysis(opts)

nwins = 0
for chrom, winstart, winend, snps in read_genotype_fn(opts.vcf_file,
                                                      opts.window_length,
                                                      opts.window_step,
                                                      opts.ind_pop_file,
                                                      opts):

    nwins += 1
    local_debug = False

    if opts.regions is not None and \
            opts.regions.amount_in_region(chrom, winstart, winend) < \
            opts.regions_min_mapped:
        print("SKIPPING WINDOW", chrom, winstart, winend,
              opts.regions.amount_in_region(chrom, winstart, winend))
        continue

    if opts.debug or local_debug:
        print
    if opts.debug or local_debug:
        print(winstart, winend, len(snps))
    if opts.debug or local_debug:
        print(opts.target_indices)

    if opts.debug or local_debug:
        for snp in snps:
            print('snp', snp['pos'], snp)

    module.run_window_analysis(chrom, winstart, winend, snps, opts)

    if opts.limit_wins > 0 and nwins >= opts.limit_wins:
        break
    if opts.progress > 0 and nwins % opts.progress == 0:
        sys.stderr.write("progress.. %s %d %d | win=%d\n" %
                         (chrom, winstart, winend, nwins))

    if opts.window_range is not None and winstart >= opts.window_range[1]:
        break

module.finish_analysis(opts)
