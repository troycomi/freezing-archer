import argparse
import sys
from collections import defaultdict
import locale
import gc

from myBedTools3 import myBedTools
import pandas as pd

# in original mybedtools, I allowed the user to specify the
# ref version with namespace.ref_version
# here I've replaced all of those with:
local_ref_version = 'b37'


class BinarySeqFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binaryseqfilegenome,
                        initialize=False, ref_version=local_ref_version)
        bt.read_as_regions(filename, file_type=myBedTools.binaryseqfilegenome)
        setattr(namespace, self.dest, bt)


class missingdict(defaultdict):
    # missingdict allows us to not add something to the
    # defaultdict just because we query it
    def __missing__(self, key):
        return self.default_factory()


class BinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binarybedfilegenome,
                        initialize=False, ref_version=local_ref_version)
        bt.read_as_regions(filename)
        setattr(namespace, self.dest, bt)


class MultBinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        def fn_outer(i):
            bt = myBedTools(myBedTools.binarybedfilegenome,
                            initialize=False, ref_version=local_ref_version)
            bt.read_as_regions(i)
            return bt
        setattr(namespace, self.dest, [fn_outer(i) for i in filenames])


class MergeBinaryBedFilesAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        sys.stderr.write("Merging %d bbg files\n" % len(filenames))
        bt = myBedTools(myBedTools.binarybedfilegenome,
                        initialize=False, ref_version=local_ref_version)
        for f in filenames:
            bt.read_as_regions(f)
        setattr(namespace, self.dest, bt)


class IntersectBinaryBedFilesAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        sys.stderr.write("intersecting %d bbg files\n" % len(filenames))
        bt = myBedTools(myBedTools.binarybedfilegenome,
                        initialize=False, ref_version=local_ref_version)
        bt.read_to_bases(myBedTools.binarybedfilegenome,
                         filenames[0], myBedTools.set_to_one)
        for f in filenames[1:]:
            bt.read_to_bases(myBedTools.binarybedfilegenome,
                             f, myBedTools.bitfn_and)
            gc.collect()
        setattr(namespace, self.dest, bt)


def munge_regions(opts):
    if opts.regions is not None:
        sys.stderr.write('Initial number of bases considered (regions): %s\n'
                         % locale.format('%d', opts.regions.total_length(),
                                         grouping=True))

    if opts.regions is not None and opts.intersect_region is not None:
        opts.regions.bases &= opts.intersect_region.bases
        opts.intersect_region = None

    if opts.regions is not None and opts.exclude_region is not None:
        sys.stderr.write(
            'Excluding %s bases.\n' %
            locale.format('%d', opts.exclude_region.total_length(),
                          grouping=True))
        opts.regions.bases &= ~opts.exclude_region.bases
        opts.exclude_region = None

    elif opts.regions is None and opts.exclude_region is not None:
        sys.stderr.write(
            'Excluding %s bases.\n' %
            locale.format('%d', opts.exclude_region.total_length(),
                          grouping=True))
        opts.regions = opts.exclude_region
        opts.exclude_region = None
        opts.regions.bases.invert()

    if opts.regions is not None:
        sys.stderr.write(
            'Final number of bases considered (regions, exclude, '
            'intersect, etc): %s\n\n' %
            locale.format('%d', opts.regions.total_length(), grouping=True))

    if opts.regions is None and not opts.vcf_is_ms_file:
        sys.stderr.write('NO MASKING GIVEN - ASSUMING WHOLE GENOME IS PERFECT '
                         '(hint: this is probably not a good assumption)\n\n')
