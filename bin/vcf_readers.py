import sys
import gzip
import pandas as pd
import bisect


class archaic_vcf:
    def __init__(self, filename, ancestral_bsg=None,
                 vcf_has_illumina_chrnums=False, regions=None):
        if filename.endswith('.gz'):
            vcffile = gzip.open(filename)
        else:
            vcffile = open(filename, 'r')

        self.vcf = {}
        self.ancestral_remap = {
            '1': '0',
            '.': '0',
            '0': '1'
        }
        self.derived_cache = None

        if filename is not None:
            sys.stderr.write("Reading VCF file %s..\n" % filename)
        c = 0
        warn_strip_chrom = False

        for line in vcffile:
            if line.strip().startswith('#'):
                continue
            try:
                (chrom, pos, _, ref, alt,
                 _, _, _, _, gt_info) = line.strip().split()
            except ValueError:
                raise ValueError(
                    "Too many columns in ARCHAIC VCF: %s?\n" % filename +
                    "Expecting one individual (10 columns):\n" + line)
            try:
                pos = int(pos)
            except ValueError:
                raise ValueError(
                    "Unable to parse position in ARCHAIC VCF: %s\n" % filename
                    + line
                )

            if vcf_has_illumina_chrnums and chrom.startswith('chr'):
                chrom = chrom[3:]
                warn_strip_chrom = True

            if (regions is not None
                    and not regions.in_region_one_based(chrom, pos)):
                continue

            if chrom in self.vcf and pos in self.vcf[chrom]:
                raise ValueError(
                    "error - duplicate position in VCF file?\n" +
                    " ".join([str(s) for s in (chrom, pos, ref, alt)])
                    + '\n' + line)

            if len(alt) > 1:
                raise ValueError(
                    "Error reading archaic VCF file:\n" +
                    filename +
                    "\nlen(alt) > 1: require only "
                    "bi-allelic SNPs in archaic VCF\n" + line)

            gt = gt_info[:3]

            anc = None
            if ancestral_bsg is not None:
                anc = ancestral_bsg.get_base_one_based(chrom, pos)

            self.add_site(chrom, pos, gt, ref, alt, anc)

            c += 1

        if filename is not None:
            sys.stderr.write(" with %d lines.\n" % c)
        if warn_strip_chrom:
            sys.stderr.write(" WARNING: REMOVED LEADING 'chr' FROM "
                             "CHROMOSOME NAMES, TO MATCH ILLUMINA VCFS\n")

    def add_site(self, chrom, pos, gt, ref, alt, ancestral):
        if chrom not in self.vcf:
            self.vcf[chrom] = {}
        if ancestral == alt:
            # 1 -> 0, . -> 0, 0 -> 1
            gt = ''.join([self.ancestral_remap[s]
                          if s in self.ancestral_remap
                          else s for s in gt])
            self.vcf[chrom][pos] = ('0' in gt, gt, ancestral, ref)
        else:
            self.vcf[chrom][pos] = ('0' in gt, gt, ref, alt)

    def has_derived(self, chrom, pos):
        if not self.has_site(chrom, pos):
            return False
        return '1' in self.vcf[chrom][pos][1]

    def get_derived(self, chrom, pos):
        if not self.has_site(chrom, pos):
            return 'N'
        return self.vcf[chrom][pos][3]

    def get_derived_count(self, chrom, pos):
        if not self.has_site(chrom, pos):
            return 0
        return self.vcf[chrom][pos][1].count('1')

    def has_site(self, chrom, pos):
        return chrom in self.vcf and pos in self.vcf[chrom]

    def get_derived_sites(self, chrom, winstart, winend):
        if self.derived_cache is None:
            self.generate_derived_cache()
        if chrom not in self.derived_cache:
            return []
        return self.derived_positions[chrom][
            bisect.bisect_left(self.derived_positions[chrom], winstart):
            bisect.bisect_left(self.derived_positions[chrom], winend)]

    def get_derived_sites_with_der_count(self, chrom, winstart, winend):
        if self.derived_cache is None:
            self.generate_derived_cache()
        if chrom not in self.derived_cache:
            return []
        return self.derived_cache[chrom][
            bisect.bisect_left(self.derived_positions[chrom], winstart):
            bisect.bisect_left(self.derived_positions[chrom], winend)]

    def generate_derived_cache(self):
        self.derived_cache = {}
        self.derived_positions = {}
        for chrom in self.vcf:
            self.derived_cache[chrom] = sorted(
                [(position, entry[1].count('1'))
                 for position, entry in self.vcf[chrom].items()
                 if '1' in entry[1]
                 ],
                key=lambda x: x[0]
            )
            self.derived_positions[chrom] = [
                entry[0] for entry in self.derived_cache[chrom]]


class ancestral_vcf(object):
    def __init__(self, filename):
        '''
        Read in vcf, storing it's chromosome, position and reference allele
        '''
        self.vcf = pd.read_csv(
            filename, sep='\t', comment='#', usecols=[0, 1, 3],
            names=['chromosome', 'position', 'ref'],
            dtype={'chromosome': str, 'position': int, 'ref': str},
            index_col=[0, 1]
        )
        if self.vcf.index.duplicated().any():
            raise ValueError(
                "error - duplicate position in VCF file?\n" +
                str(self.vcf.loc[self.vcf.index.duplicated(keep=False)])
            )

    def get_base_one_based(self, chrom, pos):
        try:
            return self.vcf.loc[chrom].loc[pos]['ref']
        except KeyError:
            return 'N'
