import read_vcf
from io import StringIO
import pytest


class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


def test_read_pop_file():
    infile = StringIO(
        u'sample\tpop\tsuper_pop\tgender\n'
        'HG00096\tGBR1\tEUR\tmale\n'
        'HG00097\tGBR2\tEUR\tfemale\n'
        'HG00099\tGBR3\tEUR2\tfemale\n'
    )
    opts = AttributeDict()
    opts.reference_individuals = ['HG00096', 'extra']
    opts.target_individuals = ['HG00097', 'extra']
    opts.exclude_individuals = ['HG00099', 'extra']
    opts.reference_populations = ['GBR1']
    opts.target_populations = ['GBR2']
    opts.exclude_populations = ['GBR3']
    opts.sample_index_in_original_file = {
        'HG00096': 10, 'HG00097': 11, 'HG00099': 12}

    read_vcf.read_1kg_ind_pop_file(infile, opts)

    assert (opts.exclude_individuals == ['HG00099']).all()
    assert opts.exclude_individuals_indexed_to_orig_file == [12]
    # target + reference list
    assert opts.get_id_from_sample_index[0] == 'HG00097'
    assert opts.get_id_from_sample_index[1] == 'HG00096'
    assert opts.get_pop_from_sample_index[0] == 'GBR2'
    assert opts.get_pop_from_sample_index[1] == 'GBR1'
    assert opts.num_reference == 1
    assert opts.num_samples == 2
    assert opts.num_target == 1
    assert (opts.reference_indices == [1]).all()
    assert (opts.reference_individuals == ['HG00096']).all()
    assert opts.reference_individuals_indexed_to_orig_file == [10]
    assert (opts.target_indices == [0]).all()
    assert (opts.target_individuals == ['HG00097']).all()
    assert opts.target_individuals_indexed_to_orig_file == [11]


def test_read_vcf_header():
    opts = AttributeDict()
    infile = StringIO(
        u'##fileformat=VCFv4.2\n'
        '##source=msprime 0.6.1\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##contig=<ID=1,length=10000000>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t'
        'FORMAT\tmsp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
    )
    read_vcf.read_vcf_header(infile, opts)
    assert opts.sample_index_in_original_file == {
        'msp_0': 0,
        'msp_1': 1,
        'msp_2': 2,
        'msp_3': 3,
        'msp_4': 4,
        'msp_5': 5,
    }
    assert infile.readline() == (
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t0|0\n')
    infile = StringIO(
        u'##fileformat=VCFv4.2\n'
        '#chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t'
    )
    with pytest.raises(ValueError) as e:
        read_vcf.read_vcf_header(infile, opts)
    assert str(e.value) == (
        'bad VCF header?\n'
        '#chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t'
    )


def test_process_vcf_line_to_genotypes(mocker):
    opts = AttributeDict()
    opts.vcf_has_illumina_chrnums = False
    opts.window_range = None
    opts.regions = None
    opts.ancestral_bsg = None
    opts.reference_individuals_indexed_to_orig_file = [1, 2]
    opts.target_individuals_indexed_to_orig_file = [0, 4, 5]
    opts.exclude_individuals_indexed_to_orig_file = [3]
    opts.num_target = 3
    opts.num_reference = 2
    opts.archaic_vcf = None
    opts.debug = True
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp == {
        'sfs_target': 1,
        'sfs_reference': 0,
        'target': True,
        'reference': False,
        'chrom': '1',
        'pos': 7,
        'anc': 'A',
        'der': 'T',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [1, 0, 0, 0, 0],
        'haplotypes_2': [0, 0, 0, 0, 0],
        'genotypes': [1, 0, 0, 0, 0],
    }
    opts.window_range = [10]
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t8\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp is None, 'snp before specified range'
    opts.window_range = [1]

    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t9\t.\tAA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp is None, 'snp biallelic'

    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t10\t.\tA\tTT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp is None, 'snp biallelic'

    opts.regions = mocker.MagicMock()
    opts.regions.in_region_one_based.return_value = False
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t11\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp is None, 'snp in region'
    opts.regions.in_region_one_based.return_value = True

    opts.ancestral_bsg = mocker.MagicMock()
    opts.ancestral_bsg.get_base_one_based.return_value = 'N'
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp is None, 'snp not in ancestral'
    opts.ancestral_bsg.get_base_one_based.return_value = 'T'  # match alt

    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t13\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n',
        opts)
    assert snp is None, 'present in excluded (after flip)'

    # test flip map
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t14\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n',
        opts)
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 3,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 14,
        'anc': 'T',
        'der': 'A',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [0, 0, 1, 1, 1],
        'haplotypes_2': [1, 0, 1, 0, 1],
        'genotypes': [1, 0, 2, 1, 2],
    }

    opts.ancestral_bsg.get_base_one_based.return_value = 'C'  # match neither
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t15\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n',
        opts)
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 1,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 15,
        'anc': 'A',
        'der': 'T',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [1, 1, 0, 0, 0],
        'haplotypes_2': [0, 1, 0, 1, 0],
        'genotypes': [1, 2, 0, 1, 0],
    }, 'match neither should be the same as matching ref'

    opts.ancestral_bsg.get_base_one_based.return_value = 'A'  # match ref
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t16\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n',
        opts)
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 1,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 16,
        'anc': 'A',
        'der': 'T',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [1, 1, 0, 0, 0],
        'haplotypes_2': [0, 1, 0, 1, 0],
        'genotypes': [1, 2, 0, 1, 0],
    }

    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t17\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t./.\n',
        opts)
    assert snp is None, 'not in target or reference'

    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t18\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1\t0|0\t1|1\t1|1\n',
        opts)
    assert snp is None, 'fixed in referece'

    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t19\t.\tA\tT\t.\tPASS\t.\tX:GT\t'
        'X:1|1\tX:1|1\tX:1|1\tX:0|0\tX:1|1\tX:1|1\n',
        opts)
    assert snp is None, 'fixed in referece, still'

    with pytest.raises(ValueError) as e:
        read_vcf.process_vcf_line_to_genotypes(
            '1\t19\t.\tA\tT\t.\tPASS\t.\tX:GGT\t'
            'X:1|1\tX:1|1\tX:1|1\tX:0|0\tX:1|1\tX:1|1\n',
            opts)

    assert str(e.value) == (
        'bad format?  GT not found\n'
        'X:GGT'
    )

    opts.archaic_vcf = mocker.MagicMock()
    opts.archaic_vcf.has_site.return_value = True
    opts.archaic_vcf.has_derived.return_value = True
    opts.archaic_vcf.get_derived.return_value = 'T'
    opts.archaic_vcf.get_derived_count.return_value = -1
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n',
        opts)
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 1,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 20,
        'anc': 'A',
        'der': 'T',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [1, 1, 0, 0, 0],
        'haplotypes_2': [0, 1, 0, 1, 0],
        'genotypes': [1, 2, 0, 1, 0],
        'arc_is_derived': True,
        'arc_der_count': -1,
        'arc_match': True,
    }

    opts.ancestral_bsg.get_base_one_based.return_value = 'T'  # match alt
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t21\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n',
        opts)
    opts.archaic_vcf.add_site.assert_not_called()
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 3,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 21,
        'anc': 'T',
        'der': 'A',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [0, 0, 1, 1, 1],
        'haplotypes_2': [1, 0, 1, 0, 1],
        'genotypes': [1, 0, 2, 1, 2],
        'arc_is_derived': False,  # has derived, no match
        'arc_der_count': -1,
        'arc_match': False,
    }

    opts.archaic_vcf.has_site.return_value = False
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t22\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n',
        opts)
    opts.archaic_vcf.add_site.assert_called_once_with(
        '1', 22, '1/1', 'T', 'A', None)
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 3,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 22,
        'anc': 'T',
        'der': 'A',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [0, 0, 1, 1, 1],
        'haplotypes_2': [1, 0, 1, 0, 1],
        'genotypes': [1, 0, 2, 1, 2],
        'arc_is_derived': False,  # has derived, no match
        'arc_der_count': -1,
        'arc_match': False,
    }

    opts.ancestral_bsg = None
    snp = read_vcf.process_vcf_line_to_genotypes(
        '1\t23\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n',
        opts)
    assert snp == {
        'sfs_target': 3,
        'sfs_reference': 1,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 23,
        'anc': 'A',
        'der': 'T',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [1, 1, 0, 0, 0],
        'haplotypes_2': [0, 1, 0, 1, 0],
        'genotypes': [1, 2, 0, 1, 0],
        'arc_der_count': -1,
        'arc_match': True,
    }


def test_vcf_to_genotypes_windowed(mocker):
    mocker.patch('read_vcf.read_vcf_header')
    mocker.patch('read_vcf.read_1kg_ind_pop_file')
    opts = AttributeDict()
    opts.vcf_has_illumina_chrnums = False
    opts.window_range = [7]
    opts.regions = mocker.MagicMock()
    opts.regions.in_region_one_based.side_effect = lambda _, pos: pos != 11
    opts.ancestral_bsg = mocker.MagicMock()
    bsg = {12: 'N', 13: 'T', 14: 'T', 15: 'C', 41: 'T', 42: 'T'}
    opts.ancestral_bsg.get_base_one_based.side_effect = (
        lambda _, pos: 'A' if pos not in bsg else bsg[pos])
    opts.reference_individuals_indexed_to_orig_file = [1, 2]
    opts.target_individuals_indexed_to_orig_file = [0, 4, 5]
    opts.exclude_individuals_indexed_to_orig_file = [3]
    opts.num_target = 3
    opts.num_reference = 2
    opts.archaic_vcf = None
    opts.debug = False
    opts.window_file = None

    vcf = StringIO(
        u'1\t6\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t9\t.\tAA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t10\t.\tA\tTT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t11\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t13\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n'
        '1\t14\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n'
        '1\t15\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n'
        '1\t16\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n'
        '1\t17\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t./.\n'
        '1\t19\t.\tA\tT\t.\tPASS\t.\tX:GT\t'
        'X:1|1\tX:1|1\tX:1|1\tX:0|0\tX:1|1\tX:1|1\n'
    )
    with pytest.raises(ValueError) as e:
        next(read_vcf.vcf_to_genotypes_windowed(vcf, 19, 100, None, opts))
    assert str(e.value) == 'Window length must be at least step size'

    snp_iter = read_vcf.vcf_to_genotypes_windowed(vcf, 19, 19, None, opts)

    expected = {
        'sfs_target': [1, 3, 3, 3],
        'sfs_reference': [0, 3, 1, 1],
        'target': [True, True, True, True],
        'reference': [False, True, True, True],
        'chrom': ['1', '1', '1', '1'],
        'pos': [7, 14, 15, 16],
        'anc': ['A', 'T', 'A', 'A'],
        'der': ['T', 'A', 'T', 'T'],
        'alt': ['T', 'T', 'T', 'T'],
        'ref': ['A', 'A', 'A', 'A'],
        'haplotypes_1': [[1, 0, 0, 0, 0], [0, 0, 1, 1, 1], [1, 1, 0, 0, 0],
                         [1, 1, 0, 0, 0]],
        'haplotypes_2': [[0, 0, 0, 0, 0], [1, 0, 1, 0, 1], [0, 1, 0, 1, 0],
                         [0, 1, 0, 1, 0]],
        'genotypes': [[1, 0, 0, 0, 0], [1, 0, 2, 1, 2], [1, 2, 0, 1, 0],
                      [1, 2, 0, 1, 0]],
    }

    chrom, start, end, snps = next(snp_iter)
    assert chrom == '1'
    assert start == 0
    assert end == 19
    for k in expected:
        actual = [s[k] for s in snps]
        assert actual == expected[k], k + " mismatch!"

    with pytest.raises(StopIteration):
        next(snp_iter)

    opts.archaic_vcf = mocker.MagicMock()
    opts.archaic_vcf.has_site.side_effect = lambda _, pos: pos != 42
    opts.archaic_vcf.has_derived.return_value = True
    opts.archaic_vcf.get_derived.return_value = 'T'
    opts.archaic_vcf.get_derived_count.return_value = -1

    vcf = StringIO(
        u'1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n'
        '1\t41\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n'
        '1\t42\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n'
    )
    snp_iter = read_vcf.vcf_to_genotypes_windowed(vcf, 19, 19, None, opts)
    chrom, start, end, snps = next(snp_iter)
    assert chrom == '1'
    assert start == 0
    assert end == 19
    assert snps == []

    chrom, start, end, snps = next(snp_iter)
    assert chrom == '1'
    assert start == 19
    assert end == 38
    assert snps == []

    expected = {
        'sfs_target': [3, 3, 3],
        'sfs_reference': [1, 3, 3],
        'target': [True, True, True],
        'reference': [True, True, True],
        'chrom': ['1', '1', '1'],
        'pos': [40, 41, 42],
        'anc': ['A', 'T', 'T'],
        'der': ['T', 'A', 'A'],
        'alt': ['T', 'T', 'T'],
        'ref': ['A', 'A', 'A'],
        'haplotypes_1': [[1, 1, 0, 0, 0], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1]],
        'haplotypes_2': [[0, 1, 0, 1, 0], [1, 0, 1, 0, 1], [1, 0, 1, 0, 1]],
        'genotypes': [[1, 2, 0, 1, 0], [1, 0, 2, 1, 2], [1, 0, 2, 1, 2]],
        'arc_is_derived': [True, False, False],
        'arc_der_count': [-1, -1, -1],
        'arc_match': [True, False, False],
    }

    chrom, start, end, snps = next(snp_iter)
    assert chrom == '1'
    assert start == 38
    assert end == 57
    for k in expected:
        actual = [s[k] for s in snps]
        assert actual == expected[k], k + " mismatch!"

    opts.archaic_vcf.add_site.assert_called_once_with(
        '1', 42, '1/1', 'T', 'A', None)

    with pytest.raises(StopIteration):
        next(snp_iter)

    opts.ancestral_bsg = None
    vcf = StringIO(
        u'1\t13\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n'
    )
    snp_iter = read_vcf.vcf_to_genotypes_windowed(vcf, 19, 19, None, opts)
    chrom, start, end, snps = next(snp_iter)
    assert chrom == '1'
    assert start == 0
    assert end == 19
    assert len(snps) == 1
    assert snps[0] == {
        'sfs_target': 3,
        'sfs_reference': 1,
        'target': True,
        'reference': True,
        'chrom': '1',
        'pos': 13,
        'anc': 'A',
        'der': 'T',
        'alt': 'T',
        'ref': 'A',
        'haplotypes_1': [1, 1, 0, 0, 0],
        'haplotypes_2': [0, 1, 0, 1, 0],
        'genotypes': [1, 2, 0, 1, 0],
        'arc_der_count': -1,
        'arc_match': True,
    }


def test_vcf_to_genotypes_windowed_window_file(mocker):
    mocker.patch('read_vcf.read_vcf_header')
    mocker.patch('read_vcf.read_1kg_ind_pop_file')
    opts = AttributeDict()
    opts.vcf_has_illumina_chrnums = False
    opts.window_range = [7]
    opts.regions = None
    opts.ancestral_bsg = None
    opts.reference_individuals_indexed_to_orig_file = [1, 2]
    opts.target_individuals_indexed_to_orig_file = [0, 4, 5]
    opts.exclude_individuals_indexed_to_orig_file = [3]
    opts.num_target = 3
    opts.num_reference = 2
    opts.archaic_vcf = None
    opts.debug = False
    opts.window_file = StringIO(
        u'1 2 3\n'
        '2 4 6\n'
        '2 3 4\n'
        '3 3 4\n'
    )
    opts.process_chromosome = 2
    snp_iter = read_vcf.vcf_to_genotypes_windowed(
        StringIO(
            u'2\t13\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n'
        ), 20, 20, None, opts, 20)  # 20's don't matter...
    assert next(snp_iter) == ('2', 4, 6, [])
    assert next(snp_iter) == ('2', 3, 4, [])
    with pytest.raises(StopIteration):
        next(snp_iter)
