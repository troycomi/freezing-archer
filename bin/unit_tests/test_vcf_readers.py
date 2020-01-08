from io import StringIO
import pytest
from vcf_readers import (ancestral_vcf, archaic_vcf)


@pytest.fixture
def default_archaic_vcf(mocker):
    vcf_input = StringIO(
        u'##fileformat=VCFv4.2\n'
        '##source=msprime 0.6.1\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##contig=<ID=1,length=10000000>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t164\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\n'
        '1\t202\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t219\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n'
        '1\t279\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\n'
        '1\t315\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '2\t315\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)
    return archaic_vcf('test.vcf')


@pytest.fixture
def empty_archaic_vcf(mocker):
    vcf_input = StringIO(
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)
    return archaic_vcf('test.vcf')


def test_archaic_init(mocker, capsys):
    vcf_input = StringIO(
        u'##fileformat=VCFv4.2\n'
        '##source=msprime 0.6.1\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##contig=<ID=1,length=10000000>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t164\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\n'
        '1\t202\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t219\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n'
        '1\t279\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\n'
        '1\t315\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '2\t315\t.\tA\tT\t.\tPASS\t.\tGT\t0|0andmorestuff\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)
    vcf = archaic_vcf('test.vcf')
    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 7 lines.\n')
    assert vcf.vcf == {
        '1': {
            7: (True, '0|0', 'A', 'T'),
            164: (False, '1|1', 'A', 'T'),
            202: (True, '0|0', 'A', 'T'),
            219: (True, '0|1', 'A', 'T'),
            279: (True, '1|0', 'A', 'T'),
            315: (True, '0|0', 'A', 'T'),
        },
        '2': {
            315: (True, '0|0', 'A', 'T'),
        }
    }


def test_archaic_init_execptions(mocker, capsys):
    # empty
    vcf_input = StringIO(
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)
    vcf = archaic_vcf('test.vcf')
    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 0 lines.\n')

    # too many columns
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\textra\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    with pytest.raises(ValueError) as e:
        archaic_vcf('test.vcf')
    assert str(e.value) == ('Too many columns in ARCHAIC VCF: test.vcf?\n'
                            'Expecting one individual (10 columns):\n'
                            '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\textra\n')

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        )

    # invalid position
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\ta\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    with pytest.raises(ValueError) as e:
        archaic_vcf('test.vcf')
    assert str(e.value) == (
        'Unable to parse position in ARCHAIC VCF: test.vcf\n'
        '1\ta\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n')

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        )

    # warn strip
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        'chr1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t8\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    vcf = archaic_vcf('test.vcf', vcf_has_illumina_chrnums=True)

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 2 lines.\n'
        " WARNING: REMOVED LEADING 'chr' FROM CHROMOSOME NAMES, "
        'TO MATCH ILLUMINA VCFS\n'
        )

    assert vcf.vcf == {
        '1': {
            7: (True, '0|0', 'A', 'T'),
            8: (True, '0|0', 'A', 'T'),
        }}

    # no warn strip and keep chr, illumina is False
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        'chr1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t8\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    vcf = archaic_vcf('test.vcf', vcf_has_illumina_chrnums=False)

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 2 lines.\n'
        )

    assert vcf.vcf == {
        'chr1': {
            7: (True, '0|0', 'A', 'T'),
        },
        '1': {
            8: (True, '0|0', 'A', 'T'),
        }
    }

    # duplicate position
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t7\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    with pytest.raises(ValueError) as e:
        archaic_vcf('test.vcf')
    assert str(e.value) == (
        'error - duplicate position in VCF file?\n'
        '1 7 G T\n'
        '1\t7\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\n'
    )

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        )

    # bi-allelic alt
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tTT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    with pytest.raises(ValueError) as e:
        archaic_vcf('test.vcf')
    assert str(e.value) == (
        'Error reading archaic VCF file:\n'
        'test.vcf\n'
        'len(alt) > 1: require only bi-allelic SNPs in archaic VCF\n'
        '1\t7\t.\tA\tTT\t.\tPASS\t.\tGT\t0|0\n'
    )

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        )


def test_archaic_vcf_init_with_regions(mocker, capsys):
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t8\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '2\t8\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t9\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)
    mock_region = mocker.MagicMock()
    mock_region.in_region_one_based.side_effect = lambda c, p: p != 8

    vcf = archaic_vcf('test.vcf', regions=mock_region)
    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 2 lines.\n'  # only reports passing
        )

    assert vcf.vcf == {
        '1': {
            7: (True, '0|0', 'A', 'T'),
            9: (True, '0|0', 'A', 'T'),
        }}


def test_archaic_vcf_init_with_bgs(mocker, capsys):
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t8\t.\tA\tA\t.\tPASS\t.\tGT\t0|0\n'
        '2\t8\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\n'
        '1\t9\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)
    mock_bsg = mocker.MagicMock()
    mock_bsg.get_base_one_based.side_effect = lambda c, p: 'T'

    vcf = archaic_vcf('test.vcf', ancestral_bsg=mock_bsg)
    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 4 lines.\n'  # only reports passing
        )

    # sites with T get messed with
    assert vcf.vcf == {
        '1': {
            7: (False, '1|1', 'T', 'A'),
            8: (True, '0|0', 'A', 'A'),
            9: (False, '1|1', 'T', 'A'),
        },
        '2': {
            8: (True, '0|0', 'A', 'G'),
        }
    }


def test_add_site(empty_archaic_vcf):
    vcf = empty_archaic_vcf

    assert vcf.vcf == {}

    # simple site
    vcf.add_site('1', 1, '1', 'r', 'a', 'c')
    assert vcf.vcf == {
        '1': {
            1: (False, '1', 'r', 'a')
        }
    }

    # existing site
    vcf.add_site('1', 1, '1', 'r2', 'a', 'c')
    assert vcf.vcf == {
        '1': {
            1: (False, '1', 'r2', 'a')
        }
    }

    # with 0 in gt
    vcf.add_site('1', 2, '10', 'r', 'a', 'c')
    assert vcf.vcf == {
        '1': {
            1: (False, '1', 'r2', 'a'),
            2: (True, '10', 'r', 'a'),
        }
    }

    # match ancestral
    vcf.add_site('2', 1, '10.', 'r', 'c', 'c')
    assert vcf.vcf == {
        '1': {
            1: (False, '1', 'r2', 'a'),
            2: (True, '10', 'r', 'a'),
        },
        '2': {
            1: (True, '010', 'c', 'r'),
        }
    }


def test_archaic_has_derived(default_archaic_vcf):
    vcf = default_archaic_vcf
    assert not vcf.has_derived('not exist', 0)
    assert not vcf.has_derived('1', 7)  # 0|0
    assert vcf.has_derived('1', 219)  # 0|1
    assert vcf.has_derived('1', 279)  # 1|0
    assert vcf.has_derived('1', 164)  # 1|1


def test_archaic_get_derived(default_archaic_vcf):
    vcf = default_archaic_vcf
    assert vcf.get_derived('not exist', 0) == 'N'
    assert vcf.get_derived('1', 219) == 'T'


def test_archaic_get_derived_count(default_archaic_vcf):
    vcf = default_archaic_vcf
    answers = {
        '1': {
            219: 1,
            164: 2,
            7: 0,
            202: 0,
            279: 1,
            315: 0,
        },
        '2': {315: 0}}
    for chrom in vcf.vcf:
        for pos in vcf.vcf[chrom]:
            assert vcf.get_derived_count(chrom, pos) == answers[chrom][pos]

    assert vcf.get_derived_count('not exist', 0) == 0


def test_archaic_has_site(default_archaic_vcf):
    vcf = default_archaic_vcf
    assert vcf.has_site('1', 7)
    assert not vcf.has_site('1', 8)
    assert not vcf.has_site('no', 7)


def test_archaic_get_derived_sites(default_archaic_vcf):
    vcf = default_archaic_vcf
    assert vcf.get_derived_sites('no', 0, 10) == []
    assert vcf.get_derived_sites('1', 0, 6) == []
    assert vcf.get_derived_sites('1', 7, 8) == []
    assert vcf.get_derived_sites('1', 219, 220) == [219]
    assert vcf.get_derived_sites('1', 220, 221) == []
    assert vcf.get_derived_sites('1', 200, 300) == [219, 279]


def test_archaic_get_derived_sites_with_der_count(default_archaic_vcf):
    vcf = default_archaic_vcf
    assert vcf.get_derived_sites_with_der_count('no', 0, 10) == []
    assert vcf.get_derived_sites_with_der_count('1', 0, 6) == []
    assert vcf.get_derived_sites_with_der_count('1', 7, 8) == []
    assert vcf.get_derived_sites_with_der_count('1', 219, 220) == [(219, 1)]
    assert vcf.get_derived_sites_with_der_count('1', 220, 221) == []
    assert vcf.get_derived_sites_with_der_count('1', 200, 300) == [(219, 1),
                                                                   (279, 1)]
    assert vcf.get_derived_sites_with_der_count('1', 100, 200) == [(164, 2)]
    assert vcf.get_derived_sites_with_der_count('2', 100, 200) == []


def test_archaic_generate_derived_cache(default_archaic_vcf):
    vcf = default_archaic_vcf
    assert vcf.derived_cache is None
    vcf.generate_derived_cache()
    assert vcf.derived_cache == {
        '1': [
            (164, 2),
            (219, 1),
            (279, 1)
        ],
        '2': [
        ]
    }


def test_ancestral_normal(mocker):
    vcf_input = StringIO(
        u'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t164\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\n'
    )

    vcf = ancestral_vcf(vcf_input)

    assert vcf.vcf.to_dict() == {
        'ref': {
            (1, 164): 'G',
            (1, 7): 'A',
        }
    }

    assert vcf.get_base_one_based(1, 7) == 'A'
    assert vcf.get_base_one_based(1, 164) == 'G'
    assert vcf.get_base_one_based(1, 8) == 'N'


def test_ancestral_chr(mocker):
    vcf_input = StringIO(
        u'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        'chr1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        'chr1\t164\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\n'
    )

    vcf = ancestral_vcf(vcf_input)

    assert vcf.vcf.to_dict() == {
        'ref': {
            ('chr1', 164): 'G',
            ('chr1', 7): 'A',
        }
    }

    assert vcf.get_base_one_based('chr1', 7) == 'A'
    assert vcf.get_base_one_based('chr1', 164) == 'G'
    assert vcf.get_base_one_based('chr1', 8) == 'N'


def test_ancestral_duplicate(mocker):
    vcf_input = StringIO(
        u'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
        '1\t7\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\n'
        '2\t7\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\n'
        '2\t8\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\n'
    )

    with pytest.raises(ValueError) as e:
        ancestral_vcf(vcf_input)

    lines = str(e.value).split('\n')
    assert len(lines) == 5
    assert lines[0] == "error - duplicate position in VCF file?"
    assert lines[1].strip() == "ref"
    assert lines[2].split() == 'chromosome position'.split()
    assert lines[3].split() == '1 7 A'.split()
    assert lines[4].split() == '7 G'.split()
