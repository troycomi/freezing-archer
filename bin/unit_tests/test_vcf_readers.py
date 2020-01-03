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
        '2\t315\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
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

    # too few columns
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
    )
    mocker.patch('vcf_readers.open',
                 return_value=vcf_input)

    vcf = archaic_vcf('test.vcf', vcf_has_illumina_chrnums=True)

    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 1 lines.\n'
        " WARNING: REMOVED LEADING 'chr' FROM CHROMOSOME NAMES, "
        'TO MATCH ILLUMINA VCFS\n'
        )

    assert vcf.vcf == {
        '1': {
            7: (True, '0|0', 'A', 'T'),
        }}


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
