import read_archaic_vcf
from io import StringIO
import pytest


@pytest.fixture
def default_vcf(mocker):
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
    mocker.patch('read_archaic_vcf.open',
                 return_value=vcf_input)
    return read_archaic_vcf.vcf_class('test.vcf')


def test_init(mocker, capsys):
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
    mocker.patch('read_archaic_vcf.open',
                 return_value=vcf_input)
    vcf = read_archaic_vcf.vcf_class('test.vcf')
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


def test_init_execptions(mocker, capsys):
    # empty
    vcf_input = StringIO(
    )
    mocker.patch('read_archaic_vcf.open',
                 return_value=vcf_input)
    vcf = read_archaic_vcf.vcf_class('test.vcf')
    assert capsys.readouterr().err == (
        'Reading VCF file test.vcf..\n'
        ' with 0 lines.\n')

    # too few columns
    vcf_input = StringIO(
        u'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmsp_0\n'
        '1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\textra\n'
    )
    mocker.patch('read_archaic_vcf.open',
                 return_value=vcf_input)

    with pytest.raises(ValueError) as e:
        read_archaic_vcf.vcf_class('test.vcf')
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
    mocker.patch('read_archaic_vcf.open',
                 return_value=vcf_input)

    with pytest.raises(ValueError) as e:
        read_archaic_vcf.vcf_class('test.vcf')
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
    mocker.patch('read_archaic_vcf.open',
                 return_value=vcf_input)

    vcf = read_archaic_vcf.vcf_class('test.vcf', vcf_has_illumina_chrnums=True)

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
