from custom_argparse import ancestral_vcf
from io import StringIO
import pytest


# these are all tests for ancetral_vcf, as the organization is off
def test_normal(mocker):
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


def test_chr(mocker):
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


def test_duplicate(mocker):
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
