import BaseLookup
import pytest
from io import StringIO


def test_init():
    bl = BaseLookup.BaseLookup()
    assert bl.basedir == './'

    bl = BaseLookup.BaseLookup(basedir='.')
    assert bl.basedir == './'
    assert bl.seq_files == {}
    assert bl.seq_files_header_len == {}
    assert bl.bases_per_line == {}
    assert bl.chr_len_hash['hg19'] == BaseLookup.BaseLookup.chr_lens
    assert bl.chr_len_hash['hg18'] == BaseLookup.BaseLookup.chr_lens_36
    assert bl.chr_len_hash['tair10'] == BaseLookup.BaseLookup.chr_lens_tair10

    assert not hasattr(bl, 'genome_len')
    assert bl.chr_offset == BaseLookup.BaseLookup.chr_offset
    assert bl.chr_lens == BaseLookup.BaseLookup.chr_lens
    assert bl.chrs == BaseLookup.BaseLookup.chrs

    bl = BaseLookup.BaseLookup(ref_version='b36')
    assert bl.genome_len == BaseLookup.BaseLookup.chr_lens_36['genome']
    assert bl.chr_offset == BaseLookup.BaseLookup.chr_offset_36

    bl = BaseLookup.BaseLookup(ref_version='tair10')
    assert bl.genome_len == BaseLookup.BaseLookup.chr_lens_tair10['genome']
    assert bl.chr_offset == BaseLookup.BaseLookup.chr_offset_tair10
    assert bl.chr_lens == BaseLookup.BaseLookup.chr_lens_tair10
    assert bl.chrs == BaseLookup.BaseLookup.chrs_tair10

    with pytest.raises(ValueError) as e:
        BaseLookup.BaseLookup(ref_version='test')
    assert ('Invalid reference version test.  '
            'Supported values are b37, b36, or tair10') in e.value


def test_constants():
    assert BaseLookup.BaseLookup.chr_nums == {
        'sim': 1,
        'chr1': 1,
        'chr2': 2,
        'chr3': 3,
        'chr4': 4,
        'chr5': 5,
        'chr6': 6,
        'chr7': 7,
        'chr8': 8,
        'chr9': 9,
        'chr10': 10,
        'chr11': 11,
        'chr12': 12,
        'chr13': 13,
        'chr14': 14,
        'chr15': 15,
        'chr16': 16,
        'chr17': 17,
        'chr18': 18,
        'chr19': 19,
        'chr20': 20,
        'chr21': 21,
        'chr22': 22,
        'chrX': 23,
        'chrY': 24,
        'chrM': 25}

    assert BaseLookup.BaseLookup.chrs == (
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
        'chrX', 'chrY', 'chrM')

    assert BaseLookup.BaseLookup.chr_lens['genome'] == 3095693981
    assert BaseLookup.BaseLookup.chr_lens_36['genome'] == 3080436051
    assert BaseLookup.BaseLookup.chr_lens_tair10['genome'] == 119667750


def test_getBase_exceptions(mocker):
    bl = BaseLookup.BaseLookup()
    with pytest.raises(ValueError) as e:
        bl.getBase('test', -1, -1)
    assert 'Invalid species "test". Options are ' in str(e.value)

    with pytest.raises(ValueError) as e:
        bl.getBase('hg18', 1, -1)
    assert 'Invalid chromosome "1". Options are ' in str(e.value)

    with pytest.raises(ValueError) as e:
        bl.getBase('hg18', 'chr22', -1)
    assert ('Base location "-1" out of range. '
            'Valid range is [1, 49691432]') in str(e.value)

    with pytest.raises(ValueError) as e:
        bl.getBase('hg18', 'chr22', 1e8)
    assert ('Base location "100000000" out of range. '
            'Valid range is [1, 49691432]') in str(e.value)

    mock_fa = StringIO(
        u'HEADER\n'
        u'this is the first line\n'
        u'abcdefghijklmnopqrstuv\n'
        u'0123456789012345678901\n'
        u'abc\ndef\n'
    )
    mock_open = mocker.patch('BaseLookup.open', return_value=mock_fa)
    result = bl.getBase('hg18', 'chr22', 10)
    mock_open.assert_called_once_with('./hg18/chr22.fa', 'r')
    bl.seq_files_header_len['hg18']['chr22'] == 7
    bl.bases_per_line['hg18']['chr22'] == 22
    assert result == 'h'

    result = bl.getBase('hg18', 'chr22', 1)
    assert result == 't'

    result = bl.getBase('hg18', 'chr22', 23)
    assert result == 'a'

    result = bl.getBase('hg18', 'chr22', 25)
    assert result == 'c'

    result = bl.getBase('hg18', 'chr22', 44)
    assert result == 'v'

    result = bl.getBase('hg18', 'chr22', 50)
    assert result == '5'

    with pytest.raises(ValueError) as e:
        bl.getBase('hg18', 'chr22', 70)
    assert (
        'bad base pos (base returned is newline): chr22 70\n'
        'abcNdef\n'
        '   ^'
    ) in str(e.value)


def test_baselookup_context(mocker):
    mock_file = mocker.MagicMock()
    mock_file.read.return_value = '1'
    # to set header and bp per line
    mock_file.readline.return_value = '1234567890'
    mock_open = mocker.patch('BaseLookup.open', return_value=mock_file)
    with BaseLookup.BaseLookup() as bl:
        for s in 'hg18 hg19'.split():
            for chrom in 'chr1 chr2'.split():
                for site in (10, 20):
                    assert bl.getBase(s, chrom, site) == '1'

    # all files opened
    assert mock_open.call_count == 4
    assert mock_open.call_args_list == [
        mocker.call('./hg18/chr1.fa', 'r'),
        mocker.call('./hg18/chr2.fa', 'r'),
        mocker.call('./hg19/chr1.fa', 'r'),
        mocker.call('./hg19/chr2.fa', 'r'),
    ]

    # and closed
    assert mock_file.close.call_count == 4
