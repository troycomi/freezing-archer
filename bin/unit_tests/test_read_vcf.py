import read_vcf
from io import StringIO


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

    assert opts.exclude_individuals == ['HG00099']
    assert opts.exclude_individuals_indexed_to_orig_file == [12]
    # target + reference list
    assert opts.get_id_from_sample_index(0) == 'HG00097'
    assert opts.get_id_from_sample_index(1) == 'HG00096'
    assert opts.get_pop_from_sample_index(0) == 'GBR2'
    assert opts.get_pop_from_sample_index(1) == 'GBR1'
    assert opts.num_reference == 1
    assert opts.num_samples == 2
    assert opts.num_target == 1
    assert opts.reference_indices == [1]
    assert opts.reference_individuals == ['HG00096']
    assert opts.reference_individuals_indexed_to_orig_file == [10]
    assert opts.target_indices == [0]
    assert opts.target_individuals == ['HG00097']
    assert opts.target_individuals_indexed_to_orig_file == [11]

    for k in sorted(opts.keys()):
        print('{k} -> {v}'.format(k=k, v=opts[k]))
