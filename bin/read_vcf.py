import pandas as pd


def vcf_to_genotypes_windowed(vcf_file, winlen, winstep,
                              vcf_ind_pop_file, opts, start=0):

    # get sample ID order from the header / results are saved to opts
    read_vcf_header(vcf_file, opts)

    # process sample ID info / results are saved to opts
    read_1kg_ind_pop_file(vcf_ind_pop_file, opts)

    snps = []
    if opts.window_file is not None:
        chrom = None
        while chrom != opts.process_chromosome:
            tokens = opts.window_file.readline().split()
            print('reading chromosomes, looking for %s: %s' %
                  (opts.process_chromosome, chrom))
            (chrom, winstart, winend) = [int(i) for i in tokens[:3]]
        print('reading chromosomes, looking for %s: %s MATCH' %
              (opts.process_chromosome, chrom))
        print('WINFILE', chrom, winstart, winend)

    else:
        chrom = None
        winstart = start
        winend = winlen

    while True:
        for line in vcf_file:
            if opts.debug:
                print('line', line, 'line')

            snp_d = process_vcf_line_to_genotypes(line, opts)

            if snp_d is None:
                if opts.debug:
                    print()
                continue

            if snp_d['pos'] <= winstart:
                continue

            chrom = snp_d['chrom']

            # if it's past the window, then yield the current set of snps
            # then adjust the window, prune snps, and check to add it
            while snp_d['pos'] > winend:
                yield (chrom, winstart, winend, snps)
                if opts.window_file is not None:
                    tokens = opts.window_file.readline().split()
                    (winchrom, winstart, winend) = [int(i) for i in tokens[:3]]
                    if winchrom != opts.process_chromosome:
                        raise StopIteration
                else:
                    winstart += winstep
                    winend += winstep
                snps = [s for s in snps if s['pos'] > winstart]

            # now it's definitely in the window
            snps.append(snp_d)

        else:
            # end of file
            break

    yield (chrom, winstart, winend, snps)


def read_vcf_header(vcf_file, opts):

    line = vcf_file.readline()
    # header starts with just one #, comments start with ##
    while line.startswith('##'):
        line = vcf_file.readline()

    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT HG00096 HG00097
    if not line.startswith('#CHROM'):
        raise ValueError("bad VCF header?\n%s" % (line))

    opts.sample_index_in_original_file = {id: i for i, id in
                                          enumerate(line.split()[9:])}


def read_1kg_ind_pop_file(ind_pop_file, opts):

    """
    Sets up important options including:
     - the indices of reference and target individuals
     - ind to pop mappings, etc.
    Necessary to set indexes to the original VCF sample order:
        opts.reference_individuals_indexed_to_orig_file
        opts.target_individuals_indexed_to_orig_file
        opts.exclude_individuals_indexed_to_orig_file
    """

    # sample  pop     super_pop       gender
    # HG00096 GBR     EUR     male
    # HG00097 GBR     EUR     female
    # HG00099 GBR     EUR     female

    pop_data = pd.read_csv(ind_pop_file,
                           sep='\t',
                           usecols=[0, 1, 2],
                           names=['sample', 'pop', 'super_pop'],
                           skiprows=1
                           )
    opts.exclude_individuals = pop_data.loc[
        (pop_data['sample'].isin(opts.exclude_individuals)) |
        (pop_data['pop'].isin(opts.exclude_populations)) |
        (pop_data['super_pop'].isin(opts.exclude_populations))
         ]['sample']
    # these are used later
    refs = pop_data.loc[
        (pop_data['sample'].isin(opts.reference_individuals)) |
        (pop_data['pop'].isin(opts.reference_populations)) |
        (pop_data['super_pop'].isin(opts.reference_populations))
         ]
    targ = pop_data.loc[
        (pop_data['sample'].isin(opts.target_individuals)) |
        (pop_data['pop'].isin(opts.target_populations)) |
        (pop_data['super_pop'].isin(opts.target_populations))
         ]

    opts.reference_individuals = refs['sample']
    opts.target_individuals = targ['sample']

    opts.exclude_individuals_indexed_to_orig_file = [
        opts.sample_index_in_original_file[ind] for ind in
        opts.exclude_individuals]
    opts.reference_individuals_indexed_to_orig_file = [
        opts.sample_index_in_original_file[ind] for ind in
        opts.reference_individuals]
    opts.target_individuals_indexed_to_orig_file = [
        opts.sample_index_in_original_file[ind] for ind in
        opts.target_individuals]

    opts.num_reference = len(opts.reference_individuals)
    opts.num_target = len(opts.target_individuals)
    opts.num_samples = opts.num_reference + opts.num_target

    # append the target and reference indivs for indexing
    sample_ids = targ.append(refs, ignore_index=True)
    opts.target_indices = sample_ids.index[:opts.num_target]
    opts.reference_indices = sample_ids.index[opts.num_target:]

    opts.get_id_from_sample_index = sample_ids['sample'].tolist()
    opts.get_pop_from_sample_index = sample_ids['pop'].tolist()


def process_vcf_line_to_genotypes(line, opts):
    # ignore comments
    # header starts with just one #, comments start with # - distinguish?
    # we've already read the header (to get ind_ids)
    if line.startswith('#'):
        return None

    snp_d = {}

    split_line = line.split()
    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT HG00096 HG00097
    snp_d['chrom'] = split_line[0]
    chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + snp_d['chrom']
    snp_d['pos'] = int(split_line[1])
    snp_d['ref'] = split_line[3]
    snp_d['alt'] = split_line[4]
    snp_d['qual'] = split_line[5]
    snp_d['filt'] = split_line[6]
    snp_d['info'] = split_line[7]
    format = split_line[8]

    # drop snps before the specified range
    if opts.window_range is not None and snp_d['pos'] < opts.window_range[0]:
        return None

    # remove non-biallelic sites
    if len(snp_d['ref']) != 1 or len(snp_d['alt']) != 1:
        if opts.debug:
            print("dropping snp because it is not biallelic, or is an indel",
                  snp_d['chrom'], snp_d['pos'], snp_d['ref'], snp_d['alt'])
        return None

    # remove sites not in regions
    if opts.regions is not None and \
            not opts.regions.in_region_one_based(chrom, snp_d['pos']):
        if opts.debug:
            print("dropping snp because it is not in our regions definition:",
                  snp_d['chrom'], snp_d['pos'])
        return None

    if opts.regions is not None and opts.debug:
        print("KEEPING SNP because it IS in our regions definition:",
              snp_d['chrom'], snp_d['pos'])

    # remove sites not in ancestral genome
    use_derived = opts.ancestral_bsg is not None
    if use_derived and \
            'N' == opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos']):
        if opts.debug:
            print("dropping snp because it is not in the ancestral genome:",
                  snp_d['chrom'], snp_d['pos'],
                  opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos']))
        return None

    # flip for derived?
    flip_for_derived = (
        use_derived and snp_d['alt'].upper() ==
        opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos']))
    if use_derived and opts.debug:
        print("checking derived from ancestral genome:",
              snp_d['chrom'], snp_d['pos'],
              'ancestral =',
              opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos']),
              'ref/alt =', snp_d['ref'], snp_d['alt'],
              'flip_for_derived =', flip_for_derived,
              'DERIVED' if flip_for_derived else '')

    # get genotypes
    genotypes = split_line[9:]

    # remove individual genotypes that aren't in the target or reference
    # reordered ind ids in read_1kg_ind_pop_file
    r_gt = [genotypes[j] for j in
            opts.reference_individuals_indexed_to_orig_file]
    t_gt = [genotypes[j] for j in
            opts.target_individuals_indexed_to_orig_file]
    e_gt = [genotypes[j] for j in
            opts.exclude_individuals_indexed_to_orig_file]

    if flip_for_derived:
        flip_map = {'./.': '1|1', '.': '1|1', '0|0': '1|1',
                    '1|0': '0|1', '0|1': '1|0', '1|1': '0|0'}
        t_gt = [flip_map[gt] for gt in t_gt]
        r_gt = [flip_map[gt] for gt in r_gt]
        e_gt = [flip_map[gt] for gt in e_gt]

        # SWAP REF AND ALT, so ref is "ancestral"
        # not sure if this is the best thing,
        # because you might still want to know what the real ref/alt bases are
        (snp_d['anc'], snp_d['der']) = (snp_d['alt'], snp_d['ref'])

    else:
        (snp_d['der'], snp_d['anc']) = (snp_d['alt'], snp_d['ref'])

    e_gt = set(e_gt)

    if '0|1' in e_gt or '1|0' in e_gt or '1|1' in e_gt:
        if opts.debug:
            print("dropping snp because it's "
                  "present in the EXCLUDE population")
        return None

    genotypes = t_gt + r_gt
    t_gt = set(t_gt)
    r_gt = set(r_gt)

    snp_d['target'] = '0|1' in t_gt or '1|0' in t_gt or '1|1' in t_gt
    snp_d['reference'] = '0|1' in r_gt or '1|0' in r_gt or '1|1' in r_gt

    # finally, remove variants that aren't in the target or reference
    if not snp_d['target'] and not snp_d['reference']:
        if opts.debug:
            print("dropping snp because it's NOT present "
                  "in either the TARGET or REFERENCE")
        return None

    # also remove variants that are fixed in target and reference
    if '1|1' in t_gt and len(t_gt) == 1 and '1|1' in r_gt and len(r_gt) == 1:
        if opts.debug:
            print("dropping snp because it's FIXED "
                  "in both the TARGET and REFERENCE")
        return None

    # get the location of genotypes, if necessary
    if ':' in format:
        gt_loc = format.split(':')
        if 'GT' not in gt_loc:
            raise ValueError(
                "bad format?  GT not found\n%s\n%s" % (format, split_line))
        gt_loc = gt_loc.index('GT')

        # strip out just genotypes
        genotypes = [gt.split(':')[gt_loc] for gt in genotypes]

    if opts.debug:
        print('reading', snp_d['pos'], genotypes)

    hap_map1 = {'./.': 0, '.': 0, '0|0': 0, '1|0': 1, '0|1': 0, '1|1': 1}
    hap_map2 = {'./.': 0, '.': 0, '0|0': 0, '1|0': 0, '0|1': 1, '1|1': 1}

    snp_d['haplotypes_1'] = [hap_map1[gt] for gt in genotypes]
    snp_d['haplotypes_2'] = [hap_map2[gt] for gt in genotypes]
    snp_d['genotypes'] = [h1 + h2 for h1, h2 in zip(snp_d['haplotypes_1'],
                                                    snp_d['haplotypes_2'])]

    if opts.debug:
        print(snp_d['genotypes'])

    snp_d['sfs_target'] = sum(snp_d['genotypes'][:opts.num_target])
    snp_d['sfs_reference'] = sum(
        snp_d['genotypes'][opts.num_target:
                           opts.num_target + opts.num_reference])

    if opts.archaic_vcf is not None:
        if flip_for_derived and not opts.archaic_vcf.has_site(snp_d['chrom'],
                                                              snp_d['pos']):
            opts.archaic_vcf.add_site(snp_d['chrom'], snp_d['pos'],
                                      '1/1', snp_d['anc'], snp_d['der'], None)

        snp_d['arc_match'] = (opts.archaic_vcf.has_derived(snp_d['chrom'],
                                                           snp_d['pos']) and
                              snp_d['der'].upper() ==
                              opts.archaic_vcf.get_derived(snp_d['chrom'],
                                                           snp_d['pos']))

        snp_d['arc_der_count'] = opts.archaic_vcf.get_derived_count(
            snp_d['chrom'], snp_d['pos'])

        if opts.ancestral_bsg is not None:
            snp_d['arc_is_derived'] = snp_d['arc_match']

        if opts.debug:
            ch = snp_d['chrom']
            pos = snp_d['pos']
            print("neand match",
                  snp_d['pos'],
                  genotypes,
                  snp_d['alt'].upper(),
                  opts.archaic_vcf.get_ref_one_based(ch, pos),
                  opts.archaic_vcf.get_alts_one_based(ch, pos),
                  opts.archaic_vcf.get_bases_one_based(ch, pos),
                  opts.archaic_vcf.has_genotype_at_site(ch, pos),
                  'DEBUG',
                  'ARCHAIC_MATCH' if snp_d['arc_match'] else '',
                  'DERIVED' if flip_for_derived else '',
                  'JACKPOT' if (flip_for_derived and
                                snp_d['arc_match'] and
                                not opts.archaic_vcf.has_ref(ch, pos)) else '',
                  'MINI_JACKPOT' if (flip_for_derived and
                                     snp_d['arc_match']) else '',
                  'NEAND_IS_HET' if (
                      len(opts.archaic_vcf.get_bases_one_based(ch, pos)) == 2)
                  else 'NEAND_IS_HOM')

    if opts.debug:
        print()

    return snp_d


def vcf_to_genotypes(vcf_file):

    snps = []

    for line in vcf_file:

        snp_d = process_vcf_line_to_genotypes(line)
        if snp_d is None:
            continue
        snps.append(snp_d)

    return snps
