from __future__ import division
import sys
import numpy as np
from collections import Counter
import arc_match_pval_tables


# TODO remove these once calc_s_dists is tested
def calc_geno_dist(gt1, gt2):
    return np.sum(np.abs(np.asarray(gt1) - np.asarray(gt2)))


def calc_s(gt1, gt2, position1, position2,
           match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000):

    dist = abs(position1 - position2)
    if dist < 10:
        return -sys.maxint

    gd = calc_geno_dist(gt1, gt2)
    if gd == 0:
        return match_bonus + dist

    elif gd <= max_mismatch:
        return mismatch_penalty

    return -sys.maxint


def calc_s_dists(genotypes, positions,
                 match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000):
    # all pairwise dists[i, j] from i to j
    pos = np.asarray(positions)
    bp = pos[None, :] - pos[:, None]
    gen = np.asarray(genotypes)
    return np.where(bp < 10,
                    np.NINF,
                    np.where(gen[None, :] == gen[:, None],
                             match_bonus + bp,
                             mismatch_penalty))


def calc_s_star(genotypes, positions,
                match_bonus, max_mismatch, mismatch_penalty):
    '''
    loop through all snps, calculating s_star given that:
        you have a set of snps that ends in k
        you are now trying to figure out the set of snps
            before k that maximizes the score with k
    try all possible j (i.e., sets that end in j), save the best one
    consider both the previous best set with j,
        and that you start over with (j,k)
    '''

    nsnps = len(genotypes)
    s_star_scores = np.zeros(nsnps)
    s_star_snps = [[]] * nsnps
    dists = calc_s_dists(genotypes, positions,
                         match_bonus, max_mismatch, mismatch_penalty)

    # TODO test adding numba here
    for k in range(nsnps):
        max_score = -sys.maxint
        max_snps = []
        for j in range(k):
            # just the set (j,k)
            score1 = dists[j, k]
            # the previous set that ends in j, plus k
            score2 = s_star_scores[j] + score1

            if max_score < score2:
                max_score = score2
                max_snps = s_star_snps[j] + [k]

            if max_score < score1:
                max_score = score1
                max_snps = [j, k]

        s_star_scores[k] = max_score
        s_star_snps[k] = max_snps

    max_ind = s_star_scores.argmax()

    return (int(s_star_scores[max_ind]), s_star_snps[max_ind])


def calc_match_pval_from_genome_tables(
        chrom, winstart, winend, ind_snps, opts):

    match_pval_debug = False

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    my_mapped_bases = (opts.window_length
                       if opts.regions is None
                       else opts.regions.amount_in_region(
                           full_chrom, winstart, winend))
    my_mapped_bases_bin = my_mapped_bases // 1000 * 1000

    # get non-ref snps
    neand_pos = [p for p in range(winstart, winend)
                 if opts.archaic_vcf.has_var_at_site(chrom, p) and
                 not opts.archaic_vcf.has_ref(chrom, p) and
                 (opts.regions is None
                  or opts.regions.in_region_one_based(full_chrom, p))]

    # get the number of sites in the region (for this ind, and for neand)
    ind_pos = set([s['pos'] for s in ind_snps])
    ind_pos2 = ind_pos.union(set(neand_pos))
    my_n_sites = len(ind_pos2)

    # count the number of matches to neanderthal
    n_match = sum([s['arc_match'] for s in ind_snps])
    test_ratio = n_match / my_n_sites

    # get the average freq
    my_sfs = int(round(sum([s['sfs_target'] for s in ind_snps])
                       / len(ind_snps)))

    # get the std freq
    my_std_dev = int(round(np.std([s['sfs_target'] for s in ind_snps])))

    # get the distribution of "null" matches

    # cycle through files
    my_counter = Counter()
    my_null = []
    for idx_factor in range(5):
        for tbl_n, tbl in enumerate(opts.match_pval_table):
            if match_pval_debug:
                print('table', tbl_n)
            for b in [x for x in range(my_mapped_bases_bin - 2000*idx_factor,
                                       my_mapped_bases_bin + 2000*idx_factor+1,
                                       1000)
                      if x in tbl.keys()]:
                if match_pval_debug:
                    print(b)
                for p in [x for x in range(my_n_sites - 3*idx_factor,
                                           my_n_sites + 3*idx_factor+1)
                          if x in tbl[b].keys()]:
                    if match_pval_debug:
                        print(b, p)
                    for sfs in [x for x in range(my_sfs - 1*idx_factor,
                                                 my_sfs + 1*idx_factor+1)
                                if x in tbl[b][p].keys()]:
                        if match_pval_debug:
                            print(b, p, sfs)
                        for sd in [x for x in range(my_std_dev - 5*idx_factor,
                                                    my_std_dev + 5*idx_factor+1
                                                    )
                                   if x in tbl[b][p][sfs].keys()]:
                            if match_pval_debug:
                                print(b, p, sfs, sd, tbl[b][p][sfs][sd])
                            for nm, c in tbl[b][p][sfs][sd].items():
                                my_null += [nm/p] * c
                                if match_pval_debug:
                                    print(b, p, sfs, sd, nm, c, my_null)
                            my_counter += tbl[b][p][sfs][sd]

            if match_pval_debug:
                print("current counter:", my_counter)
            if match_pval_debug:
                print("current null:", len(my_null), my_null)

        if match_pval_debug:
            print("finished all files with idx_factor=%d.." % idx_factor)
        if len(my_null) > 100:
            break

    if match_pval_debug:
        print("FINISHED GENERATING NULL")
    for n in sorted(my_null):
        if match_pval_debug:
            print(n)
        pass

    null_hits = sum([n >= test_ratio for n in my_null])
    total_hits = len(my_null)

    if total_hits == 0:
        pval = 'NA'
    else:
        pval = null_hits / total_hits

    return (my_mapped_bases, my_n_sites, my_sfs, my_std_dev, n_match,
            len(ind_snps), total_hits, null_hits, pval, idx_factor)


def initialize_analysis(opts):
    setattr(opts, 'pval_repeat_lookup', {})


def calc_table_match_pval(chrom, snps, opts, window_neand_pos,
                          ind, hap1_range, hap2_range):
    # get the list of snps on each haplotype
    ind_hap1 = [s for s in snps
                if s['haplotypes_1'][ind] > 0
                and hap1_range[0] <= s['pos'] <= hap1_range[1]]
    ind_hap2 = [s for s in snps
                if s['haplotypes_2'][ind] > 0
                and hap2_range[0] <= s['pos'] <= hap2_range[1]]

    ret_p = ['NA', 'NA']
    ret_m = ['NA', 'NA']
    ret_n = ['NA', 'NA']
    ret_f = ['NA', 'NA']
    ret_q = [['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'],
             ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']]

    for hapnum, ind_snps in enumerate((ind_hap1, ind_hap2)):
        if len(ind_snps) < 2:
            continue

        # get neand snps for this region
        neand_pos = [p for p in window_neand_pos
                     if ind_snps[0]['pos'] <= p <= ind_snps[-1]['pos']]

        (test_mapped_bases_bin, test_len,
         test_mh_sites, test_tot_sites,
         orig_sfs, _, test_match) = arc_match_pval_tables.get_match_stats(
             chrom, ind_snps, neand_pos, opts, pop='sfs_target')

        # HACK!  should just save sfs as an int 1-100?
        # map test_sfs to be out of 216 (108 YRI inds in table)
        test_sfs = int(orig_sfs / opts.num_target * opts.ptables_ninds_for_sfs)

        if (test_len < 10000 or
                test_mapped_bases_bin < 5000 or
                test_tot_sites < 8):
            continue

        eps_len = opts.len_eps
        eps_mapped_bases_bin = opts.mapped_eps
        eps_mh_sites = 1
        eps_sfs = 5

        factor = 1
        tot_count = 0

        match_pct = test_match / test_tot_sites

        hap_len_range = (test_len - factor * eps_len,
                         test_len + factor * eps_len)
        mapped_bases_bin = (
            test_mapped_bases_bin - factor * eps_mapped_bases_bin,
            test_mapped_bases_bin + factor * eps_mapped_bases_bin)
        mh_sites = (test_mh_sites - factor * eps_mh_sites,
                    test_mh_sites + factor * eps_mh_sites)
        sfs = (test_sfs - factor * eps_sfs, test_sfs + factor * eps_sfs)
        query_args = {
            'len': '(%f <= hap_len) & (hap_len <= %f)' % hap_len_range,
            'mapped': '(%f <= mapped_bases_bin) & (mapped_bases_bin <= %f)'
            % mapped_bases_bin,
            'mh': '(%f <= mh_sites) & (mh_sites <= %f)' % mh_sites ,
            'sfs': '(%f <= sfs) & (sfs <= %f)' % sfs}

        query = ' & '.join(query_args[x] for x in opts.table_query)

        if query in opts.pval_repeat_lookup:

            if opts.debug_pval_lookup:
                print('GETTING PVAL FROM LOOKUP',
                      len(opts.pval_repeat_lookup), query)

            ret_p[hapnum] = opts.pval_repeat_lookup[query][0]
            ret_m[hapnum] = opts.pval_repeat_lookup[query][1]
            ret_n[hapnum] = opts.pval_repeat_lookup[query][2]
            ret_f[hapnum] = opts.pval_repeat_lookup[query][3]
            ret_q[hapnum] = opts.pval_repeat_lookup[query][4]
            continue

        if test_match > 0:
            a = [x[:] for x in opts.match_pval_table.root.table.where(query)]
        else:
            a = []

        tot_count = sum(x[0] for x in a)

        match_pct_by_rows = [x[0] if x[7]/x[4] >= match_pct else 0 for x in a]
        pval = (sum(match_pct_by_rows) + 1) / (tot_count + 1)

        ret_p[hapnum] = pval
        ret_m[hapnum] = match_pct
        ret_n[hapnum] = tot_count
        ret_f[hapnum] = factor / 1.1
        ret_q[hapnum] = (test_len, test_mapped_bases_bin, test_mh_sites,
                         test_sfs, orig_sfs, test_match, test_tot_sites)

        if opts.debug_pval_lookup:
            print('SETTING PVAL FOR LOOKUP',
                  len(opts.pval_repeat_lookup), query)
        opts.pval_repeat_lookup[query] = (ret_p[hapnum], ret_m[hapnum],
                                          ret_n[hapnum], ret_f[hapnum],
                                          ret_q[hapnum])

    return (ret_p, ret_m, ret_n, ret_f, ret_q)


def calc_local_match_pval(chrom, winstart, winend, snps, opts,
                          ind_indices=None, subset_start=None,
                          subset_end=None):

    if ind_indices is None:
        ind_indices = opts.target_indices

    if isinstance(ind_indices, int):
        ind_indices = [ind_indices]

    if subset_start is None:
        subset_start = winstart
    # subset_end is *inclusive*, unlike winend
    if subset_end is None:
        subset_end = winend-1

    # get non-ref snps for this region
    neand_pos = opts.archaic_vcf.get_derived_sites(chrom, winstart, winend)

    # save two spots for each individual (one for each haplotype)
    # the code below assumes that the ref inds are first, then target
    match_pct = [None, None] * opts.num_samples

    # just look at sites in the window of interest
    subset_snps = [s for s in snps
                   if s['pos'] >= subset_start
                   and s['pos'] <= subset_end]

    # loop through all individuals (target and ref) and calculate match pct
    for ind in ind_indices + opts.reference_indices:

        # get haplotype snps (two haplotypes per person, natch)
        ind_hap1 = [s for s in subset_snps if s['haplotypes_1'][ind] > 0]
        ind_hap2 = [s for s in subset_snps if s['haplotypes_2'][ind] > 0]

        # get the number of sites to consider for this haplotype
        ind_pos1 = set([s['pos'] for s in ind_hap1])
        ind_pos1_n = ind_pos1.union(set(neand_pos))
        my_n_sites1 = len(ind_pos1_n)

        ind_pos2 = set([s['pos'] for s in ind_hap2])
        ind_pos2_n = ind_pos2.union(set(neand_pos))
        my_n_sites2 = len(ind_pos2_n)

        # count the number of matches to neanderthal
        n_match1 = sum([s['arc_match'] for s in ind_hap1])
        n_match2 = sum([s['arc_match'] for s in ind_hap2])

        # get the pct of sites that match neand
        test_ratio1 = None if my_n_sites1 == 0 else n_match1 / my_n_sites1
        test_ratio2 = None if my_n_sites2 == 0 else n_match2 / my_n_sites2

        # save those values
        match_pct[ind*2] = test_ratio1
        match_pct[ind*2+1] = test_ratio2

        if opts.debug:
            print(('debug_match_pct_fn IND=%d HAP=1 PCT=%f '
                   'MATCH_NEAND=%s NEAND_POS=%s DENOM=%s') % (
                       ind, test_ratio1, str([s['pos'] for s in ind_hap1
                                              if s['arc_match']]),
                       str(neand_pos), str(my_n_sites1)))
        if opts.debug:
            print(('debug_match_pct_fn IND=%d HAP=2 PCT=%f '
                   'MATCH_NEAND=%s NEAND_POS=%s DENOM=%s') % (
                       ind, test_ratio2, str([s['pos'] for s in ind_hap2
                                              if s['arc_match']]),
                       str(neand_pos), str(my_n_sites2)))

        if opts.debug:
            print("test_ratio1", winstart, winend, ind,
                  opts.get_pop_from_sample_index(ind), n_match1,
                  len(neand_pos), len(ind_pos1), my_n_sites1, test_ratio1)
        if opts.debug:
            print("test_ratio2", winstart, winend, ind,
                  opts.get_pop_from_sample_index(ind), n_match2,
                  len(neand_pos), len(ind_pos2), my_n_sites2, test_ratio2)

    # THIS SHOULD BE CONVERTED TO WORK WITH OPTS.TARGET_INDICES
    ref_nones = sum(p is None for p in match_pct[opts.num_target*2:])

    # make a longer pval list than you need, so you can use standard indices
    pvals = [None, None] * opts.num_samples
    # one individual's two haplotypes aren't paired with each other
    # but that shouldn't matter - we're always looking at them in aggregate
    ref_match_pct = ([match_pct[i*2] for i in opts.reference_indices] +
                     [match_pct[i*2+1] for i in opts.reference_indices])
    for ind in ind_indices:
        ps_1 = sum(match_pct[ind*2] <= p for p in ref_match_pct)
        ps_2 = sum(match_pct[ind*2+1] <= p for p in ref_match_pct)
        pvals[ind*2] = (ps_1 if ps_1 > 0 else 1) / opts.num_reference / 2
        pvals[ind*2+1] = (ps_2 if ps_2 > 0 else 1) / opts.num_reference / 2
        if opts.debug:
            print("pval1", winstart, winend, ind,
                  opts.get_pop_from_sample_index(ind),
                  match_pct[ind*2], pvals[ind*2])
        if opts.debug:
            print("pval2", winstart, winend, ind,
                  opts.get_pop_from_sample_index(ind),
                  match_pct[ind*2+1], pvals[ind*2+1])

    return (pvals, match_pct, ref_nones)


def run_window_analysis(chrom, winstart, winend, snps, opts):
    if len(snps) <= 2:
        return

    if opts.archaic_vcf is None or opts.no_pvalues:
        match_pvals2 = ['NA', 'NA'] * opts.num_samples
        match_pct2 = ['NA', 'NA'] * opts.num_samples
        ref_nones = 0
    else:
        (match_pvals2, match_pct2, ref_nones) = calc_local_match_pval(
            chrom, winstart, winend, snps, opts)
        pass

    window_neand_pos = opts.archaic_vcf.get_derived_sites(
        chrom, winstart, winend)
    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    my_mapped_bases = opts.regions.amount_in_region(
        full_chrom, winstart, winend) if opts.regions is not None \
        else opts.window_length

    # JUST LOOP OVER TARGET INDS
    for ind in opts.target_indices:

        ind1_snps_all = [s for s in snps
                         if s['genotypes'][ind] > 0]
        ind1_snps = [s for s in ind1_snps_all
                     if s['target']
                     and not s['reference']]

        # used as a count of the number of snvs in a region
        # for matching to the null model
        ind1_or_ref_snps = [s for s in snps
                            if s['genotypes'][ind] > 0 or s['reference']]

        ind1_pos = [s['pos'] for s in ind1_snps]

        # masks for which snps are on each haplotypes
        ind1_hap1 = [s['haplotypes_1'][ind] > 0 for s in ind1_snps]
        ind1_hap2 = [s['haplotypes_2'][ind] > 0 for s in ind1_snps]

        if len(ind1_snps) <= 2:
            (s_star, s_star_snps) = (0, [])

        else:
            # s_star is the score and s_star_snps are the indices
            # (wrt ind1_snps) that are selected by s_star
            (s_star, s_star_snps) = calc_s_star(
                [s['genotypes'][ind] for s in ind1_snps],
                ind1_pos,
                match_bonus=opts.s_star_match_bonus,
                max_mismatch=opts.s_star_max_mismatch,
                mismatch_penalty=opts.s_star_mismatch_penalty)

        # find the best haplotype
        n_haps1 = sum(ind1_hap1[i] for i in s_star_snps)
        n_haps2 = sum(ind1_hap2[i] for i in s_star_snps)
        all_haps = (ind1_hap1[i] + ind1_hap2[i]*2 for i in s_star_snps)

        # get S* haplotype range
        hap1_pos = [ind1_snps[i]['pos'] for i in s_star_snps if ind1_hap1[i]]
        hap1_range = (min(hap1_pos), max(hap1_pos)) \
            if len(hap1_pos) > 1 else (0, 0)
        hap2_pos = [ind1_snps[i]['pos'] for i in s_star_snps if ind1_hap2[i]]
        hap2_range = (min(hap2_pos), max(hap2_pos)) \
            if len(hap2_pos) > 1 else (0, 0)

        # currently depreciated..
        if opts.match_pval_table is not None:
            (match_table_pvals, match_table_pcts,
             match_table_Ns, match_table_factors,
             match_table_query) = calc_table_match_pval(
                 chrom, snps, opts, window_neand_pos, ind,
                 hap1_range, hap2_range)
        else:
            match_table_pvals = ('NA', 'NA')
            match_table_pcts = ('NA', 'NA')
            match_table_Ns = ('NA', 'NA')
            match_table_factors = ('NA', 'NA')
            match_table_query = [['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'],
                                 ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']]

        if len(s_star_snps) < 2 or opts.archaic_vcf is None or opts.no_pvalues:
            s_start = 'NA'
            s_stop = 'NA'
            match_pvals3 = ['NA', 'NA'] * opts.num_samples
            match_pct3 = ['NA', 'NA'] * opts.num_samples

        else:
            # get the start and stop of the S* region
            s_start = ind1_snps[min(s_star_snps)]['pos']
            s_stop = ind1_snps[max(s_star_snps)]['pos']
            (match_pvals3, match_pct3, ref_nones3) = calc_local_match_pval(
                chrom, winstart, winend, snps, opts, ind,
                subset_start=s_start, subset_end=s_stop)

        ms_fields = []
        ms_fields_labels = []
        if opts.vcf_is_ms_file:
            # count sites that are intr / S* / both
            # N: sites that are a) on that hap at all
            ms_N_sites_hap1 = [s for s in ind1_snps_all
                               if s['haplotypes_1'][ind] > 0]
            ms_N_sites_hap2 = [s for s in ind1_snps_all
                               if s['haplotypes_2'][ind] > 0]

            # S*: sites that are
            # a) in the sstar hap range, and
            # b) on that hap at all
            ms_ss_sites_hap1 = [s for s in ms_N_sites_hap1
                                if (hap1_range[0] <= s['pos']
                                    <= hap1_range[1])]
            ms_ss_sites_hap2 = [s for s in ms_N_sites_hap2
                                if (hap2_range[0] <= s['pos']
                                    <= hap2_range[1])]

            # intr and S*: sites that are
            # a) in the sstar hap range,
            # b) introgressed on that hap, and
            # c) on that hap at all (this may be redundant!)
            ms_intr_ss_sites_hap1 = [s for s in ms_ss_sites_hap1
                                     if s['haplotypes_1_intr'][ind]]
            ms_intr_ss_sites_hap2 = [s for s in ms_ss_sites_hap2
                                     if s['haplotypes_2_intr'][ind]]

            # intr: sites that are
            # a) introgressed on that hap, and
            # b) on that hap at all (this may be redundant!)
            ms_intr_sites_hap1 = [s for s in ms_N_sites_hap1
                                  if s['haplotypes_1_intr'][ind]]
            ms_intr_sites_hap2 = [s for s in ms_N_sites_hap2
                                  if s['haplotypes_2_intr'][ind]]

            ms_fields += [len(ms_intr_sites_hap1),
                          len(ms_intr_sites_hap2),
                          len(ms_ss_sites_hap1),
                          len(ms_ss_sites_hap2),
                          len(ms_intr_ss_sites_hap1),
                          len(ms_intr_ss_sites_hap2),
                          len(ms_N_sites_hap1),
                          len(ms_N_sites_hap2)]

            ms_fields_labels += ['ms_intr_sites_hap1',
                                 'ms_intr_sites_hap2',
                                 'ms_ss_sites_hap1',
                                 'ms_ss_sites_hap2',
                                 'ms_intr_ss_sites_hap1',
                                 'ms_intr_ss_sites_hap2',
                                 'ms_N_sites_hap1',
                                 'ms_N_sites_hap2']

        if opts.first_line:
            print('\t'.join([
                'chrom', 'winstart', 'winend', 'n_snps', 'n_ind_snps',
                'n_region_ind_snps', 'ind_id', 'pop', 's_star',
                'num_s_star_snps', 's_star_snps',
                'hap_1_window_pval', 'hap_2_window_pval',
                'hap_1_window_match_pct', 'hap_2_window_match_pct',
                'hap_1_window_pval_local', 'hap_2_window_pval_local',
                'hap_1_window_match_pct_local', 'hap_2_window_match_pct_local',
                'hap_1_window_pval_table', 'hap_2_window_pval_table',
                'hap_1_window_match_pct_table', 'hap_2_window_match_pct_table',
                'hap_1_window_match_N_table', 'hap_2_window_match_N_table',
                'hap_1_window_match_f_table', 'hap_2_window_match_f_table',
                'hap_1_window_match_len_table', 'hap_2_window_match_len_table',
                'hap_1_window_match_mapped_table',
                'hap_2_window_match_mapped_table',
                'hap_1_window_match_mh_table', 'hap_2_window_match_mh_table',
                'hap_1_window_match_sfs_table', 'hap_2_window_match_sfs_table',
                'hap_1_window_match_orig_sfs_table',
                'hap_2_window_match_orig_sfs_table',
                'hap_1_window_match_arc_table', 'hap_2_window_match_arc_table',
                'hap_1_window_match_tot_sites_table',
                'hap_2_window_match_tot_sites_table',
                'hap_1_s_start', 'hap_1_s_end', 'hap_2_s_start', 'hap_2_s_end',
                's_start', 's_end', 'ref_nocals',
                'n_s_star_snps_hap1', 'n_s_star_snps_hap2', 's_star_haps',
                'callable_bases'] + ms_fields_labels + opts.tag_ids))

            opts.first_line = False

        print('\t'.join(str(s) for s in [
            chrom, winstart, winend, len(snps), len(ind1_snps),
            len(ind1_or_ref_snps), opts.get_id_from_sample_index(ind),
            opts.get_pop_from_sample_index(ind), s_star, len(s_star_snps),
            ','.join(str(ind1_pos[i]) for i in s_star_snps)
            if len(s_star_snps) > 0 else '.',
            match_pvals2[ind*2], match_pvals2[ind*2+1],
            match_pct2[ind*2], match_pct2[ind*2+1],
            match_pvals3[ind*2], match_pvals3[ind*2+1],
            match_pct3[ind*2], match_pct3[ind*2+1],
            match_table_pvals[0], match_table_pvals[1],
            match_table_pcts[0], match_table_pcts[1],
            match_table_Ns[0], match_table_Ns[1],
            match_table_factors[0], match_table_factors[1],
            match_table_query[0][0], match_table_query[1][0],
            match_table_query[0][1], match_table_query[1][1],
            match_table_query[0][2], match_table_query[1][2],
            match_table_query[0][3], match_table_query[1][3],
            match_table_query[0][4], match_table_query[1][4],
            match_table_query[0][5], match_table_query[1][5],
            match_table_query[0][6], match_table_query[1][6],
            hap1_range[0], hap1_range[1],
            hap2_range[0], hap2_range[1],
            s_start, s_stop, ref_nones, n_haps1, n_haps2,
            ','.join(str(s) for s in all_haps)
            if len(s_star_snps) > 0 else '.',
            my_mapped_bases] + ms_fields + opts.tags))


def finish_analysis(opts):
    return
