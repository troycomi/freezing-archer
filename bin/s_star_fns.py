from __future__ import division
import sys
import numpy as np
from collections import Counter
import arc_match_pval_tables


def calc_s_dists(genotypes, positions,
                 match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000):
    # all pairwise dists[i, j] from i to j, i <= j
    pos = np.asarray(positions)
    bp_dist = pos[None, :] - pos[:, None]
    gen = np.asarray(genotypes)
    gen = np.abs(gen[None, :] - gen[:, None])
    return np.where(
        bp_dist < 10,
        np.NINF,
        np.where(
            gen == 0,
            match_bonus + bp_dist,
            np.where(
                gen <= max_mismatch,
                mismatch_penalty,
                np.NINF)))


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


def initialize_analysis(opts):
    setattr(opts, 'pval_repeat_lookup', {})


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

        # get haplotype snps (two haplotypes per person)
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
                  opts.get_pop_from_sample_index[ind], n_match1,
                  len(neand_pos), len(ind_pos1), my_n_sites1, test_ratio1)
        if opts.debug:
            print("test_ratio2", winstart, winend, ind,
                  opts.get_pop_from_sample_index[ind], n_match2,
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
                  opts.get_pop_from_sample_index[ind],
                  match_pct[ind*2], pvals[ind*2])
        if opts.debug:
            print("pval2", winstart, winend, ind,
                  opts.get_pop_from_sample_index[ind],
                  match_pct[ind*2+1], pvals[ind*2+1])

    return (pvals, match_pct, ref_nones)


def run_window_analysis(chrom, winstart, winend, snps, opts):
    '''
    snps is list of dictionaries with vcf keys
    '''
    if len(snps) <= 2:
        return

    if opts.archaic_vcf is None or opts.no_pvalues:
        match_pvals2 = ['NA', 'NA'] * opts.num_samples
        match_pct2 = ['NA', 'NA'] * opts.num_samples
        ref_nones = 0
    else:
        (match_pvals2, match_pct2, ref_nones) = calc_local_match_pval(
            chrom, winstart, winend, snps, opts)

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    my_mapped_bases = opts.regions.amount_in_region(
        full_chrom, winstart, winend) if opts.regions is not None \
        else opts.window_length

    # just loop over target individuals
    for ind in opts.target_indices:

        ind_snps_all = [s for s in snps if s['genotypes'][ind] > 0]
        ind_snps = [s for s in ind_snps_all
                    if s['target'] and not s['reference']]

        if len(ind_snps) <= 2:
            # TODO this is much faster, but no longer matches
            # return
            (s_star, s_star_snps) = (0, [])
            ind_pos = []

        else:
            ind_pos = [s['pos'] for s in ind_snps]
            # s_star is the score and s_star_snps are the indices
            # (wrt ind_snps) that are selected by s_star
            (s_star, s_star_snps) = calc_s_star(
                [s['genotypes'][ind] for s in ind_snps],
                ind_pos,
                match_bonus=opts.s_star_match_bonus,
                max_mismatch=opts.s_star_max_mismatch,
                mismatch_penalty=opts.s_star_mismatch_penalty)

        # TODO this is much faster, but no longer matches
        # if s_star == 0:
            # return

        # masks for which snps are on each haplotypes
        ind_hap1 = [s['haplotypes_1'][ind] > 0 for s in ind_snps]
        ind_hap2 = [s['haplotypes_2'][ind] > 0 for s in ind_snps]

        # find the best haplotype
        n_haps1 = sum(ind_hap1[i] for i in s_star_snps)
        n_haps2 = sum(ind_hap2[i] for i in s_star_snps)
        all_haps = (ind_hap1[i] + ind_hap2[i]*2 for i in s_star_snps)

        # get S* haplotype range
        hap1_pos = [ind_snps[i]['pos'] for i in s_star_snps if ind_hap1[i]]
        hap1_range = (min(hap1_pos), max(hap1_pos)) \
            if len(hap1_pos) > 1 else (0, 0)
        hap2_pos = [ind_snps[i]['pos'] for i in s_star_snps if ind_hap2[i]]
        hap2_range = (min(hap2_pos), max(hap2_pos)) \
            if len(hap2_pos) > 1 else (0, 0)

        # TODO remove from output
        if opts.match_pval_table is not None:
            raise ValueError('Match pval table is no longer supported')

        if opts.no_pvalues or len(s_star_snps) < 2 or opts.archaic_vcf is None:
            s_start = 'NA'
            s_stop = 'NA'
            match_pvals3 = ['NA', 'NA']
            match_pct3 = ['NA', 'NA']

        else:
            # get the start and stop of the S* region
            s_start = ind_snps[min(s_star_snps)]['pos']
            s_stop = ind_snps[max(s_star_snps)]['pos']
            (match_pvals3, match_pct3, _) = calc_local_match_pval(
                chrom, winstart, winend, snps, opts, ind,
                subset_start=s_start, subset_end=s_stop)
            match_pvals3 = match_pvals3[ind*2:ind*2+2]
            match_pct3 = match_pct3[ind*2:ind*2+2]

        ms_fields = []
        ms_fields_labels = []
        if opts.vcf_is_ms_file:
            # count sites that are intr / S* / both
            # N: sites that are a) on that hap at all
            ms_N_sites_hap1 = [s for s in ind_snps_all
                               if s['haplotypes_1'][ind] > 0]
            ms_N_sites_hap2 = [s for s in ind_snps_all
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

        # used as a count of the number of snvs in a region
        # for matching to the null model
        ind_or_ref_snps = [s for s in snps
                           if s['genotypes'][ind] > 0 or s['reference']]

        print('\t'.join(str(s) for s in [
            chrom, winstart, winend, len(snps), len(ind_snps),
            len(ind_or_ref_snps), opts.get_id_from_sample_index[ind],
            opts.get_pop_from_sample_index[ind], s_star, len(s_star_snps),
            ','.join(str(ind_pos[i]) for i in s_star_snps)
            if len(s_star_snps) > 0 else '.',
            match_pvals2[ind*2], match_pvals2[ind*2+1],
            match_pct2[ind*2], match_pct2[ind*2+1],
            match_pvals3[0], match_pvals3[1],
            match_pct3[0], match_pct3[1],
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            'NA', 'NA',
            hap1_range[0], hap1_range[1],
            hap2_range[0], hap2_range[1],
            s_start, s_stop, ref_nones, n_haps1, n_haps2,
            ','.join(str(s) for s in all_haps)
            if len(s_star_snps) > 0 else '.',
            my_mapped_bases] + ms_fields + opts.tags))


def finish_analysis(opts):
    return
