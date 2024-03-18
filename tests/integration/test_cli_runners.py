import os
import pytest
import subprocess
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)

@pytest.mark.integration
class TestBrokenPipeError(object):
    @pytest.mark.slow
    def test_pk_alignment_length(self):
        cmd = "pk_alignment_length -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_aln_len(self):
        cmd = "pk_aln_len -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_al(self):
        cmd = "pk_al -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_alignment_length_no_gaps(self):
        cmd = "pk_alignment_length_no_gaps -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_aln_len_no_gaps(self):
        cmd = "pk_aln_len_no_gaps -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_alng(self):
        cmd = "pk_alng -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_column_score(self):
        cmd = "pk_column_score -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_cs(self):
        cmd = "pk_cs -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_faidx(self):
        cmd = "pk_faidx -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_get_entry(self):
        cmd = "pk_get_entry -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_ge(self):
        cmd = "pk_ge -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_gc_content(self):
        cmd = "pk_gc_content -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_gc(self):
        cmd = "pk_gc -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pairwise_identity(self):
        cmd = "pk_pairwise_identity -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pairwise_id(self):
        cmd = "pk_pairwise_id -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pi(self):
        cmd = "pk_pi -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_parsimony_informative_sites(self):
        cmd = "pk_parsimony_informative_sites -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pis(self):
        cmd = "pk_pis -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_relative_composition_variability(self):
        cmd = "pk_relative_composition_variability -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rel_comp_var(self):
        cmd = "pk_rel_comp_var -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rcv(self):
        cmd = "pk_rcv -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_relative_composition_variability_taxon(self):
        cmd = "pk_relative_composition_variability_taxon -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_rel_comp_var_taxon(self):
        cmd = "pk_rel_comp_var_taxon -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rcvt(self):
        cmd = "pk_rcvt -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rename_fasta_entries(self):
        cmd = "pk_rename_fasta_entries -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rename_fasta(self):
        cmd = "pk_rename_fasta -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_sum_of_pairs_score(self):
        cmd = "pk_sum_of_pairs_score -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_sops(self):
        cmd = "pk_sops -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_sop(self):
        cmd = "pk_sop -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_variable_sites(self):
        cmd = "pk_variable_sites -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_vs(self):
        cmd = "pk_vs -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_bipartition_support_stats(self):
        cmd = "pk_bipartition_support_stats -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_bss(self):
        cmd = "pk_bss -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_branch_length_multiplier(self):
        cmd = "pk_branch_length_multiplier -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_blm(self):
        cmd = "pk_blm -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_collapse_branches(self):
        cmd = "pk_collapse_branches -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_collapse(self):
        cmd = "pk_collapse -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_cb(self):
        cmd = "pk_cb -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_covarying_evolutionary_rates(self):
        cmd = "pk_covarying_evolutionary_rates -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_cover(self):
        cmd = "pk_cover -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_dvmc(self):
        cmd = "pk_dvmc -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_degree_of_violation_of_a_molecular_clock(self):
        cmd = "pk_degree_of_violation_of_a_molecular_clock -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evo_rate(self):
        cmd = "pk_evo_rate -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evolutionary_rate(self):
        cmd = "pk_evolutionary_rate -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_hidden_paralogy_check(self):
        cmd = "pk_hidden_paralogy_check -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_clan_check(self):
        cmd = "pk_clan_check -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_internal_branch_stats(self):
        cmd = "pk_internal_branch_stats -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_ibs(self):
        cmd = "pk_ibs -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_internode_labeler(self):
        cmd = "pk_internode_labeler -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_il(self):
        cmd = "pk_il -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_last_common_ancestor_subtree(self):
        cmd = "pk_last_common_ancestor_subtree -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_lca_subtree(self):
        cmd = "pk_lca_subtree -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_lb_score(self):
        cmd = "pk_lb_score -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_long_branch_score(self):
        cmd = "pk_long_branch_score -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_lbs(self):
        cmd = "pk_lbs -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_monophyly_check(self):
        cmd = "pk_monophyly_check -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_is_monophyletic(self):
        cmd = "pk_is_monophyletic -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_nearest_neighbor_interchange(self):
        cmd = "pk_nearest_neighbor_interchange -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_nni(self):
        cmd = "pk_nni -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_patristic_distances(self):
        cmd = "pk_patristic_distances -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_pd(self):
        cmd = "pk_pd -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_polytomy_test(self):
        cmd = "pk_polytomy_test -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_polyt_test(self):
        cmd = "pk_polyt_test -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_ptt(self):
        cmd = "pk_ptt -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_polyt(self):
        cmd = "pk_polyt -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_print_tree(self):
        cmd = "pk_print_tree -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_print(self):
        cmd = "pk_print -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_pt(self):
        cmd = "pk_pt -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_prune_tree(self):
        cmd = "pk_prune_tree -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_prune(self):
        cmd = "pk_prune -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rename_tree_tips(self):
        cmd = "pk_rename_tree_tips -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rename_tree(self):
        cmd = "pk_rename_tree -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rename_tips(self):
        cmd = "pk_rename_tips -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rf_distance(self):
        cmd = "pk_rf_distance -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_robinson_foulds_distance(self):
        cmd = "pk_robinson_foulds_distance -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rf_dist(self):
        cmd = "pk_rf_dist -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rf(self):
        cmd = "pk_rf -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_root_tree(self):
        cmd = "pk_root_tree -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_root(self):
        cmd = "pk_root -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rt(self):
        cmd = "pk_rt -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_spurious_sequence(self):
        cmd = "pk_spurious_sequence -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_spurious_seq(self):
        cmd = "pk_spurious_seq -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_ss(self):
        cmd = "pk_ss -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_terminal_branch_stats(self):
        cmd = "pk_terminal_branch_stats -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tbs(self):
        cmd = "pk_tbs -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tip_labels(self):
        cmd = "pk_tip_labels -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_labels(self):
        cmd = "pk_labels -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_tip_to_tip_distance(self):
        cmd = "pk_tip_to_tip_distance -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t_dist(self):
        cmd = "pk_t2t_dist -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t(self):
        cmd = "pk_t2t -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tip_to_tip_node_distance(self):
        cmd = "pk_tip_to_tip_node_distance -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t_node_dist(self):
        cmd = "pk_t2t_node_dist -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t_nd(self):
        cmd = "pk_t2t_nd -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_tree_labels(self):
        cmd = "pk_tree_labels -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tl(self):
        cmd = "pk_tl -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_total_tree_length(self):
        cmd = "pk_total_tree_length -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tree_len(self):
        cmd = "pk_tree_len -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_treeness(self):
        cmd = "pk_treeness -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tness(self):
        cmd = "pk_tness -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_saturation(self):
        cmd = "pk_saturation -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_sat(self):
        cmd = "pk_sat -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_treeness_over_rcv(self):
        cmd = "pk_treeness_over_rcv -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_toverr(self):
        cmd = "pk_toverr -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tor(self):
        cmd = "pk_tor -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_create_concatenation_matrix(self):
        cmd = "pk_create_concatenation_matrix -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_create_concat(self):
        cmd = "pk_create_concat -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_cc(self):
        cmd = "pk_cc -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_thread_dna(self):
        cmd = "pk_thread_dna -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_pal2nal(self):
        cmd = "pk_pal2nal -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_p2n(self):
        cmd = "pk_p2n -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_compositional_bias_per_site(self):
        cmd = "pk_compositional_bias_per_site -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_comp_bias_per_site(self):
        cmd = "pk_comp_bias_per_site -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_cbps(self):
        cmd = "pk_cbps -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evolutionary_rate_per_site(self):
        cmd = "pk_evolutionary_rate_per_site -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evo_rate_per_site(self):
        cmd = "pk_evo_rate_per_site -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_erps(self):
        cmd = "pk_erps -h"
        exit_status = os.system(cmd)
        assert exit_status == 0
