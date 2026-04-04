import subprocess
import sys
import pytest

@pytest.mark.integration
class TestBrokenPipeError(object):
    @staticmethod
    def _run_cli(command_name):
        command = [sys.executable, "-m", "phykit", command_name.replace("pk_", ""), "-h"]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        return completed.returncode

    @pytest.mark.slow
    def test_pk_alignment_length(self):
        exit_status = self._run_cli("pk_alignment_length")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_aln_len(self):
        exit_status = self._run_cli("pk_aln_len")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_al(self):
        exit_status = self._run_cli("pk_al")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_alignment_length_no_gaps(self):
        exit_status = self._run_cli("pk_alignment_length_no_gaps")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_aln_len_no_gaps(self):
        exit_status = self._run_cli("pk_aln_len_no_gaps")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_alng(self):
        exit_status = self._run_cli("pk_alng")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_column_score(self):
        exit_status = self._run_cli("pk_column_score")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_cs(self):
        exit_status = self._run_cli("pk_cs")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_faidx(self):
        exit_status = self._run_cli("pk_faidx")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_get_entry(self):
        exit_status = self._run_cli("pk_get_entry")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_ge(self):
        exit_status = self._run_cli("pk_ge")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_gc_content(self):
        exit_status = self._run_cli("pk_gc_content")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_gc(self):
        exit_status = self._run_cli("pk_gc")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pairwise_identity(self):
        exit_status = self._run_cli("pk_pairwise_identity")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pairwise_id(self):
        exit_status = self._run_cli("pk_pairwise_id")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pi(self):
        exit_status = self._run_cli("pk_pi")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_parsimony_informative_sites(self):
        exit_status = self._run_cli("pk_parsimony_informative_sites")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_pis(self):
        exit_status = self._run_cli("pk_pis")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_relative_composition_variability(self):
        exit_status = self._run_cli("pk_relative_composition_variability")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rel_comp_var(self):
        exit_status = self._run_cli("pk_rel_comp_var")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rcv(self):
        exit_status = self._run_cli("pk_rcv")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_relative_composition_variability_taxon(self):
        exit_status = self._run_cli("pk_relative_composition_variability_taxon")
        assert exit_status == 0

    @pytest.mark.slow
    def test_rel_comp_var_taxon(self):
        exit_status = self._run_cli("pk_rel_comp_var_taxon")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rcvt(self):
        exit_status = self._run_cli("pk_rcvt")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rename_fasta_entries(self):
        exit_status = self._run_cli("pk_rename_fasta_entries")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rename_fasta(self):
        exit_status = self._run_cli("pk_rename_fasta")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_sum_of_pairs_score(self):
        exit_status = self._run_cli("pk_sum_of_pairs_score")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_sops(self):
        exit_status = self._run_cli("pk_sops")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_sop(self):
        exit_status = self._run_cli("pk_sop")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_variable_sites(self):
        exit_status = self._run_cli("pk_variable_sites")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_vs(self):
        exit_status = self._run_cli("pk_vs")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_bipartition_support_stats(self):
        exit_status = self._run_cli("pk_bipartition_support_stats")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_bss(self):
        exit_status = self._run_cli("pk_bss")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_branch_length_multiplier(self):
        exit_status = self._run_cli("pk_branch_length_multiplier")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_blm(self):
        exit_status = self._run_cli("pk_blm")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_collapse_branches(self):
        exit_status = self._run_cli("pk_collapse_branches")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_collapse(self):
        exit_status = self._run_cli("pk_collapse")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_cb(self):
        exit_status = self._run_cli("pk_cb")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_covarying_evolutionary_rates(self):
        exit_status = self._run_cli("pk_covarying_evolutionary_rates")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_cover(self):
        exit_status = self._run_cli("pk_cover")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_consensus_tree(self):
        exit_status = self._run_cli("pk_consensus_tree")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_consensus(self):
        exit_status = self._run_cli("pk_consensus")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_ctree(self):
        exit_status = self._run_cli("pk_ctree")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_dtt(self):
        exit_status = self._run_cli("pk_dtt")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_disparity_through_time(self):
        exit_status = self._run_cli("pk_disparity_through_time")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_dvmc(self):
        exit_status = self._run_cli("pk_dvmc")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_degree_of_violation_of_a_molecular_clock(self):
        exit_status = self._run_cli("pk_degree_of_violation_of_a_molecular_clock")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evo_rate(self):
        exit_status = self._run_cli("pk_evo_rate")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evolutionary_rate(self):
        exit_status = self._run_cli("pk_evolutionary_rate")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_hidden_paralogy_check(self):
        exit_status = self._run_cli("pk_hidden_paralogy_check")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_clan_check(self):
        exit_status = self._run_cli("pk_clan_check")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_internal_branch_stats(self):
        exit_status = self._run_cli("pk_internal_branch_stats")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_ibs(self):
        exit_status = self._run_cli("pk_ibs")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_internode_labeler(self):
        exit_status = self._run_cli("pk_internode_labeler")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_il(self):
        exit_status = self._run_cli("pk_il")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_last_common_ancestor_subtree(self):
        exit_status = self._run_cli("pk_last_common_ancestor_subtree")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_lca_subtree(self):
        exit_status = self._run_cli("pk_lca_subtree")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_lb_score(self):
        exit_status = self._run_cli("pk_lb_score")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_long_branch_score(self):
        exit_status = self._run_cli("pk_long_branch_score")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_lbs(self):
        exit_status = self._run_cli("pk_lbs")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_monophyly_check(self):
        exit_status = self._run_cli("pk_monophyly_check")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_is_monophyletic(self):
        exit_status = self._run_cli("pk_is_monophyletic")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_nearest_neighbor_interchange(self):
        exit_status = self._run_cli("pk_nearest_neighbor_interchange")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_nni(self):
        exit_status = self._run_cli("pk_nni")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_patristic_distances(self):
        exit_status = self._run_cli("pk_patristic_distances")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_pd(self):
        exit_status = self._run_cli("pk_pd")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_polytomy_test(self):
        exit_status = self._run_cli("pk_polytomy_test")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_polyt_test(self):
        exit_status = self._run_cli("pk_polyt_test")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_ptt(self):
        exit_status = self._run_cli("pk_ptt")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_polyt(self):
        exit_status = self._run_cli("pk_polyt")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_print_tree(self):
        exit_status = self._run_cli("pk_print_tree")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_print(self):
        exit_status = self._run_cli("pk_print")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_pt(self):
        exit_status = self._run_cli("pk_pt")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_prune_tree(self):
        exit_status = self._run_cli("pk_prune_tree")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_prune(self):
        exit_status = self._run_cli("pk_prune")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rename_tree_tips(self):
        exit_status = self._run_cli("pk_rename_tree_tips")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rename_tree(self):
        exit_status = self._run_cli("pk_rename_tree")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rename_tips(self):
        exit_status = self._run_cli("pk_rename_tips")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rf_distance(self):
        exit_status = self._run_cli("pk_rf_distance")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_robinson_foulds_distance(self):
        exit_status = self._run_cli("pk_robinson_foulds_distance")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rf_dist(self):
        exit_status = self._run_cli("pk_rf_dist")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_rf(self):
        exit_status = self._run_cli("pk_rf")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_root_tree(self):
        exit_status = self._run_cli("pk_root_tree")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_root(self):
        exit_status = self._run_cli("pk_root")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_rt(self):
        exit_status = self._run_cli("pk_rt")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_spurious_sequence(self):
        exit_status = self._run_cli("pk_spurious_sequence")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_spurious_seq(self):
        exit_status = self._run_cli("pk_spurious_seq")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_ss(self):
        exit_status = self._run_cli("pk_ss")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_terminal_branch_stats(self):
        exit_status = self._run_cli("pk_terminal_branch_stats")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tbs(self):
        exit_status = self._run_cli("pk_tbs")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tip_labels(self):
        exit_status = self._run_cli("pk_tip_labels")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_labels(self):
        exit_status = self._run_cli("pk_labels")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_tip_to_tip_distance(self):
        exit_status = self._run_cli("pk_tip_to_tip_distance")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t_dist(self):
        exit_status = self._run_cli("pk_t2t_dist")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t(self):
        exit_status = self._run_cli("pk_t2t")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tip_to_tip_node_distance(self):
        exit_status = self._run_cli("pk_tip_to_tip_node_distance")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t_node_dist(self):
        exit_status = self._run_cli("pk_t2t_node_dist")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_t2t_nd(self):
        exit_status = self._run_cli("pk_t2t_nd")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_tree_labels(self):
        exit_status = self._run_cli("pk_tree_labels")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tl(self):
        exit_status = self._run_cli("pk_tl")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_total_tree_length(self):
        exit_status = self._run_cli("pk_total_tree_length")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tree_len(self):
        exit_status = self._run_cli("pk_tree_len")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_treeness(self):
        exit_status = self._run_cli("pk_treeness")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tness(self):
        exit_status = self._run_cli("pk_tness")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_saturation(self):
        exit_status = self._run_cli("pk_saturation")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_sat(self):
        exit_status = self._run_cli("pk_sat")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_treeness_over_rcv(self):
        exit_status = self._run_cli("pk_treeness_over_rcv")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_toverr(self):
        exit_status = self._run_cli("pk_toverr")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_tor(self):
        exit_status = self._run_cli("pk_tor")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_create_concatenation_matrix(self):
        exit_status = self._run_cli("pk_create_concatenation_matrix")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_create_concat(self):
        exit_status = self._run_cli("pk_create_concat")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_cc(self):
        exit_status = self._run_cli("pk_cc")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_thread_dna(self):
        exit_status = self._run_cli("pk_thread_dna")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_pal2nal(self):
        exit_status = self._run_cli("pk_pal2nal")
        assert exit_status == 0
    
    @pytest.mark.slow
    def test_pk_p2n(self):
        exit_status = self._run_cli("pk_p2n")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_compositional_bias_per_site(self):
        exit_status = self._run_cli("pk_compositional_bias_per_site")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_comp_bias_per_site(self):
        exit_status = self._run_cli("pk_comp_bias_per_site")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_cbps(self):
        exit_status = self._run_cli("pk_cbps")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evolutionary_rate_per_site(self):
        exit_status = self._run_cli("pk_evolutionary_rate_per_site")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_evo_rate_per_site(self):
        exit_status = self._run_cli("pk_evo_rate_per_site")
        assert exit_status == 0

    @pytest.mark.slow
    def test_pk_erps(self):
        exit_status = self._run_cli("pk_erps")
        assert exit_status == 0
