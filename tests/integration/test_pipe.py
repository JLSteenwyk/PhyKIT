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
    def test_dna_threader_BrokenPipeError(self):
        cmd = "phykit thread_dna -p ./tests/sample_files/EOG091N44MS.fa.mafft -n ./tests/sample_files/EOG091N44MS.fa | head -n 2"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_create_concatenation_matrix_BrokenPipeError(self):
        cmd = "phykit create_concatenation_matrix -a ./tests/sample_files/alignment_list_for_create_concat_matrix.txt -p test | head -n 2"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_gc_content_BrokenPipeError(self):
        cmd = "phykit gc_content ./tests/sample_files/EOG091N44MS.fa.mafft -v | head -n 1"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_pairwise_identity_BrokenPipeError(self):
        cmd = "phykit pi ./tests/sample_files/12_YPR189W_Anc_7.546_codon_aln.fasta.clipkit -v | head -n 1"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_variable_sites_BrokenPipeError(self):
        cmd = "phykit vs ./tests/sample_files/12_YPR189W_Anc_7.546_codon_aln.fasta.clipkit -v | head -n 1"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_bipartition_support_stats_BrokenPipeError(self):
        cmd = "phykit bss ./tests/sample_files/small_Aspergillus_tre_rooted.tree -v | head -n 1"
        exit_status = os.system(cmd)
        assert exit_status == 0

