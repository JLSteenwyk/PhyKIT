import subprocess
import sys

import pytest


@pytest.mark.integration
class TestBrokenPipeError:
    @staticmethod
    def _close_producer_stdout(command, lines_to_read):
        process = subprocess.Popen(
            [sys.executable, "-u", "-m", "phykit", *command],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert process.stdout is not None
        assert process.stderr is not None

        lines = [process.stdout.readline() for _ in range(lines_to_read)]
        process.stdout.close()
        stderr = process.stderr.read()
        returncode = process.wait(timeout=30)
        process.stderr.close()

        assert all(line for line in lines)
        assert returncode == 0, stderr
        assert "BrokenPipeError" not in stderr

    @pytest.mark.slow
    def test_dna_threader_broken_pipe(self):
        self._close_producer_stdout(
            [
                "thread_dna",
                "-p",
                "./tests/sample_files/EOG091N44MS.fa.mafft",
                "-n",
                "./tests/sample_files/EOG091N44MS.fa",
            ],
            2,
        )

    @pytest.mark.slow
    def test_create_concatenation_matrix_broken_pipe(self, tmp_path):
        self._close_producer_stdout(
            [
                "create_concatenation_matrix",
                "-a",
                "./tests/sample_files/alignment_list_for_create_concat_matrix.txt",
                "-p",
                str(tmp_path / "concatenated"),
            ],
            2,
        )

    @pytest.mark.slow
    @pytest.mark.parametrize(
        ("command", "lines_to_read"),
        [
            pytest.param(
                ["gc_content", "./tests/sample_files/simple.fa", "-v"],
                1,
                id="gc-content",
            ),
            pytest.param(
                [
                    "pi",
                    "./tests/sample_files/12_YPR189W_Anc_7.546_codon_aln.fasta.clipkit",
                    "-v",
                ],
                1,
                id="pairwise-identity",
            ),
            pytest.param(
                [
                    "vs",
                    "./tests/sample_files/12_YPR189W_Anc_7.546_codon_aln.fasta.clipkit",
                ],
                0,
                id="variable-sites",
            ),
            pytest.param(
                [
                    "bss",
                    "./tests/sample_files/small_Aspergillus_tre_rooted.tree",
                    "-v",
                ],
                1,
                id="bipartitions",
            ),
        ],
    )
    def test_tabular_command_broken_pipe(self, command, lines_to_read):
        self._close_producer_stdout(command, lines_to_read)
