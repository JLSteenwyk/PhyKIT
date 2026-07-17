import json
from pathlib import Path
import sys
from unittest.mock import patch

from Bio import SeqIO
import pytest

from phykit.phykit import Phykit


SAMPLE_FILES = Path(__file__).parents[2] / "sample_files"
PROTEIN_PATH = SAMPLE_FILES / "OG0002774.aln.afa.txt"
NUCLEOTIDE_PATH = SAMPLE_FILES / "OG0002774.mrna.fa.txt"
CLIPKIT_LOG_PATH = SAMPLE_FILES / "OG0002774.aln.afa.clipkit.log.txt"
GAP_CHARS = frozenset("-?*Xx")


def thread_codons_oracle(protein, nucleotide, keep_mask, keep_terminal_stop=True):
    aligned = []
    nucleotide_offset = 0
    for index, amino_acid in enumerate(protein):
        terminal_stop = (
            keep_terminal_stop
            and amino_acid == "*"
            and index == len(protein) - 1
        )
        if amino_acid in GAP_CHARS and not terminal_stop:
            aligned.extend("---")
            continue

        codon = nucleotide[nucleotide_offset:nucleotide_offset + 3]
        nucleotide_offset += 3
        aligned.extend(codon.ljust(3, "-"))

    return "".join(
        base
        for base, keep in zip(aligned, keep_mask)
        if keep
    )


def expected_rows(protein_path, nucleotide_path, clipkit_log_path=None):
    proteins = SeqIO.to_dict(SeqIO.parse(protein_path, "fasta"))
    nucleotides = SeqIO.to_dict(SeqIO.parse(nucleotide_path, "fasta"))
    alignment_length = len(next(iter(proteins.values())).seq)
    if clipkit_log_path is None:
        keep_mask = [True] * (alignment_length * 3)
    else:
        with open(clipkit_log_path) as handle:
            keep_mask = [
                status == "keep"
                for line in handle
                for status in [line.rstrip("\n").split(" ")[1]]
                for _ in range(3)
            ]

    return [
        {
            "taxon": taxon,
            "sequence": thread_codons_oracle(
                str(record.seq),
                str(nucleotides[taxon].seq),
                keep_mask,
            ),
        }
        for taxon, record in proteins.items()
    ]


def run_command(command, protein_path, nucleotide_path, *extra_args):
    test_args = [
        "phykit",
        command,
        "-p",
        str(protein_path),
        "-n",
        str(nucleotide_path),
        *map(str, extra_args),
    ]
    with patch.object(sys, "argv", test_args), patch("builtins.print") as mocked_print:
        Phykit()
    mocked_print.assert_called_once()
    return mocked_print.call_args.args[0]


@pytest.mark.integration
class TestDNAThreader:
    @pytest.mark.parametrize("command", ["thread_dna", "p2n", "pal2nal"])
    def test_commands_match_codon_oracle_with_clipkit_mask(self, command):
        observed = run_command(
            command,
            PROTEIN_PATH,
            NUCLEOTIDE_PATH,
            "-c",
            CLIPKIT_LOG_PATH,
        )
        rows = expected_rows(PROTEIN_PATH, NUCLEOTIDE_PATH, CLIPKIT_LOG_PATH)
        expected = "\n".join(
            f">{row['taxon']}\n{row['sequence']}"
            for row in rows
        )

        assert observed == expected

    def test_command_matches_codon_oracle_without_clipkit_mask(self):
        observed = run_command("pal2nal", PROTEIN_PATH, NUCLEOTIDE_PATH)
        rows = expected_rows(PROTEIN_PATH, NUCLEOTIDE_PATH)
        expected = "\n".join(
            f">{row['taxon']}\n{row['sequence']}"
            for row in rows
        )

        assert observed == expected

    def test_json_output_matches_codon_oracle(self):
        observed = run_command(
            "thread_dna",
            PROTEIN_PATH,
            NUCLEOTIDE_PATH,
            "-c",
            CLIPKIT_LOG_PATH,
            "--json",
        )
        payload = json.loads(observed)
        rows = expected_rows(PROTEIN_PATH, NUCLEOTIDE_PATH, CLIPKIT_LOG_PATH)

        assert payload["clipkit_log_file"].endswith(CLIPKIT_LOG_PATH.name)
        assert payload["remove_stop_codon"] is True
        assert payload["rows"] == rows
        assert payload["taxa"] == rows

    def test_internal_protein_gap_preserves_codon_position(self, tmp_path):
        protein_path = tmp_path / "protein.fa"
        nucleotide_path = tmp_path / "nucleotide.fa"
        protein_path.write_text(">gene1\nA-C\n")
        nucleotide_path.write_text(">gene1\nAAACCC\n")

        observed = run_command(
            "thread_dna",
            protein_path,
            nucleotide_path,
            "-s",
            "false",
        )

        assert observed == ">gene1\nAAA---CCC"
