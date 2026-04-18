#!/usr/bin/env python3
"""Generate merged GenBank output from Python GECCO for comparison with Rust.

Writes to stdout so the Rust test can capture it.
Must match the clusters constructed in test_genbank_merged_matches_python().
"""

import sys
from io import StringIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from gecco.model import (
    Cluster,
    ClusterType,
    Domain,
    Gene,
    Protein,
    Strand,
)


def make_source(seq_id: str, dna: str) -> SeqRecord:
    """Create a minimal source SeqRecord."""
    rec = SeqRecord(Seq(dna), id=seq_id, name=seq_id)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    return rec


def main():
    source1 = make_source("seq1", "A" * 1000)
    source2 = make_source("seq2", "C" * 1000)

    domain1 = Domain(
        name="PF00394",
        start=10,
        end=50,
        hmm="Pfam",
        i_evalue=1e-10,
        pvalue=1e-14,
    )
    gene1 = Gene(
        source=source1,
        start=100,
        end=400,
        strand=Strand.Coding,
        protein=Protein(id="seq1_1", seq=Seq(""), domains=[domain1]),
        _probability=0.9,
    )

    cluster1 = Cluster(
        id="seq1_cluster_1",
        genes=[gene1],
        type=ClusterType("Polyketide"),
        type_probabilities={"Polyketide": 0.96},
    )

    domain2 = Domain(
        name="PF00109",
        start=15,
        end=60,
        hmm="Pfam",
        i_evalue=1e-10,
        pvalue=1e-14,
    )
    gene2 = Gene(
        source=source2,
        start=200,
        end=600,
        strand=Strand.Reverse,
        protein=Protein(id="seq2_1", seq=Seq(""), domains=[domain2]),
        _probability=0.85,
    )

    cluster2 = Cluster(
        id="seq2_cluster_1",
        genes=[gene2],
        type=ClusterType("NRP"),
        type_probabilities={"NRP": 0.91},
    )

    # Write merged GenBank to stdout (same as _common.write_clusters with merge=True).
    records = [c.to_seq_record() for c in [cluster1, cluster2]]
    buf = StringIO()
    SeqIO.write(records, buf, "genbank")
    sys.stdout.write(buf.getvalue())


if __name__ == "__main__":
    main()
