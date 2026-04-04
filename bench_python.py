"""Benchmark Python GECCO pipeline stages."""
import sys, os, time
sys.path.insert(0, "GECCO")

# --- 1. Gene finding ---
from gecco.orf import PyrodigalFinder
from Bio import SeqIO

genome = "GECCO/tests/test_orf/data/BGC0001737.fna"
records = list(SeqIO.parse(genome, "fasta"))

t0 = time.perf_counter()
finder = PyrodigalFinder(metagenome=True, cpus=1)
genes = list(finder.find_genes(iter(records)))
t1 = time.perf_counter()
print(f"Gene finding:  {t1-t0:.3f}s  ({len(genes)} genes)")

# --- 2. HMMER annotation ---
from gecco.hmmer import PyHMMER, HMM
from gecco.interpro import InterPro

hmm = HMM(
    id="Pfam", version="test", url="",
    path="GECCO/tests/test_hmmer/data/minipfam.hmm",
    size=10, relabel_with=r"s/(PF\d+).\d+/\1/", md5=None,
)
interpro = InterPro.load()

t0 = time.perf_counter()
genes_ann = PyHMMER(hmm, 1, None).run(genes)
t1 = time.perf_counter()
n_dom = sum(len(g.protein.domains) for g in genes_ann)
print(f"HMMER annot:   {t1-t0:.3f}s  ({n_dom} domains)")

# --- 3. CRF prediction ---
from gecco.crf import ClusterCRF
import pickle

with open("GECCO/gecco/crf/model.pkl", "rb") as f:
    crf = pickle.load(f)

# Use the BGC0001866 features for CRF benchmark (more genes)
from gecco.model import FeatureTable
ft = FeatureTable.load("GECCO/tests/test_cli/data/BGC0001866.features.tsv")
feat_genes = list(ft.to_genes())

t0 = time.perf_counter()
for _ in range(10):
    predicted = list(crf.predict_probabilities(feat_genes))
t1 = time.perf_counter()
print(f"CRF predict:   {(t1-t0)/10:.3f}s  (10 runs avg, {len(feat_genes)} genes)")

