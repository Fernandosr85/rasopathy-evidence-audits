# RASopathy evidence audits

Two audits of the same RASopathy gene panel, each examining a different point at which a
measurement fails to become clinical evidence.

| | Question | Criterion |
|---|---|---|
| **MAVE coverage audit** | Is the clinical variant matched to the right experiment at all? | PS3/BS3 |
| **Frequency controls and prospective VUS application** | Do the reference data support the evidence strength being assigned? | PP3/BP4 |

The two failure modes are independent, and both are silent: the lookup returns a result, the
estimator returns a number, neither raises an error.

They also call for different responses. Ascertainment is a tooling problem — accession-aware
matching fixes it. The resolution limit is not: no method produces benign controls that do not
exist.

---

## 1. MAVE coverage audit

**`notebooks/mave_coverage_audit.ipynb`** · **`reports/mave_report.md`**

Does multiplexed assay of variant effect (MAVE) data exist for the 40 panel variants that are
curated in ClinVar but remain VUS with no functional criterion applied?

A first pass matching HGVS protein changes against MaveDB's `hgvs_pro` column reported zero
coverage. That was wrong, and establishing why turned up five ways a plausible cross-reference
fails silently.

**Result.** Three of forty have isolated MAVE measurements — HRAS I93V and G138S, and BRAF
I208V, which an earlier manual pass had missed. Two more appear only inside multi-substitution
haplotypes. The HGVS route found one variant, and that one was a false positive.

**Usable coverage is zero.** The HRAS measurements come from a bacterial two-hybrid system
without membrane context; the BRAF measurement from a domainome screen of *domain stability in
yeast*, which is orthogonal to a gain-of-function mechanism. No covering score set carries a
published calibration.

### The five failure modes

1. **Depositor numbering.** Score sets number against the assayed construct. HRAS uses the
   mature protein; RAF1 and BRAF cover isolated domains, offsets +51 and +476.
2. **Family-level target names.** `urn:mavedb:00000057` is named `RAS`, not HRAS;
   `00000061` is `RAF`, not RAF1. Target-name metadata provides no usable HGNC key.
3. **Search misses real records.** Searching MaveDB for HRAS returns zero score sets although
   four exist. Four score sets covering panel genes were invisible to symbol search here.
4. **Haplotypes.** `p.[Leu133Ala;Met170Leu]` fails both ways: a VRS reader ignoring `members`
   drops it; a substring match on `hgvs_pro` misattributes the combination score to a single
   variant. 29,481 of 99,044 substitutions in this panel sit inside haplotypes.
5. **Isoform mismatch.** KRAS `00000115` writes `p.Met170Leu` against isoform 4B while its own
   post-mapped block maps to `NP_203524.1`, isoform 4A — where the variant is not coding at
   all (`n.*62A>T`) and residue 170 is lysine. The only false-positive coverage call here:
   seven score sets would have been recorded as measuring a variant absent from the protein
   they measured.

### Four routes, one question

| Route | Covered | Variants | False positives | Fails on |
|---|---:|---|---:|---|
| A. Symbol search only | 1/40 | I208V | 0 | FM2, FM3 |
| B. `hgvs_pro` substring | 1/40 | M170L | **1** | FM1, FM4, FM5 |
| C. Post-mapped, top level only | 3/40 | G138S, I208V, I93V | 0 | FM4 |
| D. Accession-aware | **3/40** | G138S, I208V, I93V | 0 | — |

Three distinct answer sets. A and B each return exactly one variant and share none.

### What resolves it

Matching on `(protein accession, reference residue, position, alternate residue)`, with the
clinical side from the ClinGen Allele Registry, which returns every protein representation of a
variant across isoforms — coding and non-coding alike.

MaveDB already maps assay variants to human reference sequences and retains pre- and
post-mapped VRS representations; MaveMD resolves variants through the Allele Registry so that
transcript and coordinate-system differences need not be handled by the user. None of that is
new. What this audit measures is how far a naive cross-reference goes wrong when that
information is not used.

---

## 2. Frequency controls and prospective VUS application

**`notebooks/frequency_controls_and_prospective_vus.ipynb`** ·
**`reports/preprint_skeleton.md`**

Can panel-specific computational-evidence thresholds be derived for ten RAS-MAPK genes, and at
what evidence strength?

**The limiting factor was not predictor discrimination, but empirical resolution in the benign
tail.** AUROC is 0.963 (EVEE) and 0.980 (REVEL) on a 331-pathogenic / 104-benign reference set.

### Benign tail support at published REVEL thresholds

| Threshold | Benign above | 95% CI |
|---|---:|---|
| Pejaver supporting 0.644 | 4/104 | 1.06 – 9.56% |
| RASopathy VCEP PP3 0.70 | 3/104 | 0.60 – 8.20% |
| Pejaver moderate 0.773 | 1/104 | 0.02 – 5.24% |
| **Pejaver strong 0.932** | **0/104** | **0 – 3.48%** |

### Resolution bound

Since `LR+ = TPR/FPR` and `TPR ≤ 1`, the constraint `FPR ≤ 1/LR_target` binds at any threshold
and any sensitivity. Even observing zero false positives, the sample must be large enough for
the Clopper–Pearson upper bound to fall below that ceiling.

| Target | Max admissible FPR | n needed | Resolvable at n = 104 |
|---|---:|---:|---|
| Tavtigian strong (18.7) | 5.35% | 68 | yes |
| **ClinGen strong (33.53)** | **2.98%** | **122** | **no** |
| ClinGen very strong (1,124) | 0.089% | 4,145 | no |

In the observed configuration the requirement is harsher still: at REVEL 0.932 the observed
sensitivity is 0.275, so the ClinGen strong target needs FPR ≤ 0.82% — below what a *single*
benign variant among 104 would produce (0.96%).

### Derived thresholds are operationally nonselective

Applied to 1,146 unselected ClinVar missense VUS in the same genes, the supporting threshold
derived here would qualify 83.2% of them for REVEL and 77.7% for EVEE. Published supporting
thresholds qualify 24.3% (Pejaver) and 19.5% (RASopathy VCEP). The derived values are
therefore not candidate clinical cutoffs.

### Corrections applied

Haldane–Anscombe on all four cells (the earlier version used `(TP+0.5)/(Npos+0.5)`); no
rounding before the threshold scan; VCEP frequency criteria are filtering allele frequency
≥ 0.05% (BA1) and ≥ 0.025% (BS1), not joint AF; the false-positive rate mapped onto the fixed
benign distribution rather than the resampled one. The corrections had negligible effect on the
primary estimates, which is itself worth reporting.

---

## Running

Both notebooks are self-contained and cache their API responses. Neither writes outside its own
output directory.

**MAVE audit.** Reads `blocked_vus_worksheet.csv` from the
[RASopathy Functional Evidence Atlas](https://github.com/Fernandosr85/rasopathy-functional-evidence-atlas),
downloading it automatically if not mounted. Queries `api.mavedb.org` and
`reg.clinicalgenome.org` — about 40 Allele Registry calls and 37 score-set downloads.

**Frequency controls.** Reads `1_conjunto_comparavel.csv` and `tier1_com_evee.csv` from the
EVEE notebook series, plus a cached EVEE response directory. Queries the gnomAD GraphQL API and
downloads ClinVar `variant_summary.txt.gz` (~440 MB). Step 12 makes new EVEE calls for the
unselected VUS cohort if they are not already cached.

Outputs land in `mave_coverage/` and `recomputation/` as numbered CSVs plus a `manifest.json`
recording each run's parameters and caveats.

---

## Scope and caveats

Nothing here constitutes a variant classification.

**The MAVE audit** establishes whether an assay exists and which protein it measured. Whether
the assay bears on pathogenicity depends on its design and is read from the source publication,
not from identifiers. Conflating the two is how a stability screen becomes functional evidence
for a gain-of-function disorder.

**The calibration study** measures support and stability, not clinical cutoffs. The thresholds
it derives are demonstrably more permissive than published calibrations. It concerns PP3/BP4
only; nothing in it is functional evidence.

Both are single-panel. Failure modes 2 and 5 are plausibly general — family-level naming and
isoform divergence are not specific to RAS-MAPK — but that has not been tested elsewhere.

---

## Related

[RASopathy Functional Evidence Atlas](https://github.com/Fernandosr85/rasopathy-functional-evidence-atlas)
— the curation audit these two grew out of, including the 40-variant worksheet.

## License

Code CC0. ClinVar, gnomAD, MaveDB and ClinGen Allele Registry data are subject to their own
terms.
