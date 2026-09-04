# Auditing MAVE Coverage for RASopathy VUS Reveals Silent Variant-Matching Failure Modes

*Report v4. Numbers from `mave_coverage_audit.ipynb`, retrieved 2 September 2026. Companion
to `frequency_controls_report.md`, which covers the computational-evidence side of the same
panel.*

*MaveDB is a live resource — score sets are revised and remapped, and at least one covering
this panel carries a 2026 modification date. Every count below is therefore a snapshot. The
audit was run twice, twelve days apart, and returned identical figures: 99,044 substitutions,
29,481 of them in haplotypes, 37 score sets, 4 invisible to symbol search, 3/40 covered. That
is a statement about these twelve days, not a guarantee.*

**Central claim:**

> Naive gene- and HGVS-based cross-referencing can silently misstate MAVE coverage. In this
> 40-variant RASopathy audit, an accession-aware comparison of post-mapped MAVE
> representations with ClinGen protein representations recovered three isolated measurements
> missed by HGVS matching and rejected the sole HGVS hit as an isoform-mismatch false
> positive. **None of the three recovered measurements currently provides calibrated,
> mechanism-appropriate evidence for PS3/BS3 application in this disease context.**

**The framing that survives every failure mode:**

> The question is not whether a variant string matches. It is whether the clinical variant
> exists on the protein that generated the score.

**What this is not.** MaveDB already maps assay variants to human reference sequences and
retains pre- and post-mapped VRS representations; MaveMD resolves variants through the
ClinGen Allele Registry precisely so that transcript, identifier and coordinate-system
differences need not be handled by the user. That is established, and nothing here claims
otherwise. This report measures how far a plausible-looking cross-reference goes wrong on a
real clinical panel when that mapped/VRS/CAID information is not used.

---

## 1. Motivation

Forty variants in the RASopathy panel carry a ClinVar curation record, remain classified as
VUS, and have no functional criterion applied. If functional evidence exists for any of them
and simply has not been located, that is actionable: it is the difference between "no data"
and "data not found".

A first pass looked them up in MaveDB by HGVS protein change and returned zero coverage. That
answer was wrong. Establishing why turned up five distinct ways a naive cross-reference fails
silently — and produced a measurement of how far apart four reasonable-looking routes land.

---

## 2. Results

### 2.1 Coverage

**Table 1. Four cross-referencing routes, same 40 variants.**

| Route | Covered | Variants returned | False positives | Fails on |
|---|---:|---|---:|---|
| A. Symbol search only | 1/40 | I208V | 0 | FM2, FM3 |
| B. `hgvs_pro` substring | 1/40 | M170L | **1** | FM1, FM4, FM5 |
| C. Post-mapped, top level only | 3/40 | G138S, I208V, I93V | 0 | FM4 |
| D. Accession-aware | **3/40** | G138S, I208V, I93V | 0 | — |

Three distinct answer sets across four routes. The sharpest disagreement is between A and B:
both return exactly one variant, they share none, and B's is the false positive. Two
cross-references that a reader would consider equally reasonable, run against the same forty
variants, agree on nothing.

C and D coincide here because no variant among the forty is measured only inside a haplotype
in a score set that C would otherwise reach; the distinction still matters, since 29,481 of
99,044 substitutions in the panel sit inside haplotypes.

Additionally, **2/40** are measured only inside multi-substitution haplotypes, which is
neither "measured" nor "no assay".

**Table 2. Variants with isolated MAVE measurements.**

| Gene | ClinVar | Change | Protein matched | Score sets |
|---|---:|---|---|---|
| HRAS | 40439 | I93V | NP_005334.1, NP_789765.1 | `00000057-a/b/c-1` |
| HRAS | 45304 | G138S | NP_005334.1, NP_789765.1 | `00000057-a/b/c-1` |
| BRAF | 178963 | I208V | NP_001341538.1 | `00000708-a-1` |

The BRAF hit is new — it was not recovered by the earlier manual pass. Its score set declares
the offset explicitly (`UniProt: P15056 (offset 152)`, 76-residue construct), which lets the
match be verified by hand: canonical residue 208 is construct position 56, the declared
sequence carries isoleucine there, and the score set writes `p.Ile56Val`. The original route
searched for `p.Ile208Val` and found nothing. FM1 in a single line.

None of the five relevant score sets carries a published calibration. All are
`uncalibrated`; `urn:mavedb:00000708-a-1` reports `Assay OddsPath: Not provided`. Applying
PS3/BS3 would therefore require calibrating from scratch against clinical controls.

### 2.2 Recovered is not usable

**Table 3. Assay suitability of the three recovered measurements.**

| Variants | Score sets | Assay | Why it does not support classification here |
|---|---|---|---|
| HRAS I93V, G138S | `00000057-a/b/c-1` | Bacterial two-hybrid, RAS–effector binding | No membrane context, single effector; ITC validation explicitly noisier near wild-type binding, which is where both variants sit |
| BRAF I208V | `00000708-a-1` | DHFR protein-fragment complementation in yeast; domain **stability/abundance** | Measures folding stability, not RAS-MAPK signalling. RASopathy mechanism is typically gain of function, so a normal stability score is not benign evidence and an abnormal one points the wrong way |

The BRAF score set also leaves `Assay type`, `Molecular mechanism`, `Variant consequences`
and `Model system` unspecified in its metadata, so the mechanistic mismatch is only visible
from the abstract and method text.

**Coverage is 3/40. Coverage that is directly applicable as PS3/BS3 under the evaluated
calibration and mechanism requirements is 0/40.** The distinction matters more than the count:
a table reporting three variants as "assay available" without the mechanism column would
invite exactly the wrong inference.

The phrasing is deliberate. Each of these assays is sound for the question it was designed to
answer — the BRAF screen measures domain stability, and does so well. What none of them
supplies is calibrated evidence about a gain-of-function mechanism, which is what PS3/BS3
application in this panel would require.

### 2.3 The five failure modes

**FM1 — Depositor numbering.** Score sets number positions against the assayed construct.
HRAS score sets use the mature protein; RAF1 and BRAF depositions cover isolated domains with
offsets of +51 and +476. An HGVS lookup against `hgvs_pro` returns a confident zero.
*Signature:* missed coverage, no error.

**FM2 — Family-level target names.** The target of `urn:mavedb:00000057` is named `RAS`, of
`urn:mavedb:00000061` `RAF`. Target-name metadata therefore provides no usable HGNC key for
these score sets — not by equality, not by substring, not by token.
*Signature:* missed coverage, no error.

**FM3 — Full-text search misses real records.** Searching MaveDB for `HRAS` returns zero
score sets although four exist. This follows from FM2. In this audit, **4 score sets** covering
panel genes were invisible to search and had to be supplied by URN.
*Signature:* an inventory that under-reports coverage.

**FM4 — Haplotypes.** Records such as `p.[Leu133Ala;Met170Leu]` are multi-substitution
haplotypes. A substring match on `hgvs_pro` finds `Met170Leu` inside them and attributes a
combination score to a single variant. In VRS these are `Haplotype` objects whose `location`
and `state` live one level below the top, so a naive reader of `postMapped` drops them
entirely — **29,481 of 99,044 substitutions** in this panel sit inside haplotypes.
*Signature:* fails in two directions — a VRS reader ignoring `members` drops the record; a
substring match on `hgvs_pro` misattributes the combination score to a single variant.

**FM5 — Isoform mismatch.** KRAS `urn:mavedb:00000115` writes `p.Met170Leu` in `hgvs_pro`,
numbering against isoform 4B, while its own `postMapped` block maps to `NP_203524.1`, isoform
4A. On 4A the variant is not coding: the ClinGen Allele Registry renders it
`NP_203524.1:n.*62A>T`, sixty-two bases past the stop codon. The residue at position 170 in
the assayed protein is lysine, not methionine.
*Signature:* **false-positive coverage call.** A gene-level match assigned scores from
**seven score sets** to a variant absent from the protein measured.

FM5 was the only failure mode in this audit that produced a false-positive isolated-assay
coverage call. It did not change a classification — none of the score sets involved carries a
published calibration, so nothing here could have supported PS3/BS3 either way. What it moved
was coverage ascertainment, in the more dangerous direction. It is also the hardest to notice:
nothing in the record flags the disagreement between its own two representations.

External confirmation: ClinVar records Variation ID 180859 as `NM_004985.5:c.508A>T →
NP_004976.2:p.Met170Leu`, while on `NM_033360.4` — the transcript encoding `NP_203524.1` — the
same genomic variant is `c.*62A>T`, in the 3' UTR.

### 2.4 A related observation

HRAS I93V shows the same isoform disagreement (`NP_001304983.1` non-coding). It did not
produce a false positive only because the match was keyed on protein accession. The pattern
is not a KRAS peculiarity.

`urn:mavedb:00000057-d-1` ("G12V mutated Ras") appears exclusively as haplotypes — expected,
since every measurement in that score set is the tested substitution combined with the G12V
background. Those scores are not isolated evidence for any single variant.

---

## 3. Method

### 3.1 The route that works

All five failures share a root in the *dataset metadata*: target names are free text,
`hgvs_pro` follows the depositor's construct, and search indexes that text. None of them
identifies the protein that produced the score. The field that does is
`postMapped.location.sequenceReference.label`.

Stable identifiers do exist elsewhere in the ecosystem — VRS digests and ClinGen Allele IDs,
which is what MaveMD searches on. The point is not that they are missing; it is that a user
working from dataset metadata alone will not encounter them.

The match is therefore keyed on:

```
(protein accession, reference residue, position, alternate residue)
```

The clinical side of the key comes from the ClinGen Allele Registry, which returns every
protein representation of a variant across isoforms, coding and non-coding alike. Both sides
of the key are anchored to a RefSeq protein, so isoform divergence cannot produce a spurious
hit.

This is an accession-aware audit of the same correspondence MaveMD resolves through CAIDs.

### 3.2 Diagnostic that precedes the match

Before matching, the audit flags variants that are coding on some isoforms and non-coding on
others. These are exactly the cases where a gene-level match can attribute the wrong score,
and the flag is available without the error having to occur first.

### 3.3 What the audit does not do

It establishes that an assay exists and which protein it measured. It does not judge fitness
for germline classification. For the HRAS score sets recovered here, the source study is a
bacterial two-hybrid system without membrane context, using one effector, with ITC validation
explicitly noisier near wild-type binding — which is where both recovered variants sit. Good
data for the question it was built for; not, without new calibration, for germline
classification.

---

## 4. Implications

**For anyone cross-referencing ClinVar against MAVE data at scale.** An HGVS-based lookup
produces silent missed coverage (FM1–FM3), silent misattribution (FM4), and false-positive
coverage calls (FM5). Using the mapped representation and resolving both sides to a protein
accession avoids all five. That information is available — through post-mapped VRS records,
CAIDs, and MaveMD — but it is not what a user reaches for first, and the failure is silent
when they do not.

**For a VCEP evaluating functional evidence.** Coverage in this panel is genuinely sparse —
three of forty, in two genes, no published calibration on any covering score set, and no assay
whose design fits the disease mechanism. The sparseness must nonetheless be established
through the protein-anchored route: the naive route reports one variant covered, and that one
is wrong.

**Two questions, kept apart.** *Does a measurement exist and is it correctly attributed?* is
answered by identifier handling, and this audit answers it. *Does the measurement bear on
pathogenicity in this disease?* is answered by assay design, and identifier handling cannot
touch it. Conflating them is how a stability screen becomes functional evidence for a
gain-of-function disorder.

**On the count itself.** Three of forty is not the interesting number. The interesting result
is Table 1: four reasonable-looking routes return three different answer sets, two of them
disjoint from each other, and none of the disagreements surfaces as an error.

---

## 5. Limitations

**5.1 Protein accession table is manual.** `NP_GENE` maps the ten accessions encountered here.
A score set referencing another protein is reported as unmapped rather than silently dropped,
but the table needs extending for other panels.

**5.2 Mapping completeness is not verified.** Some score sets return HTTP 404 on
`/mapped-variants`, meaning no published mapping exists. For those, the protein-anchored route
is unavailable and no coverage claim can be made either way.

**5.3 Single panel.** All five failure modes were observed in RAS-MAPK genes. FM2 and FM5 are
plausibly general — family-level naming and isoform divergence are not specific to this
pathway — but that has not been tested elsewhere.

**5.4 Calibration state is read from metadata.** A score set reported as `uncalibrated` may
have a calibration published elsewhere and not deposited.

---

## 6. Data and code

`mave_coverage_audit.ipynb`, with cached API responses and all derived tables
(`1_variants_resolved.csv` through `6_vcep_table.csv`).

Companion analysis on the computational-evidence side of the same panel:
`frequency_controls_report.md`, which examines the independent problem of
whether reference data provide sufficient empirical resolution for computational-evidence
calibration.

**How the two relate.** This audit asks whether the clinical variant has been matched to the
right experiment; the companion asks whether the available reference data support the evidence
strength being assigned. Availability and usability are separate problems, and both failures
are silent — the lookup returns a result, the estimator returns a number, neither errors.

The two call for different responses, which is worth keeping distinct. Ascertainment is a
tooling problem: accession-aware matching fixes it. The resolution limit is not a tooling
problem — no method produces benign controls that do not exist. One is addressed by changing
the pipeline, the other only by generating data.

---

## Open items

1. **~~Verify BRAF I208V~~ — done.** The score set declares `offset 152` and a 76-residue
   construct; canonical 208 maps to construct position 56, which the declared sequence shows
   as isoleucine, and the record reads `p.Ile56Val`. The attribution is confirmed. The
   consequence, however, is that the assay measures domain stability rather than signalling —
   see Table 3.

2. **Extend `NP_GENE`** if the audit is rerun on a wider panel.

3. **Add a mechanism column to the VCEP table.** `assay_available` alone is misleading when
   every available assay is mechanistically orthogonal to the disease. The audit produces the
   attribution; the mechanism judgement is manual and belongs beside it.

4. **Decide how `in combination only` is reported onward.** It is currently a distinct value
   in `measurement_context`. Collapsing it into either "measured" or "no assay" would
   misrepresent the record, and the distinction should survive into anything sent to a VCEP.
