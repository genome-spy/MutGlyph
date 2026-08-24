# Gene Annotation Tracks Implementation Plan

## Goal

Add optional, generic genomic annotation tracks to `rainfallPlot()` and
`gisticChromPlot()`. The first available track will be a scored, gene-body
track for the supported human assemblies. The track will provide:

- familiar gene symbols as navigation landmarks at broad genomic scales;
- all gene bodies as local context when the user zooms in;
- strand direction through arrow-shaped gene bodies; and
- compact widget payloads suitable for interactive RMarkdown documents.

The initial gene representation intentionally excludes exons, transcripts, and
exon-union encodings. A gene is one genomic interval with a symbol, stable
identifier, strand, and label-priority score.

The implementation must be developed on a feature branch with one coherent,
working commit per milestone:

```text
codex/gene-annotations
```

Do not create the branch or begin implementation until the plan is approved.

## Scope

Implement:

1. an R-based, reproducible preparation and scoring workflow;
2. compact assembly-specific scored gene-body annotations;
3. `GRanges` and ordinary data-frame input support;
4. a generic `annotationTracks` argument for both plots;
5. a GenomeSpy track using scored labels and `arrow-block` gene bodies;
6. synchronized rainfall and GISTIC annotation views; and
7. tests, provenance, documentation, and offline examples.

The public plotting API should be generic even though the first built-in
annotation source is scored genes:

```r
rainfallPlot(
  maf,
  ref.build = "hg19",
  annotationTracks = list(
    genes = mutglyph_gene_annotations("hg19")
  )
)
```

`gisticChromPlot()` should use the same `annotationTracks` contract. Its
existing `annotations` argument remains separate and continues to represent
labels attached to GISTIC score intervals.

Add `annotationTracks` at the end of each existing public function signature
so existing positional calls remain valid. With multiple tracks, preserve the
named-list order, give each track a compact title, and allocate a documented
height step per lane. Do not impose a silent global lane cap that could hide
gene bodies; any future cap must be explicit and tested.

Use `height = list(step = 18)` for annotation views, so the view grows by
18 pixels per occupied lane and never hides bodies through an implicit lane
limit.

## Explicitly deferred

Do not add these in the first implementation:

- exon or transcript rendering;
- representative-transcript selection through MANE, RefSeq Select, or
  GENCODE Primary;
- automatic network downloads inside a plotting function;
- a general R implementation of the GenomeSpy grammar;
- user-defined mark specifications or arbitrary per-track layout objects;
- custom JavaScript interaction beyond normal GenomeSpy behavior; or
- citation counts presented as biological importance scores.

Representative transcript sets may become relevant if MutGlyph later adds
transcript or exon tracks. They are unnecessary for one-range-per-gene bodies.

## Design decisions

### Public track input

The public argument is:

```r
annotationTracks = NULL
```

When non-`NULL`, it is a named list of annotation tables. Each table may be a
`GRanges` object or an ordinary data frame.
Names must be unique and non-empty; unnamed or duplicate tracks are rejected.

The minimum normalized contract is:

```text
chromosome / seqnames
start
end
score
```

Optional fields are:

```text
label or symbol
strand
identifier or gene_id
```

Track names come from the list names. The first renderer supports ranged
features, optional labels, optional strand direction, and scored labels. It is
a small interval-track contract, not a general visualization grammar. A finite
numeric `score` is required for every row in every track;
larger values have higher label priority. `symbol` and `gene_id` are
accepted aliases for `label` and `identifier`, respectively, and are
normalized at the boundary.

The canonical data-frame convention is 1-based, closed intervals, matching
`GRanges`. Callers supplying BED-style 0-based, half-open coordinates must
convert them before calling MutGlyph. Input aliases such as `chromosome` and
`seqnames`, and `1`/`chr1`, may be accepted at this boundary; normalized
internal data use the package's `chr`-prefixed sequence names. For human
assemblies, the initial canonical set is `chr1`--`chr22`, `chrX`, `chrY`, and
`chrM`; alternate, unlocalized, and unplaced contigs are excluded unless a
future track contract explicitly supports them.

Missing optional fields have deterministic fallbacks: absent labels produce no
label layer, and absent strand renders a directionless interval with a
rectangle mark. A missing, non-numeric, or non-finite score is an error.
Within a track, strand must be complete or absent; mixed missing and
non-missing strand values are rejected; `*` is treated as absent. Track order
follows the named-list order from top to bottom. Per-feature colors are
intentionally not in the
initial contract because their rendering semantics are not yet defined.

### Gene annotation fields

The built-in gene table must contain at least:

```text
seqnames
start
end
strand
symbol
gene_id
score
```

Here `seqnames` means the chromosome or genomic contig identifier, not the
nucleotide sequence itself. The score is used only by
`filterScoredLabels` to choose labels when space is limited.

### GenomeSpy rendering

The reusable renderer should:

- derive `start0` and `end0` with
  `linearizeGenomicCoordinate`, then `collect`/sort deterministically
  before `pileup` (which requires sorted intervals);
- use `measureText` and `filterScoredLabels` for collision-aware labels,
  passing the linearized start, end, midpoint, measured width, lane, and the
  supplied score;
- render every feature body at close zoom;
- render a fully stranded track with one `arrow` mark using
  `style = "arrow-block"`, mapping `+` to `forward` and `-` to
  `reverse`;
- render a directionless track with one `rect` mark;
- use no separate point mark for strand;
- use one opacity transition, `unitsPerPixel = c(100000, 40000)` and
  `values = c(0, 1)`, to subdue bodies at broad scales; and
- include available label, stable identifier, R-facing coordinates, strand, and
  score in tooltips.

The score must not control body color or imply biological importance. The
score-based label-priority idea is inspired by the HiGlass gene annotation
track workflow, which uses an importance column to determine which labels
remain visible at lower zoom levels:

<https://docs.higlass.io/data_preparation.html#gene-annotation-tracks>

MutGlyph will implement the workflow independently in R using NCBI
Gene2PubMed-derived citation counts as a navigation heuristic. No HiGlass code
will be copied. If an implementation adapts an existing GenomeSpy example or
specification block, record the source, pinned version or commit, and MIT
attribution beside the adapted block and in `inst/NOTICE` when required.

The upstream sort order must be deterministic so equal scores have stable
results. This plan does not promise a separate search interaction; label
selection is the initial use of the score.

### Data preparation and provenance

The reproducible workflow belongs under:

```text
data-raw/gene-annotations/
```

It should contain R scripts, a README, and a manifest. The workflow should:

- download assembly-specific gene annotations and NCBI mapping tables;
- parse gene-level features;
- retain supported canonical sequences;
- calculate citation counts from `gene2pubmed`;
- join scores through stable NCBI GeneIDs;
- retain zero-score genes;
- produce one range per gene;
- write compact generated artifacts; and
- record source URLs, releases, retrieval dates, checksums, and preparation
  decisions.

Freeze the assembly-matched UCSC RefSeq `refGene` table for hg18, hg19, and
hg38, together with the human NCBI `gene2refseq`, `gene_info`, and
`gene2pubmed` snapshots. Strip RefSeq version suffixes before joining. Retain
only `chr1`--`chr22`, `chrX`, `chrY`, and `chrM`. For each NCBI GeneID,
join `refGene.name` to the RefSeq accession column in human `gene2refseq`
to obtain GeneID, then join GeneID to `gene_info` for the symbol and to
`gene2pubmed` for PubMed IDs. The manifest must record these source-column
joins because the downloaded tables are positional rather than self-describing.
combine all retained transcript records into one inclusive gene-body range
using the minimum transcript start and maximum transcript end, and require a
single consistent strand. Convert UCSC's 0-based half-open `txStart`/`txEnd`
to the package's 1-based closed `start`/`end` fields before writing the
artifact. Retain distinct GeneIDs that share a symbol.
Calculate scores from unique `(GeneID, PubMed_ID)` pairs, retain genes with no
citations as score zero, and use the NCBI GeneID as the stable identifier.
This deliberately produces one body per gene without choosing a
representative transcript or rendering exons.

The manifest must record exact URLs, source releases or snapshot dates,
retrieval dates, checksums, taxon, chromosome filter, row counts at each
stage, and the range-collapse and score-aggregation decisions. No alternate
source may be substituted for an older assembly without naming and freezing a
compatible release.

Raw downloads and large caches must not be committed or included in the built
package.

Ship the reduced, gene-body-only tables inside `MutGlyph` for the initial
implementation. This is the simplest way to keep
`mutglyph_gene_annotations()` offline and avoids a hidden companion-package
or AnnotationHub runtime dependency. Record installed package size, serialized
resource size, and one-track widget JSON size for each assembly. Treat a
regression above 25% or total added installed-data size above 10 MB as a
release-blocking review point rather than silently moving the data to a network
resource. Plot functions and the loader must never download implicitly.

### `GRanges` support

`GRanges` is the ecosystem-facing representation for genomic intervals. It
should be supported directly when the optional Bioconductor dependencies are
installed, while ordinary data frames remain accepted without those
dependencies for lightweight custom use and tests.

Declare `GenomicRanges` and `GenomeInfoDb` in `Suggests`, not
`Imports`. Data-frame plotting must work without them. The built-in loader and
`GRanges` input path should give a clear installation error when they are
absent. The plotting path should not require `rtracklayer` or access to raw
GTF files; those belong only to data preparation. Recognize assembly aliases
`hg18`/`NCBI36`, `hg19`/`GRCh37`, and `hg38`/`GRCh38`.

Use one internal assembly canonicalizer for both plot arguments and
`GenomeInfoDb::genome()` values before comparing or storing assemblies. The
canonicalizer accepts the aliases above and rejects unknown values.

All R-facing interval inputs use 1-based, closed coordinates, matching
`GRanges`. At the normalization boundary, derive 0-based, half-open internal
coordinates for `linearizeGenomicCoordinate` and `pileup`
(`start0 = start - 1`, `end0 = end`). Retain the R-facing coordinates for
locus channels and tooltips so they align with the existing rainfall and
GISTIC plots.

Read assembly metadata only from `GenomeInfoDb::genome(gr)`. Ignore
`metadata(gr)` for assembly detection so there is one authoritative source.
If the result contains no non-missing value, the plot's `ref.build` is the
caller-supplied assertion. If it contains more than one distinct value, or an
unknown value, error. If its canonical value conflicts with `ref.build`, error
rather than guessing; no coordinate conversion between assemblies is
attempted. Data frames have no inferred assembly and are
interpreted in the plot's `ref.build`; they cannot request an alternate
assembly through the track API. Sequence aliases are normalized only when
unambiguous, `*` strand is normalized to absent, and `chrM` is the sole
mitochondrial spelling in the normalized schema.

Tests must cover one-base intervals, inclusive ends, sequence aliases,
mitochondrial naming, missing metadata, conflicting metadata, and assembly
mismatches.

Do not introduce a custom S4 annotation class in the first implementation.

## Working and re-evaluation loop

After every milestone:

1. run the smallest relevant R tests;
2. inspect the generated specification and data payload;
3. remove unused fields and premature abstractions;
4. manually render the changed track when layout or scale behavior changes;
5. compare the result with the remaining plan; and
6. commit only a coherent working increment.

Pause before proceeding if implementation requires:

- a second public track grammar;
- automatic network access during plotting;
- exon/transcript data to make the gene track work;
- custom JavaScript interaction; or
- a different coordinate or assembly contract.

## Milestone 0 — Branch and plan

Create the feature branch from the current clean `main` HEAD and commit this
plan.

```text
git switch -c codex/gene-annotations
```

Commit:

```text
docs: plan gene annotation tracks
```

Verification:

- branch name is `codex/gene-annotations`;
- the worktree is clean before implementation; and
- the plan is present under `plans/`.

## Milestone 1 — R annotation preparation workflow

Add the reproducible R scripts under `data-raw/gene-annotations/`.

Use the selected `refGene` records only to derive one bounding gene-body
interval per GeneID. Do not port the MCCA exon-union logic, because the first
MutGlyph track does not render exons.

Add tests using small local `refGene`/mapping fixtures for:

- gene parsing;
- canonical sequence filtering;
- stable-ID and symbol retention;
- citation-count aggregation;
- zero-score handling;
- strand preservation; and
- provenance manifest generation.

All tests must work without network access.

Commit:

```text
feat: add reproducible gene annotation preparation
```

## Milestone 2 — Packaged annotation data and loader

Add the compact generated artifacts for the supported assemblies and a public
loader. Add `GenomicRanges` and `GenomeInfoDb` to `Suggests`, export the
loader, and generate its NAMESPACE and documentation entries:

```r
mutglyph_gene_annotations(ref.build)
```

The loader should return `GRanges` with the documented metadata columns and a
stable assembly/source description.

This milestone must also add packaged provenance metadata. The loader returns a
`GRanges` when `GenomicRanges` is installed and otherwise reports the
installation requirement. This is intentional: ordinary data-frame tracks
remain usable without Bioconductor, while the built-in GRanges-returning loader
and `GRanges` input path require the suggested packages. It must work without
network access; no loader path may trigger a download.

Measure the installed package and serialized table sizes against the release
gates in the packaging section. Keep the resources bundled for this initial
implementation unless that gate is exceeded and the design is explicitly
revisited.

Commit:

```text
feat: add packaged scored gene-body annotations
```

Tests should verify assembly selection, field names, score preservation, and
offline loading.

## Milestone 3 — Generic annotation-track normalization

Implement one internal normalizer for `annotationTracks`.

It should:

- accept `GRanges` and data frames;
- normalize chromosome, start, end, label, strand, mandatory score, and
  identifier;
- perform explicit coordinate conversion and assembly conflict checks;
- preserve track names;
- support multiple tracks; and
- validate only the fields required by each supplied track, with score required
  for all tracks and strand either complete or absent for a whole track.

The normalized track tables must retain the supplied score for custom
annotations; the renderer must use it in `filterScoredLabels` rather than
replacing it with a gene- or GISTIC-derived value. Directionless custom
intervals use the rectangle renderer and still use their scores for labels.

Test absent labels, missing strand, missing or invalid scores, one-base
intervals, the 1-based closed input convention, assembly mismatch errors, and
custom score preservation.

Commit:

```text
feat: add generic annotation track inputs
```

Add structural tests for single and multiple tracks, including directionless
custom intervals.

## Milestone 4 — GenomeSpy gene-body renderer

Add the reusable internal GenomeSpy view for annotation tracks.

Implement the renderer invariants described above: sorted linearized
intervals, score-driven labels, `arrow-block` for complete strand,
`rect` for directionless tracks, and no point, exon, or transcript layers.

Commit:

```text
feat: render scored gene-body annotation tracks
```

Specification tests must assert the arrow mark, style, direction mapping,
scored labels (including a custom non-gene score), the rectangle fallback, and
absence of point/exon layers. Labels must use measured text width, linearized
interval position/end position, and lane in `filterScoredLabels`.

Rendering tests must also cover directionless custom intervals, deterministic
ordering, close-zoom visibility of all bodies, and the documented label-fade
or label-priority behavior at broad scales.

## Milestone 5 — Rainfall integration

Add `annotationTracks` to `rainfallPlot()` as the final argument, preserving
all existing positional arguments.

Place annotation tracks below the rainfall panel using a shared genomic x-scale
and independent y-scales. Keep the existing default specification unchanged
when `annotationTracks = NULL`.

When tracks are present, use a conditional vertical composition with one
annotation view per named track, compact lane-based heights, a shared x scale,
and independent y scales. Keep the current single-view specification path
untouched for the `NULL` case.

Verify that:

- panning and zooming remain synchronized;
- kataegis intervals remain visually prominent;
- gene bodies do not enter the rainfall distance axis; and
- widget data contain only the compact gene-body fields.

Test that `annotationTracks = NULL` preserves the current specification and
that the new argument is appended to the public signature.

Commit:

```text
feat: add annotation tracks to rainfall plots
```

Add a vignette example showing a kataegis region with gene-body context.

## Milestone 6 — GISTIC integration

Add `annotationTracks` to `gisticChromPlot()` as the final argument,
preserving all existing positional arguments.

Place annotation tracks beneath the complete GISTIC profile and chromosome
context. Preserve the existing `annotations` argument for selected score-label
annotations; do not merge its data contract with annotation tracks.

The annotation views use the existing root shared genomic x-scale and
independent y-scales. A track's supplied score is used only for its labels; it
must never be replaced by the GISTIC interval score or used as body color.

Verify that:

- the annotation track shares the genomic x-scale;
- whole-genome views remain readable;
- local zoom shows relevant gene bodies and labels; and
- existing GISTIC output is unchanged when no tracks are supplied.

Test that `annotationTracks = NULL` preserves the current specification and
that the new argument is appended to the public signature.

Commit:

```text
feat: add annotation tracks to gistic plots
```

## Milestone 7 — Documentation and release verification

Document:

- the generic `annotationTracks` contract;
- `mutglyph_gene_annotations()`;
- supported assemblies and annotation releases;
- score provenance and HiGlass inspiration;
- the fact that scores are label-priority heuristics;
- `GRanges` and data-frame usage;
- offline behavior; and
- the intentional absence of exon models.

Document that annotation-track scores are mandatory, that larger finite values
have higher label priority, that custom scores are honored in both plots, that
unstranded tracks render as rectangles, and that conflicting `GRanges`
assembly metadata is rejected.

Include one concise HiGlass acknowledgement in the user-facing documentation:

> MutGlyph’s score-based label-priority approach is inspired by the
> [HiGlass gene annotation track workflow](https://docs.higlass.io/data_preparation.html#gene-annotation-tracks).
> MutGlyph’s R implementation and data preparation are independent.

Run:

- the complete `testthat` suite;
- `R CMD check`;
- GenomeSpy specification validation;
- package-size and widget-payload checks; and
- manual browser checks at whole-genome, chromosome, and local-locus scales.

The payload checks must inspect both the compact R-to-JavaScript transport and
the decoded GenomeSpy specification, and verify that no raw source tables,
unused metadata, or per-sample fields enter the widget.

Commit:

```text
docs: document gene annotation tracks
```

## Definition of done

The feature is complete when users can supply named `GRanges` or data-frame
tracks to both plots, use the built-in scored gene-body resource for each
supported assembly, see score-prioritized symbols at broad scales, see all
nearby arrow-block gene bodies at close scales, and build examples and package
checks without network access.
