# MutGlyph development policy

MutGlyph provides interactive, GenomeSpy-based counterparts to established
cancer-genomics plots. Keep implementations small and composable, and prefer
semantic compatibility over pixel-level reproduction of another package.

## Design principles

Follow KISS. Prefer the smallest implementation that clearly supports the
common use case. Do not introduce abstractions, configuration layers,
compatibility shims, or dependencies until a concrete requirement justifies
them. Remove accidental complexity when extending existing code rather than
building around it.

Keep features composable. Data preparation, annotation, ordering, and
GenomeSpy specification construction should remain separable when that enables
reuse or testing. Prefer ordinary data frames and small, inspectable
intermediate structures over opaque objects. Make it straightforward to supply
custom mutation data, annotations, and domains where those are natural inputs.

Composability does not mean making every internal detail configurable. Keep the
common public API concise, provide sensible defaults, and expose extension
points only for demonstrated user needs. Do not build frameworks, inheritance
hierarchies, or plug-in systems for hypothetical future features.

Re-evaluate the design after each meaningful implementation step. New evidence
may reveal flaws in the original plan; simplify or explore an alternative when
necessary. If further progress requires a fundamental design decision, pause
the implementation and resolve that decision explicitly before adding more
code.

## Compatibility and public API

When MutGlyph intentionally provides an interactive counterpart to an
established plot, use the established function name and exact casing. Current
and anticipated examples include:

- `MutGlyph::oncoplot()`, corresponding to `maftools::oncoplot()`
- `MutGlyph::rainfallPlot()`, corresponding to `maftools::rainfallPlot()`
- `MutGlyph::lollipopPlot()`, corresponding to `maftools::lollipopPlot()`
- `MutGlyph::cnFreq()` and `MutGlyph::cnSpec()` if direct counterparts to the
  GenVisR functions are implemented

Use an established name only when the function is genuinely the same kind of
plot. MutGlyph-specific helpers should have descriptive, unambiguous names.

Aim for compatibility with common input types, argument names, meanings, and
semantic defaults. Representative calls from the source package should often
work after changing only the namespace. Add compatibility tests for such
calls. Do not reproduce obscure styling, output, or layout arguments when they
do not make sense for GenomeSpy. Unsupported arguments must not be silently
ignored: omit them from the API or report them clearly.

Returning a GenomeSpy htmlwidget is an intentional difference. Interactive
behavior, responsive composition, and a useful generated specification take
priority over exact visual or return-value parity.

Avoid broad imports and re-exports from maftools, GenVisR, or other source
packages. Use explicit namespace-qualified calls internally. In documentation
and examples where both packages are relevant, likewise prefer explicit calls
such as `maftools::oncoplot()` and `MutGlyph::oncoplot()` instead of relying on
package masking.

### Rationale

Familiar names and semantics make MutGlyph easy to discover and let users try
an interactive replacement with minimal changes to existing analyses. Limiting
compatibility to meaningful behavior keeps the API honest and avoids carrying
implementation-specific complexity from static plotting systems into
GenomeSpy.

## Code and idea provenance

Before copying or adapting code, verify that its license is compatible with
MutGlyph and comply with its attribution and notice requirements. Identify the
source near each adapted block with a concise comment, including enough detail
to locate the original implementation (for example, package, file or function,
and a stable URL, version, or commit when available). Preserve required
copyright and license notices in `inst/NOTICE` or another appropriate bundled
notice file. Do not describe independently written code as copied merely
because it implements the same documented behavior.

When a design or algorithm is materially informed by another package but no
code is copied, mention that influence in a nearby source comment. Also mention
it in user-facing documentation when it helps users understand compatibility,
behavior, or credit. Avoid repetitive attribution in every plot or article;
one clear contextual acknowledgement is normally enough.

## Data provenance and reproducibility

Before adding data, verify that its license, access conditions, and privacy
constraints permit redistribution in the intended form. Retain enough
provenance to explain and reproduce every packaged dataset:

- authoritative source URLs and dataset or accession identifiers
- retrieval or preparation date and, when useful, source versions, manifests,
  queries, checksums, or workflow identifiers
- all wrangling scripts and material manual decisions
- documentation of filtering, aggregation, coordinate systems, and output
  columns
- required citations and acknowledgements

Keep reproducible preparation workflows under `data-raw/` and commit them to
Git together with small manifests or metadata needed to freeze the inputs.
Commit generated package datasets under `data/`. Exclude `data-raw/` and large
download caches from the built package; do not commit downloaded source files
when a manifest and script can retrieve and verify them. Vignettes and examples
must build without network access from the packaged, generated data.

Record dataset provenance in the dataset documentation and use `inst/NOTICE`
for license text or package-level attribution when appropriate. Put a visible
acknowledgement next to an example only when the source requests it or the
context would otherwise be unclear; avoid repeating the same disclaimer on
every documentation page.
