# MutGlyph development policy

MutGlyph provides interactive, GenomeSpy-based counterparts to
established cancer-genomics plots. Keep implementations small and
composable, and prefer semantic compatibility over pixel-level
reproduction of another package.

## Design principles

Follow KISS. Prefer the smallest implementation that clearly supports
the common use case. Do not introduce abstractions, configuration
layers, compatibility shims, or dependencies until a concrete
requirement justifies them. Remove accidental complexity when extending
existing code rather than building around it.

Prioritize issues demonstrated by supported inputs, upstream behavior,
existing examples, or a plausible user workflow. Do not add validation,
branches, configuration, or abstractions solely for contrived states
with no realistic path into the package. Fix correctness problems that
affect the documented contract, but otherwise record or defer
theoretical edge cases until evidence makes them concrete.

Keep features composable. Data preparation, annotation, ordering, and
GenomeSpy specification construction should remain separable when that
enables reuse or testing. Prefer ordinary data frames and small,
inspectable intermediate structures over opaque objects. Make it
straightforward to supply custom mutation data, annotations, and domains
where those are natural inputs.

Add a brief rationale comment beside non-obvious heuristics,
workarounds, optimizations, compatibility behavior, and safety
boundaries. Explain why the code exists, the assumption it relies on, or
why the simpler-looking approach is unsuitable. Do not merely narrate
the syntax, and do not add comments to obvious code. When a heuristic
can be replaced by a small deterministic solution, prefer fixing it over
documenting it.

Composability does not mean making every internal detail configurable.
Keep the common public API concise, provide sensible defaults, and
expose extension points only for demonstrated user needs. Do not build
frameworks, inheritance hierarchies, or plug-in systems for hypothetical
future features.

Re-evaluate the design after each meaningful implementation step. New
evidence may reveal flaws in the original plan; simplify or explore an
alternative when necessary. If further progress requires a fundamental
design decision, pause the implementation and resolve that decision
explicitly before adding more code.

## Testing

Tests should protect public behavior and documented invariants, not
restate the implementation. Tests that assert exact nested-spec layout,
list positions, internal helper names, transform order, or incidental
styling without those details being part of the contract are effectively
useless: they duplicate the implementation and reject valid refactors
without finding user-facing bugs.

Prefer tests through exported functions and semantic assertions. For
generated GenomeSpy specifications, check user-visible composition such
as which views exist, how axes and scales are shared, and whether
annotations are placed in the intended relationship to the main plot.
Locate relevant components by meaning rather than fixed list positions
when practical. Keep focused unit tests for pure normalization and
preparation rules, and use integration tests for the final widget
behavior. Do not use string searches in a bundled runtime as a
substitute for exercising the behavior it is meant to provide.

## Release notes

`NEWS.md` is release history. Add entries under the release they belong
to; do not insert unreleased feature work into an existing version’s
section. Keep unreleased changes in the pull request description or a
separately marked unreleased section, and update `NEWS.md` when
preparing the corresponding release.

## User-facing documentation

Keep user-facing documentation focused on information needed to use the
package, interpret its output, or understand its current public
behavior. A detail is not documentation-worthy merely because it was
discovered or considered during implementation; include it only when it
is relevant to a plausible user decision or workflow.

State visible behavior directly. Do not explain renderer mechanics,
internal layering, serialization choices, workarounds, or other
implementation rationale unless users need the information to use the
feature correctly. Do not explain historical behavior or justify current
behavior by comparing it with a former implementation or appearance.
Historical changes belong in release notes, not current guides or
function reference pages.

## Compatibility and public API

### Naming and semantics

Familiar names and semantics make MutGlyph easy to discover and let
users try an interactive replacement with minimal changes to existing
analyses. Limiting compatibility to meaningful behavior keeps the API
honest and avoids carrying implementation-specific complexity from
static plotting systems into GenomeSpy.

When MutGlyph intentionally provides an interactive counterpart to an
established plot, use the established function name and exact casing.
Current and anticipated examples include:

- [`MutGlyph::oncoplot()`](https://genomespy.app/MutGlyph/reference/oncoplot.md),
  corresponding to
  [`maftools::oncoplot()`](https://rdrr.io/pkg/maftools/man/oncoplot.html)
- [`MutGlyph::rainfallPlot()`](https://genomespy.app/MutGlyph/reference/rainfallPlot.md),
  corresponding to
  [`maftools::rainfallPlot()`](https://rdrr.io/pkg/maftools/man/rainfallPlot.html)
- [`MutGlyph::lollipopPlot()`](https://genomespy.app/MutGlyph/reference/lollipopPlot.md),
  corresponding to
  [`maftools::lollipopPlot()`](https://rdrr.io/pkg/maftools/man/lollipopPlot.html)

Use an established name only when the function is genuinely the same
kind of plot. MutGlyph-specific helpers should have descriptive,
unambiguous names.

Aim for compatibility with common input types, argument names, meanings,
and semantic defaults. Representative calls from the source package
should often work after changing only the namespace. Add compatibility
tests for such calls. Do not reproduce obscure styling, output, or
layout arguments when they do not make sense for GenomeSpy. Unsupported
arguments must not be silently ignored: omit them from the API or report
them clearly.

Returning a GenomeSpy htmlwidget is an intentional difference.
Interactive behavior, responsive composition, and a useful generated
specification take priority over exact visual or return-value parity.

### Upstream reuse

Checking upstream capabilities first avoids parallel implementations
that can silently diverge as maftools evolves. Reusing stable public
behavior reduces MutGlyph’s correctness and maintenance burden; owning a
local implementation is justified only when it provides a concrete
benefit that upstream reuse cannot.

Before implementing behavior that overlaps with maftools, inspect its
current public API, documentation, and relevant implementation to
determine whether it already provides the needed parsing, normalization,
annotation, aggregation, or other functionality. Prefer explicit,
namespace-qualified calls to stable public maftools functions and
accessors when they meet MutGlyph’s needs. Do not reimplement upstream
functionality merely because it appears small: duplicate code adds
correctness, compatibility, testing, and maintenance burdens.

Implement equivalent behavior locally only when using maftools directly
would conflict with MutGlyph’s interactive or composable design, package
constraints, performance needs, or documented correctness requirements.
Keep such local implementations minimal, record the reason in a nearby
rationale comment, and add compatibility or parity tests for the
supported behavior. Avoid depending on undocumented maftools internals
as a shortcut; if no suitable public API exists, prefer a small owned
implementation with an explicit contract.

### Namespace discipline

Avoid broad imports and re-exports from maftools or other source
packages. Use explicit namespace-qualified calls internally. In
documentation and examples where both packages are relevant, likewise
prefer explicit calls such as
[`maftools::oncoplot()`](https://rdrr.io/pkg/maftools/man/oncoplot.html)
and
[`MutGlyph::oncoplot()`](https://genomespy.app/MutGlyph/reference/oncoplot.md)
instead of relying on package masking.

## Code and idea provenance

Before copying or adapting code, verify that its license is compatible
with MutGlyph and comply with its attribution and notice requirements.
Identify the source near each adapted block with a concise comment,
including enough detail to locate the original implementation (for
example, package, file or function, and a stable URL, version, or commit
when available). Preserve required copyright and license notices in
`inst/NOTICE` or another appropriate bundled notice file. Do not
describe independently written code as copied merely because it
implements the same documented behavior.

When a design or algorithm is materially informed by another package but
no code is copied, mention that influence in a nearby source comment.
Also mention it in user-facing documentation when it helps users
understand compatibility, behavior, or credit. Avoid repetitive
attribution in every plot or article; one clear contextual
acknowledgement is normally enough.

## Data provenance and reproducibility

Before adding data, verify that its license, access conditions, and
privacy constraints permit redistribution in the intended form. Retain
enough provenance to explain and reproduce every packaged dataset:

- authoritative source URLs and dataset or accession identifiers
- retrieval or preparation date and, when useful, source versions,
  manifests, queries, checksums, or workflow identifiers
- all wrangling scripts and material manual decisions
- documentation of filtering, aggregation, coordinate systems, and
  output columns
- required citations and acknowledgements

Keep reproducible preparation workflows under `data-raw/` and commit
them to Git together with small manifests or metadata needed to freeze
the inputs. Commit generated package datasets under `data/`. Exclude
`data-raw/` and large download caches from the built package; do not
commit downloaded source files when a manifest and script can retrieve
and verify them. Vignettes and examples must build without network
access from the packaged, generated data.

Record dataset provenance in the dataset documentation and use
`inst/NOTICE` for license text or package-level attribution when
appropriate. Put a visible acknowledgement next to an example only when
the source requests it or the context would otherwise be unclear; avoid
repeating the same disclaimer on every documentation page.
