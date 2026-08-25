---
name: r-package-release
description: Prepare and publish a versioned R package release by auditing changes since the previous release, drafting NEWS for approval, running release checks, and creating the corresponding Git tag and GitHub release. Use for explicit release, versioning, tagging, or publishing requests; not for ordinary development changes.
---

# R Package Release

Use this skill for a deliberate release workflow. A release changes project
history and may publish an external artifact, so separate preparation from
execution.

## Phase 1: prepare the changelog for approval

- Read the repository's `AGENTS.md` and release-related project files.
- Confirm the target version, release branch, and previous release tag. Do not
  guess a version from the current branch name.
- Require a clean working tree before starting, or clearly separate and
  preserve unrelated user changes.
- Audit every commit after the previous release tag. Categorize user-facing
  features, fixes, compatibility changes, documentation, dependency/runtime
  changes, and release-only maintenance. Do not write release notes from the
  most recent feature alone.
- Draft a new version section at the top of `NEWS.md`, based on the complete
  commit audit. Preserve existing release sections verbatim. Distinguish
  user-visible changes from internal policy, test cleanup, and development-only
  files; mention dependency/runtime upgrades when they can affect users.
- Report the commits and user-facing changes included in the draft, and the
  changes intentionally omitted.

Stop after this phase and ask the user to approve or revise the changelog.
Do not change `DESCRIPTION`, commit, push, tag, create a GitHub release, or
submit to CRAN/Bioconductor before approval. The draft may remain as an
uncommitted working-tree change for review.

`NEWS.md` is release history: never put unreleased work into an existing
version's section. If the user rejects the draft, revise only the proposed
section and ask again.

## Phase 2: execute the approved release

Proceed only after explicit approval of the changelog and target version.

- Update the package version in `DESCRIPTION` to the approved version.
- Run the repository's strongest practical local release checks. For an R
  package, normally include the full testthat suite, JavaScript tests and
  generated-artifact checks when applicable, `R CMD build .`, and
  `R CMD check --as-cran` on the resulting source tarball.
- If checks need network access or external credentials, report that boundary
  instead of silently weakening release verification. Inspect the built
  tarball for accidental development files, generated caches, or oversized
  artifacts when the package includes bundled data or runtime assets.
- Commit the approved version and NEWS changes with a release-specific message.
- Push the release commit to the intended release branch and verify the remote
  commit.
- Create an annotated `vX.Y.Z` tag on that exact commit and push the tag.
- Create and publish the GitHub Release for that tag, using the approved NEWS
  section as its notes. Prefer the available GitHub MCP connector for GitHub
  operations; otherwise provide the exact UI or CLI fallback.
- Verify the release/tag and any release-triggered CI or documentation
  deployment. Report failures without claiming the release completed.

Do not merge, delete branches, or submit to CRAN/Bioconductor unless the user
explicitly requests those additional actions. A GitHub/package release and a
CRAN or Bioconductor submission are separate workflows.
