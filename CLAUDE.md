# DRAKKAR — Agent instructions

## Release process

When the user asks to release a version (e.g. "release 1.7.2"), follow every
step below in order. Do not skip any step.

### 1. Update CHANGELOG.md

Move the pending items from `## [Unreleased]` into a new versioned section
directly above the previous release, using today's date:

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added / Changed / Fixed
- …
```

Reset the Unreleased block to:

```markdown
## [Unreleased]

### Added

- No unreleased changes yet.
```

If there are no unreleased items, write a minimal entry that describes what
the release contains (e.g. documentation or maintenance changes).

### 2. Bump the version in two files

- `pyproject.toml` — the `version = "…"` field under `[project]`
- `drakkar/__init__.py` — the `__version__ = "…"` line

Both must match the new version exactly.

### 3. Commit

Stage `CHANGELOG.md`, `pyproject.toml`, `drakkar/__init__.py`, and any other
files changed as part of the release (e.g. `drakkar/workflow/config.yaml`).
Commit with the message:

```
Release vX.Y.Z
```

### 4. Tag

Create an annotated or lightweight tag on the release commit:

```bash
git tag vX.Y.Z
```

### 5. Push branch and tag

Both must be pushed for the GitHub Actions release workflow to fire:

```bash
git push origin main
git push origin vX.Y.Z
```

Pushing the tag triggers `.github/workflows/release.yml`, which builds the
package, creates the GitHub release with auto-generated notes, and publishes
to PyPI. **If the tag is not pushed, no GitHub release is created.**

### Common mistakes to avoid

- Committing without pushing the tag — the GitHub release will not be created.
- Pushing the tag before pushing the branch — the tag may point to a commit
  that is not yet on the remote branch.
- Forgetting to bump one of the two version files.
- Leaving the `[Unreleased]` section non-empty after cutting a release.
