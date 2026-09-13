# Bioconda Release Automation

The `Sync Bioconda release` workflow prepares a Bioconda recipe PR after a
stable GitHub release is published. It does not publish to PyPI, merge upstream
PRs, or upload Conda packages. Bioconda controls review and publication.

## One-Time Setup

1. Keep a fork of `bioconda/bioconda-recipes`. The default is
   `JLSteenwyk/bioconda-recipes`; override it with repository variable
   `BIOCONDA_FORK=OWNER/bioconda-recipes` if needed.
2. Add repository Actions secret `BIOCONDA_GITHUB_TOKEN` using GitHub's settings
   or `gh secret set BIOCONDA_GITHUB_TOKEN --repo JLSteenwyk/PhyKIT`. Enter the
   credential interactively; never put it in a file or command-line argument.
   The credential must be able to push branches to that fork and create/edit
   PRs in the public upstream repository. A dedicated account's classic PAT
   with `public_repo` scope supports this cross-owner fork workflow. That scope
   grants access to public repositories beyond this fork, so use a dedicated
   account, an expiry, and rotate it. Fine-grained PATs have limitations for
   contributions to repositories owned by others; do not assume a fork-only
   token can also open an upstream PR. A suitably authorized GitHub App/user
   token is another option. The built-in `GITHUB_TOKEN` is repository-scoped
   and cannot substitute for this credential.
3. Ensure this workflow and `scripts/bioconda/` exist in future release tags.
   For the manual **Run workflow** button, the workflow must also exist on the
   repository's default branch. At implementation time development is on
   `main`, but GitHub's default branch is `master`: merge the automation into
   the default branch or deliberately switch the default to `main` before
   relying on manual dispatch. The automation does not change branch settings.
4. Run a dry run for an existing PyPI release before enabling write operations.

## Release and Rerun

Publish the package to PyPI using the existing release procedure, then publish
the corresponding stable GitHub release with tag `X.Y.Z` or `vX.Y.Z`. Either
publication order works if both PyPI artifacts become available within the
retry window (30 attempts, 30 seconds between attempts, plus network time).
Drafts do not trigger publication and prereleases are skipped. Version inputs
containing prerelease, development, or local-version suffixes are rejected.

From Actions, select **Sync Bioconda release**, then **Run workflow**. Supply
the exact published version. `dry_run` defaults to `true`: the recipe is built
and tested, but nothing is pushed. To publish the PR, rerun with `dry_run=false`.

```sh
gh workflow run bioconda-release.yml --repo JLSteenwyk/PhyKIT \
  --ref main -f version=2.7.0 -f dry_run=true
```

The default-branch prerequisite above applies even when supplying `--ref main`.
If a future release workflow creates releases with its built-in `GITHUB_TOKEN`,
do not rely on the resulting release event to start another workflow: explicitly
dispatch this workflow after PyPI publication, or use an appropriate App/PAT
credential for the release event. Concurrent sync runs are serialized; rerun a
version if GitHub replaces a pending run with a newer one.

## What Gets Synchronized

- The exact release's universal wheel and source archive are downloaded from
  PyPI and checked against the SHA256 digests returned by PyPI.
- Wheel metadata supplies the version, console scripts, Python constraints, and
  runtime requirements. Nothing is imported from the development checkout.
- Dependency names use the explicit `DEPENDENCIES` mapping in `sync.py`.
  For example, PyPI `matplotlib` maps to Conda `matplotlib-base`. Unknown names,
  markers, extras, direct URLs, and unsupported version constraints fail for
  manual review rather than being dropped or guessed.
- Existing recipe metadata, build script, run exports, maintainers, and tests
  are retained. Obsolete help checks for removed entry points are removed.
  A managed `phykit_release_test.py` is added alongside existing tests.
- A new version starts at build 0; changes to the current version increment
  its upstream build number. A synchronized recipe is a no-op. Rerunning an
  unmerged automation PR uses the upstream recipe as the baseline, so it does
  not repeatedly increment the build number.
- Conda builds and tests the generated package using a Python version allowed
  by the recipe. Python 3.12 is the workflow tooling interpreter, not a pin on
  the built noarch package. Checks cover every entry point and `--help`, runtime
  imports, a known tree length, and headless PNG rendering.

## Publishing and Failures

The build job has no Bioconda write credential. A fresh publish job receives
only the generated artifact and uses the secret for GitHub operations; it
never installs or executes the generated recipe. Actions are pinned to commit
SHAs. Update those pins deliberately when maintaining the workflow.

Publishing uses `phykit-release/X.Y.Z` in the fork, normal commits, and
non-force pushes. Reruns update the same open PR. Do not manually repurpose
these automation branches. Unrelated changes on the branch are rejected.

Before pushing, the publisher checks open PRs mentioning `phykit` in their
titles and verifies their changed paths. A conflicting maintainer/autobump PR
causes a failure containing its URL; resolve that PR and rerun rather than
opening a duplicate. It will also detect the still-open packaging follow-up
PR #69188. A PR whose title never mentions PhyKIT cannot be found by this check.

If the upstream PhyKIT recipe changed during the build, publication stops and
requires a rerun. Merge conflicts and rejected pushes also fail without force
updates. Review the Actions logs and summary; generated recipes remain in the
`bioconda-recipe` artifact for 14 days, including when a build fails. Fix source
packaging or an explicit mapping and rerun the same version. No write credential
is needed for a dry run; a real publish fails clearly if the secret is absent.

## Local Verification

```sh
python -m pip install -r scripts/bioconda/requirements.txt
python -m pytest tests/unit/test_bioconda_release.py
python scripts/bioconda/sync.py --version 2.7.0 \
  --recipe /path/to/bioconda-recipes/recipes/phykit --output /tmp/sync-result.json
conda build /path/to/bioconda-recipes/recipes/phykit \
  --override-channels -c conda-forge -c bioconda --no-anaconda-upload
```

The sync command edits the supplied recipe, so use a separate checkout. Its
second run against unchanged inputs must report `"changed": false`. Tests mock
PyPI/GitHub and exercise publishing with temporary local Git repositories;
they never push to GitHub.

References: [GitHub workflow triggers](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows),
[GitHub token scope](https://docs.github.com/en/actions/concepts/security/github_token),
[Bioconda recipe guidelines](https://bioconda.github.io/contributor/guidelines.html).
