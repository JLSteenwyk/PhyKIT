#!/usr/bin/env python3
"""Publish an already-tested recipe from a fresh, credential-isolated runner."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


UPSTREAM = "bioconda/bioconda-recipes"
FILES = ("meta.yaml", "phykit_release_test.py")


def run(*args, cwd=None):
    return subprocess.check_output(args, cwd=cwd, text=True).strip()


def gh_json(*args):
    return json.loads(run("gh", *args))


def check_pull_requests(fork, branch):
    """Do not duplicate or overwrite a maintainer/autobump recipe PR."""
    pulls = gh_json(
        "pr", "list", "--repo", UPSTREAM, "--state", "open",
        "--search", "phykit in:title", "--limit", "1000",
        "--json", "number,url,headRefName,headRepositoryOwner,headRepository",
    )
    own = None
    for pull in pulls:
        repository = pull.get("headRepository") or {}
        owner = pull.get("headRepositoryOwner") or {}
        if (pull["headRefName"] == branch
                and owner.get("login") == fork.split("/")[0]
                and repository.get("name") == fork.split("/")[1]):
            own = pull
            continue
        pages = gh_json("api", "--paginate", "--slurp",
                        f"repos/{UPSTREAM}/pulls/{pull['number']}/files?per_page=100")
        if any(f["filename"].startswith("recipes/phykit/") for page in pages for f in page):
            raise ValueError(f"Another PhyKIT recipe PR is open: {pull['url']}. Resolve it, then rerun.")
    return own


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--fork", default="JLSteenwyk/bioconda-recipes")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9-]+/bioconda-recipes", args.fork):
        raise ValueError("Expected OWNER/bioconda-recipes fork")
    if not os.environ.get("GH_TOKEN"):
        raise ValueError("Configure the BIOCONDA_GITHUB_TOKEN repository secret")
    manifest = json.loads((args.artifact / "result.json").read_text())
    version = manifest["version"]
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        raise ValueError("Invalid manifest version")
    if not manifest["changed"]:
        print("Upstream recipe is already synchronized")
        return
    branch = f"phykit-release/{version}"
    own = check_pull_requests(args.fork, branch)
    fork = gh_json("api", f"repos/{args.fork}")
    if not fork.get("fork") or fork.get("parent", {}).get("full_name") != UPSTREAM:
        raise ValueError("Configured repository is not a fork of bioconda-recipes")
    checkout = args.checkout.resolve()
    run("git", "clone", "--filter=blob:none", "--sparse",
        f"https://github.com/{UPSTREAM}.git", str(checkout))
    git = lambda *values: run("git", *values, cwd=checkout)
    git("sparse-checkout", "set", "recipes/phykit")
    # Unrelated Bioconda updates are fine; changes to this recipe require a rebuild.
    if git("rev-parse", "origin/master:recipes/phykit") != manifest["recipe_tree"]:
        raise ValueError("Upstream PhyKIT recipe changed during the build; rerun to rebuild")
    git("config", "user.name", "PhyKIT release automation")
    git("config", "user.email", "41898282+github-actions[bot]@users.noreply.github.com")
    git("remote", "add", "fork", f"https://github.com/{args.fork}.git")
    remote_branch = git("ls-remote", "--heads", "fork", f"refs/heads/{branch}")
    if remote_branch:
        git("fetch", "fork", f"refs/heads/{branch}:refs/remotes/fork/{branch}")
        git("switch", "-c", branch, f"fork/{branch}")
        changed = git("diff", "--name-only", "origin/master...HEAD").splitlines()
        if set(changed) - {f"recipes/phykit/{name}" for name in FILES}:
            raise ValueError("Automation branch contains unrelated changes; refusing to overwrite")
        git("merge", "--no-edit", "origin/master")
    else:
        git("switch", "-c", branch, "origin/master")
    for name in FILES:
        source = args.artifact / "recipe" / name
        if source.is_symlink() or not source.is_file():
            raise ValueError(f"Missing regular artifact file: {name}")
        shutil.copyfile(source, checkout / "recipes" / "phykit" / name)
    git("add", *(f"recipes/phykit/{name}" for name in FILES))
    if git("diff", "--cached", "--name-only"):
        git("commit", "-m", f"Update phykit {version} release metadata and package checks")
    # Check again immediately before writing, including autobump PRs opened during setup.
    own = check_pull_requests(args.fork, branch)
    run("gh", "auth", "setup-git")
    git("push", "fork", f"HEAD:refs/heads/{branch}")
    body = (
        f"Synchronize the published PyPI PhyKIT {version} release.\n\n"
        f"- Verify wheel and source SHA256 checksums.\n"
        f"- Synchronize {manifest['entry_points']} console scripts, Python constraints, "
        "and explicitly mapped runtime dependencies.\n"
        "- Preserve upstream recipe tests and add installed-command, import, tree-length, "
        "and headless plotting checks.\n"
        "- Build and test the generated Conda package before opening this PR.\n\n"
        f"Source SHA256: `{manifest['source_sha256']}`\n\n"
        "Generated by PhyKIT release automation. Bioconda review and publication are still required."
    )
    if own:
        run("gh", "pr", "edit", str(own["number"]), "--repo", UPSTREAM, "--body", body)
        url = own["url"]
    else:
        url = run("gh", "pr", "create", "--repo", UPSTREAM, "--base", "master",
                  "--head", f"{args.fork.split('/')[0]}:{branch}",
                  "--title", f"Update phykit {version} release metadata", "--body", body)
    print(url)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as summary:
            summary.write(f"Bioconda recipe pull request: {url}\n")


if __name__ == "__main__":
    main()
