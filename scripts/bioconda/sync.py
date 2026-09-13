#!/usr/bin/env python3
"""Synchronize a Bioconda recipe from verified, published PyPI artifacts.

Does not import PhyKIT or execute code from its release artifacts.
"""

from __future__ import annotations

import argparse
import configparser
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import urlopen
from zipfile import ZipFile

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name
from packaging.version import Version
from ruamel.yaml import YAML


DEPENDENCIES = {
    "biopython": "biopython",
    "matplotlib": "matplotlib-base",
    "numpy": "numpy",
    "scipy": "scipy",
    "scikit-learn": "scikit-learn",
    "tqdm": "tqdm",
    "umap-learn": "umap-learn",
}
TEST_NAME = "phykit_release_test.py"
VERSION_LINE = re.compile(r'(?m)^{% set version = ["\']([^"\']+)["\'] %}$')


def release_version(value):
    value = value.removeprefix("v")
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", value):
        raise ValueError("Expected a stable X.Y.Z version or vX.Y.Z tag")
    return str(Version(value))


def read_url(url):
    with urlopen(url, timeout=60) as response:
        data = response.read(50 * 1024 * 1024 + 1)
    if len(data) > 50 * 1024 * 1024:
        raise ValueError(f"Download exceeds 50 MiB: {url}")
    return data


def wait_for_release(version, attempts=30, delay=30):
    if attempts < 1 or delay < 0:
        raise ValueError("Attempts must be positive and delay nonnegative")
    for attempt in range(attempts):
        try:
            release = json.loads(read_url(f"https://pypi.org/pypi/phykit/{version}/json"))
            if release["info"]["version"] != version:
                raise ValueError("PyPI returned a different version")
            files = [f for f in release["urls"] if not f.get("yanked")]
            wheels = [f for f in files if f["filename"].endswith("-py3-none-any.whl")]
            sources = [f for f in files if f["packagetype"] == "sdist"]
            if len(wheels) == len(sources) == 1:
                return wheels[0], sources[0]
        except HTTPError as error:
            if error.code not in (404, 429, 500, 502, 503, 504):
                raise
        except URLError:
            pass
        if attempt + 1 < attempts:
            print(f"Waiting for both PyPI artifacts for {version} ({attempt + 1}/{attempts})", flush=True)
            time.sleep(delay)
    raise ValueError(f"PyPI release {version} is unavailable, incomplete, or yanked")


def verified_download(artifact):
    parsed = urlparse(artifact["url"])
    if parsed.scheme != "https" or parsed.netloc != "files.pythonhosted.org":
        raise ValueError("Expected a files.pythonhosted.org HTTPS artifact")
    digest = artifact["digests"]["sha256"]
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError("Invalid SHA256 digest")
    data = read_url(artifact["url"])
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Checksum mismatch: {artifact['filename']}")
    return data


def conda_specifiers(value):
    result = []
    for specifier in sorted(SpecifierSet(value), key=str):
        if specifier.operator not in ("==", "!=", ">=", "<=", ">", "<"):
            raise ValueError(f"Unsupported Conda version constraint: {specifier}")
        if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)*(?:\.\*)?", specifier.version):
            raise ValueError(f"Unsupported Conda version: {specifier.version}")
        result.append(str(specifier))
    return ",".join(result)


def wheel_metadata(data, version):
    with ZipFile(io.BytesIO(data)) as archive:
        paths = [p for p in archive.namelist() if p.endswith(".dist-info/METADATA")]
        if len(paths) != 1:
            raise ValueError("Expected one wheel METADATA file")
        metadata = BytesParser().parsebytes(archive.read(paths[0]))
        if canonicalize_name(metadata["Name"]) != "phykit" or metadata["Version"] != version:
            raise ValueError("Wheel name/version does not match the requested release")
        config = configparser.ConfigParser(interpolation=None)
        config.optionxform = str
        config.read_string(archive.read(paths[0].replace("METADATA", "entry_points.txt")).decode())
    entries = dict(config["console_scripts"])
    if entries.get("phykit") != "phykit.phykit:main":
        raise ValueError("Missing or invalid phykit entry point")
    for name, target in entries.items():
        if not re.fullmatch(r"(?:phykit|pk_[A-Za-z0-9_]+)", name):
            raise ValueError(f"Unexpected console script name: {name}")
        if not re.fullmatch(r"phykit\.phykit:[A-Za-z_][A-Za-z0-9_]*", target):
            raise ValueError(f"Unexpected console script target: {target}")
    python = metadata.get("Requires-Python")
    if not python:
        raise ValueError("Missing Requires-Python")
    python = "python " + conda_specifiers(python)
    dependencies = []
    for value in metadata.get_all("Requires-Dist", []):
        requirement = Requirement(value)
        name = canonicalize_name(requirement.name)
        if name not in DEPENDENCIES:
            raise ValueError(f"Unknown dependency mapping: {name}; update DEPENDENCIES explicitly")
        if requirement.marker or requirement.extras or requirement.url:
            raise ValueError(f"Conditional, extra, or URL dependency needs manual mapping: {value}")
        constraint = conda_specifiers(str(requirement.specifier))
        dependencies.append((DEPENDENCIES[name] + " " + constraint).strip())
    if not dependencies:
        raise ValueError("Missing runtime dependency metadata")
    return entries, python, sorted(set(dependencies))


def parse_recipe(text):
    """Mask Jinja tokens before round-trip YAML parsing; never execute templates."""
    versions = VERSION_LINE.findall(text)
    if len(versions) != 1:
        raise ValueError("Expected one literal Jinja version assignment")
    tokens = {}

    def mask(match):
        token = f"__PHYKIT_JINJA_{len(tokens)}__"
        tokens[token] = match.group()
        return token

    headers = re.findall(r"(?m)^{%.*%}\n?", text)
    if any(not re.fullmatch(r'{% set (?:name|version) = ["\'][^"\']+["\'] %}\n?', h) for h in headers):
        raise ValueError("Unsupported Jinja statement in recipe")
    body = re.sub(r"(?m)^{%.*%}\n?", "", text).lstrip("\n")
    body = re.sub(r"{{.*?}}", mask, body)
    if "{%" in body or "{{" in body or "{#" in body:
        raise ValueError("Unsupported Jinja template in recipe")
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.width = 120
    yaml.indent(mapping=2, sequence=4, offset=2)
    return yaml, yaml.load(body), headers, tokens, versions[0]


def update_recipe(text, version, source, entries, python, dependencies, test_changed):
    yaml, recipe, headers, tokens, previous = parse_recipe(text)
    if Version(version) < Version(previous):
        raise ValueError(f"Refusing to downgrade recipe {previous} to {version}")
    if recipe["build"].get("noarch") != "python":
        raise ValueError("Expected a noarch: python recipe")
    if not isinstance(recipe["source"], dict):
        raise ValueError("Multiple recipe sources require manual review")
    # Preserve the upstream order for existing commands; append new ones deterministically.
    old_entries = recipe["build"].get("entry_points", [])
    names = [value.split(" = ")[0] for value in old_entries]
    order = list(dict.fromkeys([n for n in names if n in entries] + sorted(entries)))
    desired_entries = [f"{name} = {entries[name]}" for name in order]
    run = sorted(dependencies + [python])
    host = [python if re.match(r"^python(?:\s|$)", p) else p for p in recipe["requirements"]["host"]]
    if python not in host:
        raise ValueError("Missing host Python requirement")
    test = recipe.setdefault("test", {})
    files = list(test.get("files", []))
    if TEST_NAME not in files:
        files.append(TEST_NAME)
    commands = [c for c in test.get("commands", []) if not (
        c.endswith(" --help") and c[:-7] in names and c[:-7] not in entries
    )]
    if f"python {TEST_NAME}" not in commands:
        commands.append(f"python {TEST_NAME}")
    changed = any((
        previous != version,
        recipe["source"].get("sha256") != source["digests"]["sha256"],
        old_entries != desired_entries,
        sorted(recipe["requirements"]["run"]) != run,
        list(recipe["requirements"]["host"]) != host,
        list(test.get("files", [])) != files,
        list(test.get("commands", [])) != commands,
        test_changed,
    ))
    if not changed:
        return text
    # Retain a version-templated PyPI URL so Bioconda's autobump can still update it.
    version_token = next((key for key, value in tokens.items() if value == "{{ version }}"),
                         "__PHYKIT_VERSION__")
    tokens[version_token] = "{{ version }}"
    recipe["source"]["url"] = (
        f"https://pypi.org/packages/source/p/phykit/phykit-{version_token}.tar.gz"
    )
    recipe["source"]["sha256"] = source["digests"]["sha256"]
    build = recipe["build"]
    build["number"] = 0 if previous != version else int(build["number"]) + 1
    build["entry_points"] = desired_entries
    recipe["requirements"]["host"] = host
    recipe["requirements"]["run"] = run
    test["files"] = files
    test["commands"] = commands
    output = io.StringIO()
    yaml.dump(recipe, output)
    body = output.getvalue()
    for token, value in tokens.items():
        body = body.replace(token, value)
    header = "".join(headers)
    header = VERSION_LINE.sub(lambda m: m.group().replace(previous, version), header)
    return header.rstrip() + "\n\n" + body


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--attempts", type=int, default=30)
    parser.add_argument("--delay", type=float, default=30)
    parser.add_argument("--output", type=Path, required=True, help="JSON result manifest")
    args = parser.parse_args()
    version = release_version(args.version)
    wheel, source = wait_for_release(version, args.attempts, args.delay)
    entries, python, dependencies = wheel_metadata(verified_download(wheel), version)
    verified_download(source)
    path = args.recipe / "meta.yaml"
    before = path.read_text()
    test_path = args.recipe / TEST_NAME
    test = Path(__file__).with_name("package_test.py").read_text()
    test_changed = not test_path.exists() or test_path.read_text() != test
    after = update_recipe(before, version, source, entries, python, dependencies, test_changed)
    path.write_text(after)
    test_path.write_text(test)
    result = {"version": version, "changed": before != after or test_changed,
              "entry_points": len(entries), "source_sha256": source["digests"]["sha256"]}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
