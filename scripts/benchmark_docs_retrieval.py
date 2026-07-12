#!/usr/bin/env python3
"""Evaluate gold-answer support and top-five documentation retrieval."""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
BENCHMARK = DOCS / "_data" / "retrieval_benchmark.json"
TOKEN = re.compile(r"[a-z0-9_]+")
STOPWORDS = {
    "a", "an", "and", "are", "as", "at", "be", "by", "can", "do", "does",
    "for", "from", "how", "i", "in", "is", "it", "my", "of", "on", "or",
    "should", "the", "this", "to", "use", "what", "where", "which", "with",
}


def tokens(text: str) -> list[str]:
    return [token for token in TOKEN.findall(text.lower()) if token not in STOPWORDS]


def documents() -> dict[str, str]:
    paths = [
        DOCS / "getting_started" / "index.rst",
        DOCS / "formats" / "index.rst",
        DOCS / "glossary" / "index.rst",
        DOCS / "troubleshooting" / "index.rst",
        DOCS / "frequently_asked_questions" / "index.rst",
        *sorted((DOCS / "reference" / "commands").glob("*.rst")),
        *sorted((DOCS / "reference" / "categories").glob("*.rst")),
        *sorted((DOCS / "tutorials" / "pages").glob("*.rst")),
    ]
    result = {}
    for path in paths:
        key = str(path.relative_to(DOCS))
        content = path.read_text()
        include = re.search(r"\.\. include:: /_generated/commands/([^\s]+)", content)
        if include:
            content += "\n" + (DOCS / "_generated" / "commands" / include.group(1)).read_text()
        result[key] = content
    return result


def rank(query: str, corpus: dict[str, str]) -> list[str]:
    query_terms = tokens(query)
    tokenized = {key: tokens(text) for key, text in corpus.items()}
    document_frequency = Counter(
        term for words in tokenized.values() for term in set(words)
    )
    count = len(tokenized)
    average_length = sum(map(len, tokenized.values())) / count
    scores = {}
    for key, words in tokenized.items():
        frequencies = Counter(words)
        path_terms = set(tokens(key))
        lines = corpus[key].splitlines()
        headings = [
            line for index, line in enumerate(lines[:-1])
            if line.strip()
            and lines[index + 1].strip()
            and len(set(lines[index + 1].strip())) == 1
            and lines[index + 1].strip()[0] in "=-*#^"
        ]
        title_terms = set(tokens("\n".join(headings)))
        score = 0.0
        for term in query_terms:
            if frequencies[term]:
                inverse_frequency = math.log(
                    1 + (count - document_frequency[term] + 0.5)
                    / (document_frequency[term] + 0.5)
                )
                frequency = frequencies[term]
                normalization = frequency + 1.5 * (
                    1 - 0.75 + 0.75 * len(words) / average_length
                )
                score += inverse_frequency * frequency * 2.5 / normalization
                if term in title_terms:
                    score += inverse_frequency * 1.5
                if term in path_terms:
                    score += inverse_frequency * 3
        scores[key] = score
    return sorted(scores, key=lambda key: (-scores[key], key))


def main() -> int:
    benchmark = json.loads(BENCHMARK.read_text())
    questions = benchmark["questions"]
    if len(questions) < 50:
        raise SystemExit("retrieval benchmark must contain at least 50 questions")
    corpus = documents()
    retrieval_passes = 0
    support_passes = 0
    failures = []
    for item in questions:
        missing_sources = [source for source in item["sources"] if source not in corpus]
        if missing_sources:
            failures.append({"id": item["id"], "missing_sources": missing_sources})
            continue
        top_five = rank(item["question"], corpus)[:5]
        retrieved = any(source in top_five for source in item["sources"])
        source_text = " ".join(
            " ".join(corpus[source].lower().split()) for source in item["sources"]
        )
        supported = all(
            " ".join(term.lower().split()) in source_text
            for term in item["answer_terms"]
        )
        retrieval_passes += int(retrieved)
        support_passes += int(supported)
        if not retrieved or not supported:
            failures.append(
                {"id": item["id"], "retrieved": retrieved, "supported": supported, "top_five": top_five}
            )
    retrieval_accuracy = retrieval_passes / len(questions)
    factual_accuracy = support_passes / len(questions)
    print(f"questions: {len(questions)}")
    print(f"gold-answer factual support: {factual_accuracy:.1%}")
    print(f"canonical source in top five: {retrieval_accuracy:.1%}")
    if failures:
        print(json.dumps(failures, indent=2))
    return 0 if factual_accuracy >= 0.95 and retrieval_accuracy >= 0.95 else 1


if __name__ == "__main__":
    raise SystemExit(main())
