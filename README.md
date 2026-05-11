# clothoid-halley-coq

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20128416.svg)](https://doi.org/10.5281/zenodo.20128416)

A robust Halley solver for the chord-length parameter of the clothoid
(Euler-spiral) $G^1$ Hermite interpolation problem, with a Coq
formalisation of the underlying derivative identities, three production
implementations (C# / Java / TypeScript), and a 9,058-record real-world
benchmark on the ProRail Spoorgeometrie dataset.

> **Repository status: public — proprietary (source-available).**
> The full source is open for reading and academic evaluation. It is
> **not open source**: see the [License](#license) section below before
> running, modifying, or redistributing any of it.
> The accompanying paper is at
> [`docs/mathematics/Clothoid_L_Halley_Solver.pdf`](docs/mathematics/Clothoid_L_Halley_Solver.pdf).

---

## Contents

### `coq/` — Formal verification

Three Coq files; all compile cleanly under Coq 8.13.1 and Coq 8.20.1 with
Coquelicot 3.x, pass `coqchk` with no `type-in-type`, no unsafe
(co)fixpoints, no positivity holes, and **no `Admitted` or `Axiom`**
beyond the four standard axioms used by Coquelicot itself (classical
logic, decidable Dedekind reals, functional extensionality).

- **`Clothoid.v`** — Bertolazzi–Frego $G^1$ residual $f(A) = Y_0(2A, \delta - A, \varphi_0)$. Proves $f'(A) = \int(t^2-t)\cos\varphi$ and $f''(A) = -\int(t^2-t)^2\sin\varphi$ via parameter-differentiation under the integral.
- **`Clothoid_L.v`** — chord-length residual $f(L) = L^2(P^2 + Q^2) - d^2$. Proves the four integral-level derivatives $P' = -T$, $Q' = R$, $R' = -S_{2s}$, $T' = S_{2c}$, the first-derivative composite, and the closed (no `Admitted`) second-derivative composite via `auto_derive` + `Derive` rewrites + `ring`.
- **`ClothoidPolish.v`** — moment-notation polish: proves $f'(A) = X_2 - X_1$ and the integral identity $-(Y_4 - 2 Y_3 + Y_2) = \int -(t^2-t)^2 \sin\varphi$.

```bash
cd coq && make
```

### `python/` — Reference implementation, derivations, benchmarks

- **Symbolic derivations** (`clothoid_*_derive*.py`): SymPy-verified $f'$, $f''$ for the $A$- and $L$-formulations.
- **Reference solver** (`clothoid_halleyL_bench.py`): the canonical Halley / Newton solver every other implementation is bit-compared against.
- **Data pipeline**:
  - `fetch_prorail_clothoids.py` → fetches all `Overgangsboog` records from the ProRail Spoorgeometrie ArcGIS service into `data/prorail_clothoids.json.gz` (CC BY 4.0).
  - `build_golden_vectors.py` → filters monotone-branch + runs both solvers → `data/golden_vectors.json` (9,058 cases).
  - `bench_corpus.py` + `run_all_benches.py` → drives the four-language benchmark, writes `data/benchmark_results.json`.

### `csharp/` — C# (.NET 8) implementation

`Clothoid.Halley` library + xUnit golden-vector tests + benchmark harness. See [csharp/README.md](csharp/README.md). Bit-identical to the Python reference on every case in the 9,058-record corpus (chord-length agreement within $10^{-9}$ m, iteration counts match exactly). Median Halley solve: **0.59 µs**.

### `java/` — Java 21 implementation

Maven project producing a shaded JAR + JUnit 5 golden-vector tests + benchmark. See [java/README.md](java/README.md). Same numerical agreement guarantees. Median Halley solve: **1.38 µs**.

### `typescript/` — TypeScript / Node.js implementation

Zero-runtime-dependency ESM module + `node:test` golden-vector suite + benchmark. See [typescript/README.md](typescript/README.md). Same numerical agreement guarantees. Median Halley solve: **0.88 µs**.

### `docs/mathematics/` — Academic paper

- `Clothoid_L_Halley_Solver.tex` — the LaTeX source.
- `Clothoid_L_Halley_Solver.pdf` — the rendered paper (7 pages, includes the formal-verification section and the cross-language benchmark table).
- `references.bib` — the bibliography (Bertolazzi & Frego 2015 / 2018, Coquelicot, Householder, Vázquez-Méndez & Casal).
- `generate_benchmark_graphs.py` — renders the bar charts from `data/benchmark_results.json`.

### `data/` — Datasets

- `prorail_clothoids.json.gz` — raw snapshot of 9,058 ProRail clothoid transitions (CC BY 4.0 ProRail Spoorgeometrie; see `data/LICENSE_DATA.txt`).
- `golden_vectors.json` — filtered + solved corpus, consumed by every language's test and benchmark harness.
- `benchmark_results.json` — measured per-language Halley / Newton numbers.

---

## End-to-end reproduction

```bash
# 1. Re-fetch the ProRail snapshot (skip if you trust the committed copy).
python python/fetch_prorail_clothoids.py

# 2. Rebuild golden vectors from the snapshot.
python python/build_golden_vectors.py

# 3. Run the four-language tests.
dotnet test  csharp/Clothoid.Halley.Tests
mvn -f java/pom.xml test
cd typescript && npm install && npm test && cd ..

# 4. Run the cross-language benchmark.
python python/run_all_benches.py

# 5. Regenerate the bar charts.
python docs/mathematics/generate_benchmark_graphs.py

# 6. Rebuild the PDF (requires pdflatex / texlive).
cd docs/mathematics && pdflatex Clothoid_L_Halley_Solver.tex \
                    && bibtex   Clothoid_L_Halley_Solver \
                    && pdflatex Clothoid_L_Halley_Solver.tex \
                    && pdflatex Clothoid_L_Halley_Solver.tex
```

Toolchain used for the committed numbers: Coq 8.20.1 + Coquelicot 3.x;
Python 3.14 + NumPy + SciPy + Matplotlib; .NET 8 / .NET 10 SDK; OpenJDK 21
(Corretto); Node.js 22; TeX Live 2026.

---

## License

This repository is **public** but **not open source**.

- All code, documentation, build configuration, and the manuscript text
  are covered by a proprietary licence — see [LICENSE](LICENSE) — that
  grants narrow permissions for *reading*, *academic citation*, and
  *unmodified execution to reproduce the paper's results*. All other
  uses (commercial use, redistribution, derivative works) require a
  separate written licence from Merkator Group.
- The two data files `data/prorail_clothoids.json.gz` and
  `data/golden_vectors.json` are derived from ProRail Spoorgeometrie
  and are redistributed under **CC BY 4.0** (see
  [data/LICENSE_DATA.txt](data/LICENSE_DATA.txt)). That licence applies
  to those files only; it is not extended to the surrounding code.

### No AI / ML training

Use of this repository — in whole or in part, in any form, by any
party, including automated crawlers — to train, fine-tune, retrieval-
augment, distil, evaluate, benchmark, or otherwise develop any
artificial-intelligence or machine-learning system, or to assemble any
dataset for those purposes, is **expressly prohibited**. This is an
explicit text-and-data-mining reservation within the meaning of
Article 4(3) of the EU DSM Directive (Directive (EU) 2019/790) and the
corresponding opt-outs in other jurisdictions. See the **"NO USE FOR
TRAINING OF AI / MACHINE-LEARNING SYSTEMS"** section of
[LICENSE](LICENSE) for the full terms, including the explicit
preservation of access logs as evidence in subsequent proceedings.

Licensing enquiries: <jeroen.bloemscheer@merkator.com>.
