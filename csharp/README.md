# Clothoid.Halley — C# (.NET 8) reference implementation

Halley/Newton solver for the clothoid chord-length residual. Bit-for-bit
match (within `1e-9` m) with the Python reference on the 9,058-record
ProRail corpus.

`Clothoid.Halley/Solver.cs` and `Clothoid.Halley/GaussLegendre.cs` are
BSD-3-Clause, copyright (c) 2026 Merkator Group. Copy both files together
into NetTopologySuite; they have no other dependencies in this repo.
Every other file here remains EUPL-1.2. The library's NuGet licence
expression is `BSD-3-Clause`.

## Build and test

```bash
dotnet test -c Release           # runs xUnit golden-vector test
```

## Run the benchmark

```bash
dotnet run -c Release --project Clothoid.Halley.Bench
```

Prints one JSON line with `halley_us`, `newton_us`, iteration means,
and runtime info. Consumed by `python/run_all_benches.py`.

## API

```csharp
using Clothoid.Halley;

SolverResult r = ClothoidSolver.SolveHalleyL(
    p0: new[] { 0.0, 0.0 },
    p1: new[] { 100.0, 0.0 },
    k0: 0.0,
    k1: 0.01);
// r.L = arc length in metres; r.Iterations = steps taken
```
