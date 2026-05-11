# Clothoid.Halley — C# (.NET 8) reference implementation

Halley/Newton solver for the clothoid chord-length residual. Bit-for-bit
match (within `1e-9` m) with the Python reference on the 9,058-record
ProRail corpus.

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
