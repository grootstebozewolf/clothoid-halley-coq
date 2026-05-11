# clothoid-halley — TypeScript / Node.js reference implementation

Halley/Newton solver for the clothoid chord-length residual. Bit-for-bit
match (within `1e-9` m) with the Python reference on the 9,058-record
ProRail corpus. Zero runtime dependencies; only `typescript` and
`@types/node` are needed to build.

## Build and test

```bash
npm install
npm test                         # tsc + node --test on golden vectors
```

## Run the benchmark

```bash
npm run bench
```

Prints one JSON object with `halley_us`, `newton_us`, iteration means,
and runtime info. Consumed by `python/run_all_benches.py`.

## API

```ts
import { solveHalleyL } from "clothoid-halley";

const r = solveHalleyL([0, 0], [100, 0], 0, 0.01);
// r.L = arc length in metres; r.iterations = steps taken
```
