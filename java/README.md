# clothoid-halley — Java 21 reference implementation

Halley/Newton solver for the clothoid chord-length residual. Bit-for-bit
match (within `1e-9` m) with the Python reference on the 9,058-record
ProRail corpus.

## Build and test

```bash
mvn package                      # runs JUnit 5 golden-vector test, then builds shaded JAR
```

## Run the benchmark

```bash
java -jar target/clothoid-halley-1.0.0.jar
```

Prints one JSON object with `halley_us`, `newton_us`, iteration means,
and runtime info. Consumed by `python/run_all_benches.py`.

## API

```java
import com.grootstebozewolf.clothoid.ClothoidSolver;
import com.grootstebozewolf.clothoid.ClothoidSolver.Result;

Result r = ClothoidSolver.solveHalleyL(
    new double[] { 0.0, 0.0 },
    new double[] { 100.0, 0.0 },
    0.0,
    0.01);
// r.L() = arc length in metres; r.iterations() = steps taken
```
