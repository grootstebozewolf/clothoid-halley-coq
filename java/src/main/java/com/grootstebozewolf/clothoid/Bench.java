package com.grootstebozewolf.clothoid;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Benchmark harness: warmup + 50 timed repetitions of the full golden
 * vector corpus. Emits the same JSON shape as the C# and TypeScript
 * harnesses so the aggregator can merge them into one Table 1.
 */
public final class Bench {

    private static final int WARMUPS = 5;
    private static final int REPEATS = 50;

    public static void main(String[] args) throws IOException {
        Case[] cases = loadGoldenCases();

        // Warmup
        for (int w = 0; w < WARMUPS; w++) {
            runHalley(cases);
            runNewton(cases);
        }

        double[] halleyMs = new double[REPEATS];
        double[] newtonMs = new double[REPEATS];
        long halleyIters = 0, newtonIters = 0;
        for (int r = 0; r < REPEATS; r++) {
            long t0 = System.nanoTime();
            halleyIters = runHalley(cases);
            long t1 = System.nanoTime();
            halleyMs[r] = (t1 - t0) / 1_000_000.0;

            t0 = System.nanoTime();
            newtonIters = runNewton(cases);
            t1 = System.nanoTime();
            newtonMs[r] = (t1 - t0) / 1_000_000.0;
        }
        Arrays.sort(halleyMs);
        Arrays.sort(newtonMs);
        double halleyMedianMs = halleyMs[REPEATS / 2];
        double newtonMedianMs = newtonMs[REPEATS / 2];
        double halleyMinMs    = halleyMs[0];
        double newtonMinMs    = newtonMs[0];

        Map<String, Object> result = new LinkedHashMap<>();
        result.put("language",         "Java 21");
        result.put("runtime",          System.getProperty("java.version"));
        result.put("cases",            cases.length);
        result.put("halley_us",        1000.0 * halleyMedianMs / cases.length);
        result.put("newton_us",        1000.0 * newtonMedianMs / cases.length);
        result.put("halley_us_min",    1000.0 * halleyMinMs    / cases.length);
        result.put("newton_us_min",    1000.0 * newtonMinMs    / cases.length);
        result.put("halley_iter_mean", (double) halleyIters / cases.length);
        result.put("newton_iter_mean", (double) newtonIters / cases.length);
        System.out.println(new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(result));
    }

    private static long runHalley(Case[] cs) {
        long iters = 0;
        for (Case c : cs) {
            ClothoidSolver.Result r = ClothoidSolver.solveHalleyL(c.p0, c.p1, c.k0, c.k1);
            iters += r.iterations();
        }
        return iters;
    }

    private static long runNewton(Case[] cs) {
        long iters = 0;
        for (Case c : cs) {
            ClothoidSolver.Result r = ClothoidSolver.solveNewtonL(c.p0, c.p1, c.k0, c.k1);
            iters += r.iterations();
        }
        return iters;
    }

    private static Case[] loadGoldenCases() throws IOException {
        File path = findRepoRelative("data/golden_vectors.json");
        JsonNode root = new ObjectMapper().readTree(path);
        JsonNode cases = root.get("cases");
        Case[] out = new Case[cases.size()];
        for (int i = 0; i < cases.size(); i++) {
            JsonNode c = cases.get(i);
            out[i] = new Case(
                pair(c.get("P0")),
                pair(c.get("P1")),
                c.get("k0").asDouble(),
                c.get("k1").asDouble()
            );
        }
        return out;
    }

    private static double[] pair(JsonNode arr) {
        return new double[] { arr.get(0).asDouble(), arr.get(1).asDouble() };
    }

    static File findRepoRelative(String relative) {
        File dir = new File("").getAbsoluteFile();
        while (dir != null) {
            File candidate = new File(dir, relative.replace('/', File.separatorChar));
            if (candidate.isFile()) return candidate;
            dir = dir.getParentFile();
        }
        throw new IllegalStateException("could not locate '" + relative + "' walking up from cwd");
    }

    private record Case(double[] p0, double[] p1, double k0, double k1) {}
}
