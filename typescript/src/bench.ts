import { readFile } from "node:fs/promises";
import { existsSync } from "node:fs";
import * as path from "node:path";
import { fileURLToPath } from "node:url";
import { solveHalleyL, solveNewtonL } from "./solver.js";

const WARMUPS = 5;
const REPEATS = 50;

interface GoldenCase {
    P0:           [number, number];
    P1:           [number, number];
    k0:           number;
    k1:           number;
}

interface GoldenFile {
    cases: GoldenCase[];
}

async function main(): Promise<void> {
    const cases = await loadGoldenCases();

    // Warmup
    for (let w = 0; w < WARMUPS; w++) {
        runHalley(cases);
        runNewton(cases);
    }

    const halleyMs = new Float64Array(REPEATS);
    const newtonMs = new Float64Array(REPEATS);
    let halleyIters = 0, newtonIters = 0;
    for (let r = 0; r < REPEATS; r++) {
        let t0 = process.hrtime.bigint();
        halleyIters = runHalley(cases);
        let t1 = process.hrtime.bigint();
        halleyMs[r] = Number(t1 - t0) / 1_000_000;

        t0 = process.hrtime.bigint();
        newtonIters = runNewton(cases);
        t1 = process.hrtime.bigint();
        newtonMs[r] = Number(t1 - t0) / 1_000_000;
    }
    const sorted = (a: Float64Array) => Array.from(a).sort((x, y) => x - y);
    const halleyArr = sorted(halleyMs);
    const newtonArr = sorted(newtonMs);
    const halleyMedianMs = halleyArr[Math.floor(REPEATS / 2)]!;
    const newtonMedianMs = newtonArr[Math.floor(REPEATS / 2)]!;
    const halleyMinMs    = halleyArr[0]!;
    const newtonMinMs    = newtonArr[0]!;

    const result = {
        language:         "TypeScript (Node.js)",
        runtime:          process.version,
        cases:            cases.length,
        halley_us:        (1000 * halleyMedianMs) / cases.length,
        newton_us:        (1000 * newtonMedianMs) / cases.length,
        halley_us_min:    (1000 * halleyMinMs)    / cases.length,
        newton_us_min:    (1000 * newtonMinMs)    / cases.length,
        halley_iter_mean: halleyIters / cases.length,
        newton_iter_mean: newtonIters / cases.length,
    };
    console.log(JSON.stringify(result, null, 2));
}

function runHalley(cs: GoldenCase[]): number {
    let iters = 0;
    for (const c of cs) {
        iters += solveHalleyL(c.P0, c.P1, c.k0, c.k1).iterations;
    }
    return iters;
}

function runNewton(cs: GoldenCase[]): number {
    let iters = 0;
    for (const c of cs) {
        iters += solveNewtonL(c.P0, c.P1, c.k0, c.k1).iterations;
    }
    return iters;
}

async function loadGoldenCases(): Promise<GoldenCase[]> {
    const p = findRepoRelative("data/golden_vectors.json");
    const buf = await readFile(p, "utf8");
    const doc = JSON.parse(buf) as GoldenFile;
    return doc.cases;
}

export function findRepoRelative(rel: string): string {
    let dir = path.dirname(fileURLToPath(import.meta.url));
    while (dir && dir !== path.dirname(dir)) {
        const candidate = path.join(dir, rel.split("/").join(path.sep));
        if (existsSync(candidate)) return candidate;
        dir = path.dirname(dir);
    }
    throw new Error(`could not locate '${rel}' walking up from import.meta.url`);
}

void main();
