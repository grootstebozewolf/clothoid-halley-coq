// SPDX-FileCopyrightText: 2026 Merkator Group
// SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { solveHalleyL, solveNewtonL } from "../src/solver.js";
import { findRepoRelative } from "../src/bench.js";

interface GoldenCase {
    objectid:    number;
    P0:          [number, number];
    P1:          [number, number];
    k0:          number;
    k1:          number;
    L:           number;
    iter_halley: number;
    iter_newton: number;
}

interface GoldenFile { cases: GoldenCase[] }

const AGREE_TOL = 1e-9;

test("Halley and Newton agree with Python reference on all ProRail cases", async () => {
    const path = findRepoRelative("data/golden_vectors.json");
    const doc  = JSON.parse(await readFile(path, "utf8")) as GoldenFile;
    assert.ok(doc.cases.length > 1000, `expected >1000 cases, got ${doc.cases.length}`);

    let halleyIterMatch = 0;
    let newtonIterMatch = 0;
    for (const c of doc.cases) {
        const rh = solveHalleyL(c.P0, c.P1, c.k0, c.k1);
        const rn = solveNewtonL(c.P0, c.P1, c.k0, c.k1);

        assert.ok(Math.abs(rh.L - c.L) < AGREE_TOL,
            `Halley L mismatch for OID ${c.objectid}: ts=${rh.L}, py=${c.L}`);
        assert.ok(Math.abs(rn.L - c.L) < AGREE_TOL,
            `Newton L mismatch for OID ${c.objectid}: ts=${rn.L}, py=${c.L}`);

        if (rh.iterations === c.iter_halley) halleyIterMatch++;
        if (rn.iterations === c.iter_newton) newtonIterMatch++;
    }
    const halleyAgreement = halleyIterMatch / doc.cases.length;
    const newtonAgreement = newtonIterMatch / doc.cases.length;
    assert.ok(halleyAgreement >= 0.99, `Halley iter agreement ${halleyAgreement} < 99%`);
    assert.ok(newtonAgreement >= 0.99, `Newton iter agreement ${newtonAgreement} < 99%`);
});
