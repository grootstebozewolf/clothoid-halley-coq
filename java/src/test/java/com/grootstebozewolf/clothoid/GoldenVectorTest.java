// SPDX-FileCopyrightText: 2026 Merkator Group
// SPDX-License-Identifier: EUPL-1.2

package com.grootstebozewolf.clothoid;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.jupiter.api.Assertions.assertTrue;

public final class GoldenVectorTest {

    /** 1 nanometre — well below any realistic measurement precision. */
    private static final double AGREE_TOL = 1e-9;

    @Test
    void halleyAndNewton_agreeWithPythonReference_onAllProRailCases() throws IOException {
        File path = Bench.findRepoRelative("data/golden_vectors.json");
        JsonNode root = new ObjectMapper().readTree(path);
        JsonNode cases = root.get("cases");
        assertTrue(cases.size() > 1000, "expected >1000 cases, got " + cases.size());

        int halleyIterMatch = 0;
        int newtonIterMatch = 0;
        for (int i = 0; i < cases.size(); i++) {
            JsonNode c = cases.get(i);
            double[] P0 = { c.get("P0").get(0).asDouble(), c.get("P0").get(1).asDouble() };
            double[] P1 = { c.get("P1").get(0).asDouble(), c.get("P1").get(1).asDouble() };
            double k0 = c.get("k0").asDouble();
            double k1 = c.get("k1").asDouble();
            double Lpy = c.get("L").asDouble();
            int iterH = c.get("iter_halley").asInt();
            int iterN = c.get("iter_newton").asInt();
            long oid = c.get("objectid").asLong();

            ClothoidSolver.Result rh = ClothoidSolver.solveHalleyL(P0, P1, k0, k1);
            ClothoidSolver.Result rn = ClothoidSolver.solveNewtonL(P0, P1, k0, k1);

            assertTrue(Math.abs(rh.L() - Lpy) < AGREE_TOL,
                    "Halley L mismatch for OID " + oid + ": java=" + rh.L() + " python=" + Lpy);
            assertTrue(Math.abs(rn.L() - Lpy) < AGREE_TOL,
                    "Newton L mismatch for OID " + oid + ": java=" + rn.L() + " python=" + Lpy);

            if (rh.iterations() == iterH) halleyIterMatch++;
            if (rn.iterations() == iterN) newtonIterMatch++;
        }

        double halleyAgreement = (double) halleyIterMatch / cases.size();
        double newtonAgreement = (double) newtonIterMatch / cases.size();
        assertTrue(halleyAgreement >= 0.99,
                "Halley iteration agreement " + halleyAgreement + " < 99%");
        assertTrue(newtonAgreement >= 0.99,
                "Newton iteration agreement " + newtonAgreement + " < 99%");
    }
}
