// SPDX-FileCopyrightText: 2026 Merkator Group
// SPDX-License-Identifier: EUPL-1.2

using System.IO;
using System.Text.Json;
using System.Text.Json.Serialization;
using Xunit;

namespace Clothoid.Halley.Tests;

public sealed class GoldenVectorTests
{
    /// <summary>1 nanometre — comfortably below any realistic measurement precision.</summary>
    private const double AgreeTol = 1e-9;

    [Fact]
    public void HalleyAndNewton_AgreeWithPythonReference_OnAllProRailCases()
    {
        var golden = LoadGoldenVectors();
        Assert.True(golden.Count > 1000, $"expected >1000 cases, got {golden.Count}");

        int halleyIterMatch = 0;
        int newtonIterMatch = 0;
        foreach (var c in golden)
        {
            var rh = ClothoidSolver.SolveHalleyL(c.P0, c.P1, c.K0, c.K1);
            var rn = ClothoidSolver.SolveNewtonL(c.P0, c.P1, c.K0, c.K1);

            Assert.True(System.Math.Abs(rh.L - c.L) < AgreeTol,
                $"Halley L mismatch for OID {c.Objectid}: csharp={rh.L:R}, python={c.L:R}, dL={rh.L - c.L:R}");
            Assert.True(System.Math.Abs(rn.L - c.L) < AgreeTol,
                $"Newton L mismatch for OID {c.Objectid}: csharp={rn.L:R}, python={c.L:R}, dL={rn.L - c.L:R}");

            if (rh.Iterations == c.IterHalley) halleyIterMatch++;
            if (rn.Iterations == c.IterNewton) newtonIterMatch++;
        }

        // Iteration counts should match the Python reference on >99 % of cases.
        // A small allowance covers harmless platform-dependent rounding in the
        // last bit that flips the convergence test one iteration early/late.
        double halleyAgreement = (double)halleyIterMatch / golden.Count;
        double newtonAgreement = (double)newtonIterMatch / golden.Count;
        Assert.True(halleyAgreement >= 0.99,
            $"Halley iteration agreement {halleyAgreement:P2} < 99% ({halleyIterMatch}/{golden.Count})");
        Assert.True(newtonAgreement >= 0.99,
            $"Newton iteration agreement {newtonAgreement:P2} < 99% ({newtonIterMatch}/{golden.Count})");
    }

    private static List<GoldenCase> LoadGoldenVectors()
    {
        var path = FindRepoRelative("data/golden_vectors.json");
        using var fs = File.OpenRead(path);
        var doc = JsonSerializer.Deserialize<GoldenFile>(fs, JsonOpts)
            ?? throw new InvalidDataException("failed to parse golden_vectors.json");
        return doc.Cases;
    }

    private static readonly JsonSerializerOptions JsonOpts = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        ReadCommentHandling  = JsonCommentHandling.Skip,
        AllowTrailingCommas  = true,
    };

    private static string FindRepoRelative(string relative)
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            var candidate = Path.Combine(dir.FullName, relative.Replace('/', Path.DirectorySeparatorChar));
            if (File.Exists(candidate)) return candidate;
            dir = dir.Parent;
        }
        throw new FileNotFoundException(
            $"could not locate '{relative}' walking up from {AppContext.BaseDirectory}");
    }

    private sealed class GoldenFile
    {
        [JsonPropertyName("cases")]    public List<GoldenCase> Cases { get; set; } = new();
        [JsonPropertyName("count")]    public int Count { get; set; }
    }

    public sealed class GoldenCase
    {
        [JsonPropertyName("objectid")]    public long      Objectid   { get; set; }
        [JsonPropertyName("P0")]          public double[]  P0         { get; set; } = System.Array.Empty<double>();
        [JsonPropertyName("P1")]          public double[]  P1         { get; set; } = System.Array.Empty<double>();
        [JsonPropertyName("k0")]          public double    K0         { get; set; }
        [JsonPropertyName("k1")]          public double    K1         { get; set; }
        [JsonPropertyName("L_design")]    public double    LDesign    { get; set; }
        [JsonPropertyName("d")]           public double    D          { get; set; }
        [JsonPropertyName("L")]           public double    L          { get; set; }
        [JsonPropertyName("iter_halley")] public int       IterHalley { get; set; }
        [JsonPropertyName("iter_newton")] public int       IterNewton { get; set; }
    }
}
