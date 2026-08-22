(* SPDX-FileCopyrightText: 2026 NetTopologySuite.Proofs contributors *)
(* SPDX-License-Identifier: BSD-3-Clause *)

(* ============================================================================
   NetTopologySuite.Proofs.Distance
   Vendored into clothoid-halley-coq/coq/ under BSD-3-Clause only
   (see coq/LICENSE_BSD-3-Clause.txt). Not licensed as EUPL-1.2.
   ----------------------------------------------------------------------------
   Foundational properties of Euclidean distance on the 2D plane.

   The optimisation tracked in JTS PR #1111 (and the gap noted in NTS epic #828
   for LineStringSnapper) replaces `distance(p, q) < tol` with
   `distance_sq(p, q) < tol * tol`, saving the square root call.

   The justification rests on a small theorem: on the non-negative reals,
   x ≤ y  iff  x² ≤ y².  Since Euclidean distance is non-negative, the squared
   form of any distance comparison against a non-negative threshold is
   equivalent.  This file proves it cleanly.

   Author: NetTopologySuite.Proofs contributors
   License: BSD-3-Clause (see coq/LICENSE_BSD-3-Clause.txt)
   ========================================================================== *)

Require Import Reals.
Require Import Lra.
Open Scope R_scope.

(* -------------------------------------------------------------------------- *)
(* A 2-D point is a pair of reals.                                            *)
(* -------------------------------------------------------------------------- *)

Record Point : Type := mkPoint { px : R; py : R }.

Definition dist_sq (p q : Point) : R :=
  (px p - px q) * (px p - px q) + (py p - py q) * (py p - py q).

Definition dist (p q : Point) : R := sqrt (dist_sq p q).

(* -------------------------------------------------------------------------- *)
(* Foundational properties of squared distance.                               *)
(* -------------------------------------------------------------------------- *)

Lemma sqr_nonneg : forall x : R, 0 <= x * x.
Proof.
  intros x. pose proof (Rle_0_sqr x) as H. unfold Rsqr in H. exact H.
Qed.

Lemma sqr_eq_zero : forall x : R, x * x = 0 -> x = 0.
Proof.
  intros x H. apply Rsqr_0_uniq. unfold Rsqr. exact H.
Qed.

Lemma dist_sq_nonneg : forall p q, 0 <= dist_sq p q.
Proof.
  intros p q. unfold dist_sq.
  pose proof (sqr_nonneg (px p - px q)) as Hx.
  pose proof (sqr_nonneg (py p - py q)) as Hy.
  lra.
Qed.

Lemma dist_sq_sym : forall p q, dist_sq p q = dist_sq q p.
Proof.
  intros p q. unfold dist_sq.
  ring.
Qed.

Lemma dist_sq_zero_iff_eq : forall p q,
  dist_sq p q = 0 <-> (px p = px q /\ py p = py q).
Proof.
  intros p q. unfold dist_sq. split.
  - intros H.
    pose proof (sqr_nonneg (px p - px q)) as Hx.
    pose proof (sqr_nonneg (py p - py q)) as Hy.
    assert (Hxz : (px p - px q) * (px p - px q) = 0) by lra.
    assert (Hyz : (py p - py q) * (py p - py q) = 0) by lra.
    apply sqr_eq_zero in Hxz.
    apply sqr_eq_zero in Hyz.
    split; lra.
  - intros [Hx Hy]. rewrite Hx, Hy. ring.
Qed.

(* -------------------------------------------------------------------------- *)
(* The headline theorem: on non-negative reals, the square is monotone.       *)
(* This is the formal justification for replacing distance comparisons with   *)
(* squared-distance comparisons.                                              *)
(* -------------------------------------------------------------------------- *)

Lemma sq_monotone_nonneg : forall x y, 0 <= x -> 0 <= y ->
  (x <= y <-> x * x <= y * y).
Proof.
  intros x y Hx Hy.
  split; intros H.
  - apply Rmult_le_compat; lra.
  - destruct (Rle_or_lt x y) as [Hle | Hlt].
    + exact Hle.
    + exfalso.
      assert (y * y < x * x).
      { apply Rmult_le_0_lt_compat; lra. }
      lra.
Qed.

(* -------------------------------------------------------------------------- *)
(* The result NTS would cite: comparing a distance to a non-negative          *)
(* threshold is equivalent to comparing squared distance to the squared       *)
(* threshold.                                                                 *)
(*                                                                            *)
(*   For all points p q and threshold t with t >= 0,                          *)
(*     dist(p, q) <= t   iff   dist_sq(p, q) <= t * t.                        *)
(*                                                                            *)
(* This is the formal "no, the optimisation cannot lie" guarantee.            *)
(* -------------------------------------------------------------------------- *)

Theorem dist_le_iff_dist_sq_le : forall p q t,
  0 <= t -> (dist p q <= t <-> dist_sq p q <= t * t).
Proof.
  intros p q t Ht.
  unfold dist.
  pose proof (dist_sq_nonneg p q) as Hd.
  pose proof (sqrt_pos (dist_sq p q)) as Hsd.
  rewrite (sq_monotone_nonneg (sqrt (dist_sq p q)) t Hsd Ht).
  rewrite sqrt_sqrt by exact Hd.
  reflexivity.
Qed.

(* -------------------------------------------------------------------------- *)
(* Distance (with sqrt) properties.                                           *)
(* -------------------------------------------------------------------------- *)

Lemma dist_nonneg : forall p q, 0 <= dist p q.
Proof.
  intros p q. unfold dist. apply sqrt_pos.
Qed.

Lemma dist_refl : forall p, dist p p = 0.
Proof.
  intros p. unfold dist.
  replace (dist_sq p p) with 0.
  - apply sqrt_0.
  - unfold dist_sq. ring.
Qed.

Lemma dist_sym : forall p q, dist p q = dist q p.
Proof.
  intros p q. unfold dist. rewrite (dist_sq_sym p q). reflexivity.
Qed.

Theorem dist_eq_zero_iff : forall p q,
  dist p q = 0 <-> (px p = px q /\ py p = py q).
Proof.
  intros p q. unfold dist. split.
  - intros H.
    apply sqrt_eq_0 in H; [| apply dist_sq_nonneg].
    apply dist_sq_zero_iff_eq. exact H.
  - intros Hxy.
    replace (dist_sq p q) with 0.
    + apply sqrt_0.
    + symmetry. apply dist_sq_zero_iff_eq. exact Hxy.
Qed.

(* -------------------------------------------------------------------------- *)
(* Squared-distance properties (translation invariance, scaling, strictness). *)
(* -------------------------------------------------------------------------- *)

Definition pt_translate (p : Point) (vx vy : R) : Point :=
  mkPoint (px p + vx) (py p + vy).

Lemma dist_sq_translation_invariant : forall p q vx vy,
  dist_sq (pt_translate p vx vy) (pt_translate q vx vy) = dist_sq p q.
Proof.
  intros p q vx vy. unfold dist_sq, pt_translate. simpl. ring.
Qed.

Definition pt_scale (c : R) (p : Point) : Point :=
  mkPoint (c * px p) (c * py p).

Lemma dist_sq_scale : forall c p q,
  dist_sq (pt_scale c p) (pt_scale c q) = c * c * dist_sq p q.
Proof.
  intros c p q. unfold dist_sq, pt_scale. simpl. ring.
Qed.

Lemma dist_sq_pos_iff_distinct : forall p q,
  0 < dist_sq p q <-> ~ (px p = px q /\ py p = py q).
Proof.
  intros p q. split.
  - intros H Hxy. apply dist_sq_zero_iff_eq in Hxy. lra.
  - intros H.
    pose proof (dist_sq_nonneg p q) as Hnn.
    destruct (Req_dec (dist_sq p q) 0) as [Heq | Hne].
    + exfalso. apply H. apply dist_sq_zero_iff_eq. exact Heq.
    + lra.
Qed.

(* -------------------------------------------------------------------------- *)
(* Coordinate decompositions and identities.                                  *)
(* -------------------------------------------------------------------------- *)

Lemma dist_sq_xy_decompose : forall p q,
  dist_sq p q = (px p - px q) * (px p - px q) + (py p - py q) * (py p - py q).
Proof. intros. unfold dist_sq. reflexivity. Qed.

Lemma dist_sq_x_only : forall p q,
  py p = py q -> dist_sq p q = (px p - px q) * (px p - px q).
Proof.
  intros p q H. unfold dist_sq. rewrite H. ring.
Qed.

Lemma dist_sq_y_only : forall p q,
  px p = px q -> dist_sq p q = (py p - py q) * (py p - py q).
Proof.
  intros p q H. unfold dist_sq. rewrite H. ring.
Qed.

Lemma dist_sq_neg_x : forall p q,
  dist_sq (mkPoint (- px p) (py p)) (mkPoint (- px q) (py q)) = dist_sq p q.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_neg_y : forall p q,
  dist_sq (mkPoint (px p) (- py p)) (mkPoint (px q) (- py q)) = dist_sq p q.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_self_zero : forall p, dist_sq p p = 0.
Proof. intros p. unfold dist_sq. ring. Qed.

Lemma dist_sq_at_origin_x : forall a, dist_sq (mkPoint 0 0) (mkPoint a 0) = a * a.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_at_origin_y : forall b, dist_sq (mkPoint 0 0) (mkPoint 0 b) = b * b.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_general : forall a b c d,
  dist_sq (mkPoint a b) (mkPoint c d) = (a - c) * (a - c) + (b - d) * (b - d).
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_pythagorean : forall a b,
  dist_sq (mkPoint 0 0) (mkPoint a b) = a * a + b * b.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_negate_both : forall p q,
  dist_sq (mkPoint (- px p) (- py p)) (mkPoint (- px q) (- py q)) = dist_sq p q.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_swap_xy : forall p q,
  dist_sq (mkPoint (py p) (px p)) (mkPoint (py q) (px q)) = dist_sq p q.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_triangular_sq_form : forall a b c d e f,
  dist_sq (mkPoint a b) (mkPoint e f) + dist_sq (mkPoint c d) (mkPoint e f)
  - 2 * ((a - e) * (c - e) + (b - f) * (d - f))
  = dist_sq (mkPoint a b) (mkPoint c d).
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_le_sum_xy : forall p q,
  dist_sq p q <= 2 * ((px p - px q) * (px p - px q) + (py p - py q) * (py p - py q)).
Proof.
  intros. unfold dist_sq.
  pose proof (sqr_nonneg (px p - px q)).
  pose proof (sqr_nonneg (py p - py q)). lra.
Qed.

Lemma dist_sq_zero_at_same_coord : forall a b,
  dist_sq (mkPoint a b) (mkPoint a b) = 0.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_diff_x_only_zero : forall a b,
  dist_sq (mkPoint a 0) (mkPoint a b) = b * b.
Proof. intros. unfold dist_sq. simpl. ring. Qed.

Lemma dist_sq_pos_when_x_diff : forall a c b,
  a <> c ->
  0 < dist_sq (mkPoint a b) (mkPoint c b).
Proof.
  intros a c b H. unfold dist_sq. simpl.
  assert (Hne : a - c <> 0) by lra.
  pose proof (Rle_0_sqr (a - c)). unfold Rsqr in *.
  destruct (Req_dec ((a - c) * (a - c)) 0) as [Heq | Hne0].
  - exfalso. apply H. assert (a - c = 0) by (apply Rsqr_0_uniq; unfold Rsqr; exact Heq).
    lra.
  - assert ((b - b) * (b - b) = 0) by (assert (b - b = 0) by ring; rewrite H1; ring).
    lra.
Qed.

Lemma dist_sq_eq_dist_eq : forall p q r s,
  px p = px r -> py p = py r -> px q = px s -> py q = py s ->
  dist_sq p q = dist_sq r s.
Proof.
  intros p q r s Hpx Hpy Hqx Hqy. unfold dist_sq.
  rewrite Hpx, Hpy, Hqx, Hqy. reflexivity.
Qed.

(* -------------------------------------------------------------------------- *)
(* Distance as a genuine metric: the sqrt bridge, the strict comparison form,  *)
(* positivity, and the triangle inequality.                                    *)
(*                                                                            *)
(* Together with dist_nonneg / dist_eq_zero_iff / dist_sym (above), the        *)
(* triangle inequality makes (Point, dist) a metric space — the foundational   *)
(* fact every distance-based predicate (nearest-point, snapping, buffering)    *)
(* silently relies on.                                                        *)
(* -------------------------------------------------------------------------- *)

(* The sqrt bridge: dist squared (the value, not dist_sq) recovers dist_sq.    *)
Lemma dist_mul_self : forall p q, dist p q * dist p q = dist_sq p q.
Proof.
  intros p q. unfold dist. apply sqrt_sqrt. apply dist_sq_nonneg.
Qed.

(* Strict counterpart of sq_monotone_nonneg. *)
Lemma sq_strict_monotone_nonneg : forall x y, 0 <= x -> 0 <= y ->
  (x < y <-> x * x < y * y).
Proof.
  intros x y Hx Hy. split; intros H.
  - apply Rmult_le_0_lt_compat; lra.
  - destruct (Rle_or_lt y x) as [Hle | Hlt].
    + exfalso. assert (y * y <= x * x) by (apply Rmult_le_compat; lra). lra.
    + exact Hlt.
Qed.

(* The strict form of the squared-distance fast path: matches the `distance <  *)
(* tol` comparison in the JTS PR #1111 / NTS #828 snapping optimisation        *)
(* (the file header's motivating case), which the existing `<=` row does not   *)
(* cover.                                                                      *)
Theorem dist_lt_iff_dist_sq_lt : forall p q t,
  0 <= t -> (dist p q < t <-> dist_sq p q < t * t).
Proof.
  intros p q t Ht. unfold dist.
  pose proof (dist_sq_nonneg p q) as Hd.
  pose proof (sqrt_pos (dist_sq p q)) as Hsd.
  rewrite (sq_strict_monotone_nonneg (sqrt (dist_sq p q)) t Hsd Ht).
  rewrite sqrt_sqrt by exact Hd. reflexivity.
Qed.

(* Positivity at the distance level (companion to dist_sq_pos_iff_distinct). *)
Lemma dist_pos_iff_distinct : forall p q,
  0 < dist p q <-> ~ (px p = px q /\ py p = py q).
Proof.
  intros p q. split.
  - intros H Hxy. apply dist_eq_zero_iff in Hxy. lra.
  - intros H. pose proof (dist_nonneg p q) as Hnn.
    destruct (Req_dec (dist p q) 0) as [Heq | Hne].
    + exfalso. apply H. apply dist_eq_zero_iff. exact Heq.
    + lra.
Qed.

(* 2-D Cauchy–Schwarz: |u·v|² ≤ |u|²|v|², from the Lagrange identity           *)
(* (a²+b²)(c²+d²) − (ac+bd)² = (ad−bc)² ≥ 0.                                   *)
Lemma cauchy_schwarz_2d : forall a b c d,
  (a * c + b * d) * (a * c + b * d) <= (a * a + b * b) * (c * c + d * d).
Proof.
  intros a b c d. pose proof (sqr_nonneg (a * d - b * c)). nra.
Qed.

(* The triangle inequality.  Proof: with u = p−q, v = q−r, the dot product     *)
(* dot = u·v satisfies dist_sq p r = dist_sq p q + dist_sq q r + 2·dot, and    *)
(* Cauchy–Schwarz gives dot ≤ √(dist_sq p q · dist_sq q r); squaring the       *)
(* (non-negative) target reduces the goal to that bound.                       *)
Theorem dist_triangle : forall p q r,
  dist p r <= dist p q + dist q r.
Proof.
  intros p q r.
  pose proof (dist_sq_nonneg p q) as HX.
  pose proof (dist_sq_nonneg q r) as HY.
  pose proof (dist_sq_nonneg p r) as HZ.
  set (dot := (px p - px q) * (px q - px r) + (py p - py q) * (py q - py r)).
  assert (Hdec : dist_sq p r = dist_sq p q + dist_sq q r + 2 * dot)
    by (unfold dist_sq, dot; ring).
  assert (HcsS : dot * dot <= dist_sq p q * dist_sq q r).
  { unfold dot, dist_sq.
    pose proof (sqr_nonneg ((px p - px q) * (py q - py r)
                            - (py p - py q) * (px q - px r))). nra. }
  pose proof (sqrt_pos (dist_sq p q * dist_sq q r)) as HSpos.
  assert (HSsq : sqrt (dist_sq p q * dist_sq q r)
                 * sqrt (dist_sq p q * dist_sq q r)
                 = dist_sq p q * dist_sq q r).
  { apply sqrt_sqrt. apply Rmult_le_pos; assumption. }
  assert (Hdot_le : dot <= sqrt (dist_sq p q * dist_sq q r)).
  { destruct (Rle_or_lt dot 0) as [Hle | Hgt].
    - lra.
    - apply (proj2 (sq_monotone_nonneg dot
               (sqrt (dist_sq p q * dist_sq q r)) (Rlt_le _ _ Hgt) HSpos)).
      rewrite HSsq. exact HcsS. }
  unfold dist.
  pose proof (sqrt_pos (dist_sq p q)) as HsX.
  pose proof (sqrt_pos (dist_sq q r)) as HsY.
  apply (proj2 (sq_monotone_nonneg (sqrt (dist_sq p r))
            (sqrt (dist_sq p q) + sqrt (dist_sq q r))
            (sqrt_pos _) (Rplus_le_le_0_compat _ _ HsX HsY))).
  rewrite (sqrt_sqrt (dist_sq p r) HZ).
  replace ((sqrt (dist_sq p q) + sqrt (dist_sq q r))
           * (sqrt (dist_sq p q) + sqrt (dist_sq q r)))
    with (sqrt (dist_sq p q) * sqrt (dist_sq p q)
          + sqrt (dist_sq q r) * sqrt (dist_sq q r)
          + 2 * (sqrt (dist_sq p q) * sqrt (dist_sq q r))) by ring.
  rewrite (sqrt_sqrt (dist_sq p q) HX).
  rewrite (sqrt_sqrt (dist_sq q r) HY).
  rewrite <- (sqrt_mult (dist_sq p q) (dist_sq q r) HX HY).
  rewrite Hdec. lra.
Qed.

(* -------------------------------------------------------------------------- *)
(* Assumption audit. The proofs above rely only on the constructions of the   *)
(* standard library's classical real arithmetic.  Run with `make` or          *)
(* `rocq compile theories/Distance.v` and inspect the `Print Assumptions`     *)
(* output: any axiom not part of the classical Reals stdlib is a red flag.    *)
(* -------------------------------------------------------------------------- *)

Print Assumptions dist_le_iff_dist_sq_le.
Print Assumptions dist_triangle.
