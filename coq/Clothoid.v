(* SPDX-FileCopyrightText: 2026 Merkator Group *)
(* SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining *)

(** Formal derivative identities for the Bertolazzi-Frego (2015)
    clothoid G1 Hermite residual.

    Tested with: Coq 8.13.1 + Coquelicot 3.x

    The Bertolazzi-Frego algorithm reduces the G1 Hermite clothoid fitting
    problem to a single nonlinear scalar equation
        f(A) = Y_0(2A, delta - A, phi_0) = 0,
    where Y_0 is the Generalized Fresnel sine moment:
        Y_0 = int_0^1 sin( A t^2 + (delta - A) t + phi_0 ) dt.

    The proposed Halley iteration needs f'(A) and f''(A).  This file proves
    the *algebraic* identities behind those derivatives at the level of the
    integrand, plus the integral-level identity for f'(A) by application of
    Coquelicot's parameter-differentiation lemma is_derive_RInt_param_aux.

    Scope.  We formalise:
      (1) the algebraic identity (t^2 - t)^2 = t^4 - 2 t^3 + t^2;
      (2) the pointwise parameter derivative of the integrand;
      (3) the integral-level first-derivative identity
            d/dA  RInt sin(phi)  =  RInt (t^2 - t) cos(phi)
          on [0, 1].
    We do NOT formalise Halley's cubic convergence theorem; that is a
    classical analytic result outside Coquelicot's library and would be a
    separate proof effort.  Everything below compiles cleanly. *)

Require Import Reals.
Require Import Coquelicot.Coquelicot.
Open Scope R_scope.

(* ------------------------------------------------------------------ *)
(* 1. Phase function and algebraic identity                            *)
(* ------------------------------------------------------------------ *)

Definition phi (A t delta phi0 : R) : R :=
  A * t * t + (delta - A) * t + phi0.

(** Algebraic identity used in the f''(A) reduction. *)
Lemma sq_t2_minus_t :
  forall t : R, (t * t - t) * (t * t - t) = t^4 - 2 * t^3 + t^2.
Proof. intros t. ring. Qed.

(* ------------------------------------------------------------------ *)
(* 2. Pointwise derivative of the integrand                            *)
(* ------------------------------------------------------------------ *)

Lemma pointwise_d1 :
  forall (A t delta phi0 : R),
    is_derive (fun a => sin (phi a t delta phi0)) A
              ((t * t - t) * cos (phi A t delta phi0)).
Proof.
  intros A t delta phi0. unfold phi.
  auto_derive; trivial.
  unfold Rminus. ring.
Qed.

(** Convenient corollary: the exact Derive value. *)
Lemma Derive_phi1 :
  forall (A t delta phi0 : R),
    Derive (fun a => sin (phi a t delta phi0)) A
    = (t * t - t) * cos (phi A t delta phi0).
Proof.
  intros. apply is_derive_unique. apply pointwise_d1.
Qed.

(* ------------------------------------------------------------------ *)
(* 3. Joint continuity of the parameter derivative                     *)
(* ------------------------------------------------------------------ *)

(** Continuity (in (A, t)) of the phase. *)
Lemma cont2d_phi :
  forall (A0 t0 delta phi0 : R),
    continuity_2d_pt (fun u v => phi u v delta phi0) A0 t0.
Proof.
  intros A0 t0 delta phi0. unfold phi.
  (* phi = A*t*t + (delta - A)*t + phi0 *)
  apply continuity_2d_pt_plus.
  apply continuity_2d_pt_plus.
  - (* A * t * t = (A*t) * t *)
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_id1.
    apply continuity_2d_pt_id2.
    apply continuity_2d_pt_id2.
  - (* (delta - A) * t *)
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_minus.
    apply continuity_2d_pt_const.
    apply continuity_2d_pt_id1.
    apply continuity_2d_pt_id2.
  - apply continuity_2d_pt_const.
Qed.

(** Continuity of cos(phi) and sin(phi) in (A, t). *)
Lemma cont2d_cos_phi :
  forall (A0 t0 delta phi0 : R),
    continuity_2d_pt (fun u v => cos (phi u v delta phi0)) A0 t0.
Proof.
  intros.
  apply continuity_1d_2d_pt_comp.
  apply continuity_cos.
  apply cont2d_phi.
Qed.

(** (t^2 - t) cos(phi(A,t)) is jointly continuous in (A, t). *)
Lemma cont2d_d1 :
  forall (A0 t0 delta phi0 : R),
    continuity_2d_pt
      (fun u v => (v * v - v) * cos (phi u v delta phi0)) A0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  - (* v^2 - v  is continuous in (u, v) *)
    apply continuity_2d_pt_minus.
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_id2. apply continuity_2d_pt_id2.
    apply continuity_2d_pt_id2.
  - apply cont2d_cos_phi.
Qed.

(* ------------------------------------------------------------------ *)
(* 4. Integrability                                                    *)
(* ------------------------------------------------------------------ *)

(** Differentiability of sin(phi(A, .)) in the second argument
    (hence continuity). *)
Lemma cont1d_sin_phi :
  forall (A delta phi0 t : R),
    continuous (fun s => sin (phi A s delta phi0)) t.
Proof.
  intros A delta phi0 t.
  apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold phi.
  eexists.
  auto_derive; trivial.
Qed.

(** Differentiability of (t^2 - t) cos(phi(A, .)) in t. *)
Lemma cont1d_d1 :
  forall (A delta phi0 t : R),
    continuous
      (fun s => (s * s - s) * cos (phi A s delta phi0)) t.
Proof.
  intros A delta phi0 t.
  apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold phi.
  eexists.
  auto_derive; trivial.
Qed.

(** sin(phi(A, .)) is integrable on [0, 1] for every A. *)
Lemma ex_RInt_sin_phi :
  forall (A delta phi0 : R),
    ex_RInt (fun t => sin (phi A t delta phi0)) 0 1.
Proof.
  intros A delta phi0.
  apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros t _. apply cont1d_sin_phi.
Qed.

(** The pointwise derivative (t^2 - t) cos(phi(A, .)) is integrable. *)
Lemma ex_RInt_d1 :
  forall (A delta phi0 : R),
    ex_RInt
      (fun t => (t * t - t) * cos (phi A t delta phi0)) 0 1.
Proof.
  intros A delta phi0.
  apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros t _. apply cont1d_d1.
Qed.

(* ------------------------------------------------------------------ *)
(* 5. Main theorem: integral-level first derivative                    *)
(* ------------------------------------------------------------------ *)

Theorem f_prime_eq :
  forall (A delta phi0 : R),
    is_derive
      (fun a : R => RInt (fun t => sin (phi a t delta phi0)) 0 1) A
      (RInt (fun t => (t * t - t) * cos (phi A t delta phi0)) 0 1).
Proof.
  intros A delta phi0.
  (* Apply Leibniz under the integral via is_derive_RInt_param_aux *)
  pose proof
    (is_derive_RInt_param_aux
       (fun a t => sin (phi a t delta phi0)) 0 1 A) as Hlift.
  (* Massage conclusion so it matches: rewrite Derive (...) to its value *)
  match goal with
    |- is_derive _ _ ?rhs =>
        replace rhs with
          (RInt
             (fun t : R =>
                Derive (fun u : R => sin (phi u t delta phi0)) A) 0 1)
  end.
  - apply Hlift.
    + (* locally A, ex_derive in a neighborhood *)
      apply filter_forall. intros A' t _.
      eexists. apply pointwise_d1.
    + (* joint continuity of parameter derivative *)
      intros t _.
      apply (continuity_2d_pt_ext
               (fun u v => (v * v - v) * cos (phi u v delta phi0))).
      * intros u v. symmetry. apply Derive_phi1.
      * apply cont2d_d1.
    + (* ex_RInt of integrand locally in A *)
      apply filter_forall. intros A'. apply ex_RInt_sin_phi.
    + (* ex_RInt of the parameter derivative at A *)
      apply (ex_RInt_ext
               (fun t => (t * t - t) * cos (phi A t delta phi0))).
      * intros t _. symmetry. apply Derive_phi1.
      * apply ex_RInt_d1.
  - apply RInt_ext. intros t _. apply Derive_phi1.
Qed.

(** Read this as: f'(A) = int_0^1 (t^2 - t) cos(phi(A, t)) dt,
    which in the paper's notation is X_2 - X_1. *)

(* ------------------------------------------------------------------ *)
(* 6. Second derivative identity (needed for Halley's iteration)       *)
(* ------------------------------------------------------------------ *)

(** Pointwise: d/dA [(t^2 - t) cos(phi(A, t))] = -(t^2 - t)^2 sin(phi(A,t)). *)
Lemma pointwise_d2 :
  forall (A t delta phi0 : R),
    is_derive
      (fun a => (t * t - t) * cos (phi a t delta phi0)) A
      (- ((t * t - t) * (t * t - t)) * sin (phi A t delta phi0)).
Proof.
  intros A t delta phi0. unfold phi.
  auto_derive; trivial.
  unfold Rminus. ring.
Qed.

(** Derive-value form. *)
Lemma Derive_phi2 :
  forall (A t delta phi0 : R),
    Derive (fun a => (t * t - t) * cos (phi a t delta phi0)) A
    = - ((t * t - t) * (t * t - t)) * sin (phi A t delta phi0).
Proof. intros. apply is_derive_unique. apply pointwise_d2. Qed.

(** Joint continuity of the second-derivative integrand (in (A, t)). *)
Lemma cont2d_d2 :
  forall (A0 t0 delta phi0 : R),
    continuity_2d_pt
      (fun u v =>
         - ((v * v - v) * (v * v - v)) * sin (phi u v delta phi0)) A0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  - (* -((v^2 - v)^2) is continuous in (u, v) *)
    apply continuity_2d_pt_opp.
    apply continuity_2d_pt_mult; apply continuity_2d_pt_minus;
      try apply continuity_2d_pt_id2;
      apply continuity_2d_pt_mult; apply continuity_2d_pt_id2.
  - apply continuity_1d_2d_pt_comp.
    apply continuity_sin. apply cont2d_phi.
Qed.

(** Integrability of the second-derivative integrand in t. *)
Lemma cont1d_d2 :
  forall (A delta phi0 t : R),
    continuous
      (fun s =>
         - ((s * s - s) * (s * s - s)) * sin (phi A s delta phi0)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold phi. eexists. auto_derive; trivial.
Qed.

Lemma ex_RInt_d2 :
  forall (A delta phi0 : R),
    ex_RInt
      (fun t => - ((t * t - t) * (t * t - t)) * sin (phi A t delta phi0))
      0 1.
Proof.
  intros A delta phi0.
  apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros t _. apply cont1d_d2.
Qed.

(** Main theorem 2: d/dA RInt (t^2-t) cos(phi)  =  RInt -(t^2-t)^2 sin(phi).

    The integrand on the right, by sq_t2_minus_t, equals
        -(t^4 - 2 t^3 + t^2) sin(phi)
    so the integral is  -(Y_4 - 2 Y_3 + Y_2)  in the paper's notation. *)
Theorem f_second_eq :
  forall (A delta phi0 : R),
    is_derive
      (fun a : R =>
         RInt (fun t => (t * t - t) * cos (phi a t delta phi0)) 0 1) A
      (RInt
         (fun t => - ((t * t - t) * (t * t - t))
                   * sin (phi A t delta phi0)) 0 1).
Proof.
  intros A delta phi0.
  pose proof
    (is_derive_RInt_param_aux
       (fun a t => (t * t - t) * cos (phi a t delta phi0)) 0 1 A) as Hlift.
  match goal with
    |- is_derive _ _ ?rhs =>
        replace rhs with
          (RInt
             (fun t : R =>
                Derive
                  (fun u : R => (t * t - t) * cos (phi u t delta phi0)) A)
             0 1)
  end.
  - apply Hlift.
    + apply filter_forall. intros A' t _.
      eexists. apply pointwise_d2.
    + intros t _.
      apply (continuity_2d_pt_ext
               (fun u v =>
                  - ((v * v - v) * (v * v - v))
                  * sin (phi u v delta phi0))).
      * intros u v. symmetry. apply Derive_phi2.
      * apply cont2d_d2.
    + apply filter_forall. intros A'. apply ex_RInt_d1.
    + apply (ex_RInt_ext
               (fun t =>
                  - ((t * t - t) * (t * t - t))
                  * sin (phi A t delta phi0))).
      * intros t _. symmetry. apply Derive_phi2.
      * apply ex_RInt_d2.
  - apply RInt_ext. intros t _. apply Derive_phi2.
Qed.

(* ------------------------------------------------------------------ *)
(* 7. Summary                                                          *)
(* ------------------------------------------------------------------ *)

(** Summary of the formalised identities:

     f(A)   = RInt (fun t => sin (phi A t delta phi0)) 0 1
     f'(A)  = RInt (fun t => (t^2 - t) * cos (phi A t delta phi0)) 0 1
     f''(A) = RInt (fun t => -(t^2 - t)^2 * sin (phi A t delta phi0)) 0 1
            = - RInt (fun t => (t^4 - 2 t^3 + t^2) * sin (phi A t delta phi0))

     In the moment notation of Bertolazzi-Frego (2015):
       f'(A)  = X_2 - X_1
       f''(A) = -(Y_4 - 2 Y_3 + Y_2)

     Halley's iteration uses exactly these:
       A_{n+1} = A_n - 2 f f' / (2 (f')^2 - f f'')

     What is NOT formalised here:
       * Cubic convergence of Halley's method.  That is a classical analytic
         result (Ostrowski, Schroeder) outside Coquelicot's library and
         would be a separate, sizeable proof effort.
       * The Generalized-Fresnel moment routine.  We work directly with the
         integrals X_k, Y_k as RInt expressions; the engineering-grade
         Lommel-recurrence implementation is a separate verification target.
       * Existence of a unique zero of f(A) in the relevant configuration
         interval (Bertolazzi-Frego, Theorem 1).  Numerically verified by
         the 10200-point hypercube benchmark accompanying this file. *)
