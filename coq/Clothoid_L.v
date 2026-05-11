(* SPDX-FileCopyrightText: 2026 Merkator Group *)
(* SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining *)

(** Formal derivative identities for the chord-length clothoid residual.

    Given fixed (P0, P1, kappa_0, kappa_1) with d = |P1 - P0|, the
    chord-length residual viewed as a function of the unknown arc length L
    is
        f(L) = L^2 * (P(L)^2 + Q(L)^2) - d^2
    where, with psi(tau) = kappa_0 * tau + (kappa_1 - kappa_0) * tau^2 / 2,
        P(L) = int_0^1 cos(L * psi(tau)) dtau
        Q(L) = int_0^1 sin(L * psi(tau)) dtau.

    This file proves the derivative identities used by Halley-on-L:
        P'(L)  = - T(L)
        Q'(L)  =   R(L)
        R'(L)  = - S2s(L)
        T'(L)  =   S2c(L)
        f'(L)  = 2 L (P^2 + Q^2) + 2 L^2 (Q R - P T)
        f''(L) = 2 (P^2 + Q^2) + 8 L (Q R - P T)
                 + 2 L^2 (R^2 + T^2 - P S2c - Q S2s)
    with R, T, S2c, S2s the first- and second-moment integrals.

    Tested with Coq 8.13.1 / 8.20.1 + Coquelicot 3.x. *)

Require Import Reals.
Require Import Coquelicot.Coquelicot.
Open Scope R_scope.

(* ==================================================================== *)
(* 1. Definitions                                                       *)
(* ==================================================================== *)

Definition psi (tau k0 k1 : R) : R :=
  k0 * tau + (k1 - k0) * tau * tau / 2.

Definition P_int (L k0 k1 : R) : R :=
  RInt (fun tau => cos (L * psi tau k0 k1)) 0 1.

Definition Q_int (L k0 k1 : R) : R :=
  RInt (fun tau => sin (L * psi tau k0 k1)) 0 1.

Definition R_int (L k0 k1 : R) : R :=
  RInt (fun tau => psi tau k0 k1 * cos (L * psi tau k0 k1)) 0 1.

Definition T_int (L k0 k1 : R) : R :=
  RInt (fun tau => psi tau k0 k1 * sin (L * psi tau k0 k1)) 0 1.

Definition S2c_int (L k0 k1 : R) : R :=
  RInt (fun tau =>
          psi tau k0 k1 * psi tau k0 k1 * cos (L * psi tau k0 k1)) 0 1.

Definition S2s_int (L k0 k1 : R) : R :=
  RInt (fun tau =>
          psi tau k0 k1 * psi tau k0 k1 * sin (L * psi tau k0 k1)) 0 1.

(** The chord-length residual.  d2 is the squared chord length.

    Written with Coquelicot's [minus] / [plus] / [mult] / [opp] on R so the
    is_derive_minus / is_derive_plus / is_derive_mult composition lemmas
    apply by syntactic unification (avoids Rminus vs minus unification
    failures).  Convertible to the usual a*b - c form on R. *)
Definition fL (L k0 k1 d2 : R) : R :=
  minus (mult (mult L L)
              (plus (mult (P_int L k0 k1) (P_int L k0 k1))
                    (mult (Q_int L k0 k1) (Q_int L k0 k1))))
        d2.

(* ==================================================================== *)
(* 2. Pointwise derivatives of the integrands w.r.t. L                  *)
(* ==================================================================== *)

(** d/dL cos(L * psi) = - psi * sin(L * psi). *)
Lemma dL_cos :
  forall (tau L k0 k1 : R),
    is_derive (fun u => cos (u * psi tau k0 k1)) L
              (- psi tau k0 k1 * sin (L * psi tau k0 k1)).
Proof.
  intros. auto_derive; trivial. ring.
Qed.

(** d/dL sin(L * psi) = psi * cos(L * psi). *)
Lemma dL_sin :
  forall (tau L k0 k1 : R),
    is_derive (fun u => sin (u * psi tau k0 k1)) L
              (psi tau k0 k1 * cos (L * psi tau k0 k1)).
Proof.
  intros. auto_derive; trivial. ring.
Qed.

(** d/dL psi * cos(L * psi) = - psi^2 * sin(L * psi). *)
Lemma dL_psi_cos :
  forall (tau L k0 k1 : R),
    is_derive (fun u => psi tau k0 k1 * cos (u * psi tau k0 k1)) L
              (- (psi tau k0 k1 * psi tau k0 k1)
               * sin (L * psi tau k0 k1)).
Proof.
  intros. auto_derive; trivial. ring.
Qed.

(** d/dL psi * sin(L * psi) = psi^2 * cos(L * psi). *)
Lemma dL_psi_sin :
  forall (tau L k0 k1 : R),
    is_derive (fun u => psi tau k0 k1 * sin (u * psi tau k0 k1)) L
              ((psi tau k0 k1 * psi tau k0 k1)
               * cos (L * psi tau k0 k1)).
Proof.
  intros. auto_derive; trivial. ring.
Qed.

(** Derive-value forms (used to discharge `Derive` in is_derive_RInt_param_aux). *)
Lemma Derive_dL_cos : forall tau L k0 k1,
  Derive (fun u => cos (u * psi tau k0 k1)) L
  = - psi tau k0 k1 * sin (L * psi tau k0 k1).
Proof. intros. apply is_derive_unique. apply dL_cos. Qed.

Lemma Derive_dL_sin : forall tau L k0 k1,
  Derive (fun u => sin (u * psi tau k0 k1)) L
  = psi tau k0 k1 * cos (L * psi tau k0 k1).
Proof. intros. apply is_derive_unique. apply dL_sin. Qed.

Lemma Derive_dL_psi_cos : forall tau L k0 k1,
  Derive (fun u => psi tau k0 k1 * cos (u * psi tau k0 k1)) L
  = - (psi tau k0 k1 * psi tau k0 k1) * sin (L * psi tau k0 k1).
Proof. intros. apply is_derive_unique. apply dL_psi_cos. Qed.

Lemma Derive_dL_psi_sin : forall tau L k0 k1,
  Derive (fun u => psi tau k0 k1 * sin (u * psi tau k0 k1)) L
  = (psi tau k0 k1 * psi tau k0 k1) * cos (L * psi tau k0 k1).
Proof. intros. apply is_derive_unique. apply dL_psi_sin. Qed.

(* ==================================================================== *)
(* 3. Joint continuity in (L, tau)                                       *)
(* ==================================================================== *)

(** psi is jointly continuous in (L, tau).  (Only tau matters, but we use
    the 2D form because is_derive_RInt_param_aux wants it.) *)
Lemma cont2d_psi :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt (fun _ v => psi v k0 k1) L0 t0.
Proof.
  intros L0 t0 k0 k1. unfold psi.
  apply continuity_2d_pt_plus.
  - apply continuity_2d_pt_mult.
    apply continuity_2d_pt_const. apply continuity_2d_pt_id2.
  - (* (k1 - k0) * v * v / 2  -- written as (...) * /2 *)
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_mult.
    apply continuity_2d_pt_const.
    apply continuity_2d_pt_id2.
    apply continuity_2d_pt_id2.
    apply continuity_2d_pt_const.
Qed.

(** L * psi is jointly continuous. *)
Lemma cont2d_Lpsi :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt (fun u v => u * psi v k0 k1) L0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  apply continuity_2d_pt_id1. apply cont2d_psi.
Qed.

(** cos(L * psi) and sin(L * psi) are jointly continuous. *)
Lemma cont2d_cos_Lpsi :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt (fun u v => cos (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_1d_2d_pt_comp.
  apply continuity_cos. apply cont2d_Lpsi.
Qed.

Lemma cont2d_sin_Lpsi :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt (fun u v => sin (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_1d_2d_pt_comp.
  apply continuity_sin. apply cont2d_Lpsi.
Qed.

(** The four parameter-derivative integrands are jointly continuous. *)
Lemma cont2d_neg_psi_sin :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt
      (fun u v => - psi v k0 k1 * sin (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  apply continuity_2d_pt_opp. apply cont2d_psi.
  apply cont2d_sin_Lpsi.
Qed.

Lemma cont2d_psi_cos :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt
      (fun u v => psi v k0 k1 * cos (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  apply cont2d_psi. apply cont2d_cos_Lpsi.
Qed.

Lemma cont2d_psi_sin :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt
      (fun u v => psi v k0 k1 * sin (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  apply cont2d_psi. apply cont2d_sin_Lpsi.
Qed.

Lemma cont2d_neg_psi2_sin :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt
      (fun u v =>
         - (psi v k0 k1 * psi v k0 k1) * sin (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  apply continuity_2d_pt_opp.
  apply continuity_2d_pt_mult. apply cont2d_psi. apply cont2d_psi.
  apply cont2d_sin_Lpsi.
Qed.

Lemma cont2d_psi2_cos :
  forall (L0 t0 k0 k1 : R),
    continuity_2d_pt
      (fun u v =>
         (psi v k0 k1 * psi v k0 k1) * cos (u * psi v k0 k1)) L0 t0.
Proof.
  intros. apply continuity_2d_pt_mult.
  apply continuity_2d_pt_mult. apply cont2d_psi. apply cont2d_psi.
  apply cont2d_cos_Lpsi.
Qed.

(* ==================================================================== *)
(* 4. Integrability                                                      *)
(* ==================================================================== *)

(** 1-D continuity in tau (at fixed L). Each integrand we need. *)
Lemma cont1d_cos_Lpsi : forall L k0 k1 t,
  continuous (fun tau => cos (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

Lemma cont1d_sin_Lpsi : forall L k0 k1 t,
  continuous (fun tau => sin (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

Lemma cont1d_psi_cos : forall L k0 k1 t,
  continuous (fun tau => psi tau k0 k1 * cos (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

Lemma cont1d_psi_sin : forall L k0 k1 t,
  continuous (fun tau => psi tau k0 k1 * sin (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

Lemma cont1d_neg_psi_sin : forall L k0 k1 t,
  continuous (fun tau => - psi tau k0 k1 * sin (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

Lemma cont1d_psi2_cos : forall L k0 k1 t,
  continuous
    (fun tau => (psi tau k0 k1 * psi tau k0 k1) * cos (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

Lemma cont1d_neg_psi2_sin : forall L k0 k1 t,
  continuous
    (fun tau =>
       - (psi tau k0 k1 * psi tau k0 k1) * sin (L * psi tau k0 k1)) t.
Proof.
  intros. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
  unfold psi. eexists. auto_derive; trivial.
Qed.

(** Integrability on [0, 1]. *)
Lemma ex_RInt_cos_Lpsi : forall L k0 k1,
  ex_RInt (fun tau => cos (L * psi tau k0 k1)) 0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_cos_Lpsi.
Qed.

Lemma ex_RInt_sin_Lpsi : forall L k0 k1,
  ex_RInt (fun tau => sin (L * psi tau k0 k1)) 0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_sin_Lpsi.
Qed.

Lemma ex_RInt_psi_cos : forall L k0 k1,
  ex_RInt (fun tau => psi tau k0 k1 * cos (L * psi tau k0 k1)) 0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_psi_cos.
Qed.

Lemma ex_RInt_psi_sin : forall L k0 k1,
  ex_RInt (fun tau => psi tau k0 k1 * sin (L * psi tau k0 k1)) 0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_psi_sin.
Qed.

Lemma ex_RInt_neg_psi_sin : forall L k0 k1,
  ex_RInt (fun tau => - psi tau k0 k1 * sin (L * psi tau k0 k1)) 0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_neg_psi_sin.
Qed.

Lemma ex_RInt_psi2_cos : forall L k0 k1,
  ex_RInt
    (fun tau => (psi tau k0 k1 * psi tau k0 k1) * cos (L * psi tau k0 k1))
    0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_psi2_cos.
Qed.

Lemma ex_RInt_neg_psi2_sin : forall L k0 k1,
  ex_RInt
    (fun tau =>
       - (psi tau k0 k1 * psi tau k0 k1) * sin (L * psi tau k0 k1))
    0 1.
Proof.
  intros. apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
  intros tau _. apply cont1d_neg_psi2_sin.
Qed.

(* ==================================================================== *)
(* 5. Integral-level derivatives                                         *)
(* ==================================================================== *)

(** Helper: RInt of negated R-valued function. *)
Lemma RInt_Ropp_R :
  forall (f : R -> R) (a b : R),
    ex_RInt f a b ->
    RInt (fun x => - f x) a b = - RInt f a b.
Proof.
  intros f a b H.
  change (- RInt f a b) with (opp (RInt f a b)).
  rewrite <- (RInt_opp f a b H).
  apply RInt_ext. intros. reflexivity.
Qed.

(** P'(L) = - T(L). *)
Theorem is_derive_P :
  forall (L k0 k1 : R),
    is_derive (fun u => P_int u k0 k1) L (- T_int L k0 k1).
Proof.
  intros L k0 k1.
  unfold P_int, T_int.
  pose proof
    (is_derive_RInt_param_aux
       (fun u tau => cos (u * psi tau k0 k1)) 0 1 L) as Hlift.
  (* Massage RHS *)
  replace (- RInt (fun tau => psi tau k0 k1 * sin (L * psi tau k0 k1)) 0 1)
     with (RInt
             (fun tau =>
                Derive (fun u => cos (u * psi tau k0 k1)) L) 0 1).
  - apply Hlift.
    + apply filter_forall. intros L' tau _.
      eexists. apply dL_cos.
    + intros tau _.
      apply (continuity_2d_pt_ext
               (fun u v => - psi v k0 k1 * sin (u * psi v k0 k1))).
      * intros u v. symmetry. apply Derive_dL_cos.
      * apply cont2d_neg_psi_sin.
    + apply filter_forall. intros L'. apply ex_RInt_cos_Lpsi.
    + apply (ex_RInt_ext
               (fun tau => - psi tau k0 k1 * sin (L * psi tau k0 k1))).
      * intros tau _. symmetry. apply Derive_dL_cos.
      * apply ex_RInt_neg_psi_sin.
  - transitivity
      (RInt (fun tau => - (psi tau k0 k1 * sin (L * psi tau k0 k1))) 0 1).
    + apply RInt_ext. intros tau _. rewrite Derive_dL_cos.
      apply Ropp_mult_distr_l_reverse.
    + apply RInt_Ropp_R. apply ex_RInt_psi_sin.
Qed.

(** Q'(L) = R(L). *)
Theorem is_derive_Q :
  forall (L k0 k1 : R),
    is_derive (fun u => Q_int u k0 k1) L (R_int L k0 k1).
Proof.
  intros L k0 k1.
  unfold Q_int, R_int.
  pose proof
    (is_derive_RInt_param_aux
       (fun u tau => sin (u * psi tau k0 k1)) 0 1 L) as Hlift.
  replace (RInt (fun tau => psi tau k0 k1 * cos (L * psi tau k0 k1)) 0 1)
     with (RInt
             (fun tau =>
                Derive (fun u => sin (u * psi tau k0 k1)) L) 0 1).
  - apply Hlift.
    + apply filter_forall. intros L' tau _.
      eexists. apply dL_sin.
    + intros tau _.
      apply (continuity_2d_pt_ext
               (fun u v => psi v k0 k1 * cos (u * psi v k0 k1))).
      * intros u v. symmetry. apply Derive_dL_sin.
      * apply cont2d_psi_cos.
    + apply filter_forall. intros L'. apply ex_RInt_sin_Lpsi.
    + apply (ex_RInt_ext
               (fun tau => psi tau k0 k1 * cos (L * psi tau k0 k1))).
      * intros tau _. symmetry. apply Derive_dL_sin.
      * apply ex_RInt_psi_cos.
  - apply RInt_ext. intros tau _. apply Derive_dL_sin.
Qed.

(** R'(L) = - S2s(L). *)
Theorem is_derive_R :
  forall (L k0 k1 : R),
    is_derive (fun u => R_int u k0 k1) L (- S2s_int L k0 k1).
Proof.
  intros L k0 k1.
  unfold R_int, S2s_int.
  pose proof
    (is_derive_RInt_param_aux
       (fun u tau => psi tau k0 k1 * cos (u * psi tau k0 k1)) 0 1 L) as Hlift.
  replace
    (- RInt
        (fun tau =>
           psi tau k0 k1 * psi tau k0 k1 * sin (L * psi tau k0 k1)) 0 1)
    with
      (RInt
         (fun tau =>
            Derive (fun u => psi tau k0 k1 * cos (u * psi tau k0 k1)) L)
         0 1).
  - apply Hlift.
    + apply filter_forall. intros L' tau _.
      eexists. apply dL_psi_cos.
    + intros tau _.
      apply (continuity_2d_pt_ext
               (fun u v =>
                  - (psi v k0 k1 * psi v k0 k1) * sin (u * psi v k0 k1))).
      * intros u v. symmetry. apply Derive_dL_psi_cos.
      * apply cont2d_neg_psi2_sin.
    + apply filter_forall. intros L'. apply ex_RInt_psi_cos.
    + apply (ex_RInt_ext
               (fun tau =>
                  - (psi tau k0 k1 * psi tau k0 k1)
                  * sin (L * psi tau k0 k1))).
      * intros tau _. symmetry. apply Derive_dL_psi_cos.
      * apply ex_RInt_neg_psi2_sin.
  - transitivity
      (RInt
         (fun tau =>
            - (psi tau k0 k1 * psi tau k0 k1 * sin (L * psi tau k0 k1)))
         0 1).
    + apply RInt_ext. intros tau _.
      rewrite Derive_dL_psi_cos.
      apply Ropp_mult_distr_l_reverse.
    + apply RInt_Ropp_R.
      apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
      intros tau _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
      unfold psi. eexists. auto_derive; trivial.
Qed.

(** Convenience: (P_int L)' as is_derive in the precise form we want. *)
Lemma is_derive_P_at :
  forall (L k0 k1 : R),
    is_derive (fun u => P_int u k0 k1) L (- T_int L k0 k1).
Proof. apply is_derive_P. Qed.

Lemma is_derive_Q_at :
  forall (L k0 k1 : R),
    is_derive (fun u => Q_int u k0 k1) L (R_int L k0 k1).
Proof. apply is_derive_Q. Qed.

(** T'(L) = S2c(L). *)
Theorem is_derive_T :
  forall (L k0 k1 : R),
    is_derive (fun u => T_int u k0 k1) L (S2c_int L k0 k1).
Proof.
  intros L k0 k1.
  unfold T_int, S2c_int.
  pose proof
    (is_derive_RInt_param_aux
       (fun u tau => psi tau k0 k1 * sin (u * psi tau k0 k1)) 0 1 L) as Hlift.
  replace
    (RInt
       (fun tau =>
          psi tau k0 k1 * psi tau k0 k1 * cos (L * psi tau k0 k1)) 0 1)
    with
      (RInt
         (fun tau =>
            Derive (fun u => psi tau k0 k1 * sin (u * psi tau k0 k1)) L)
         0 1).
  - apply Hlift.
    + apply filter_forall. intros L' tau _.
      eexists. apply dL_psi_sin.
    + intros tau _.
      apply (continuity_2d_pt_ext
               (fun u v =>
                  (psi v k0 k1 * psi v k0 k1) * cos (u * psi v k0 k1))).
      * intros u v. symmetry. apply Derive_dL_psi_sin.
      * apply cont2d_psi2_cos.
    + apply filter_forall. intros L'. apply ex_RInt_psi_sin.
    + apply (ex_RInt_ext
               (fun tau =>
                  (psi tau k0 k1 * psi tau k0 k1)
                  * cos (L * psi tau k0 k1))).
      * intros tau _. symmetry. apply Derive_dL_psi_sin.
      * apply ex_RInt_psi2_cos.
  - apply RInt_ext. intros tau _. apply Derive_dL_psi_sin.
Qed.

(* ==================================================================== *)
(* 6. f'(L): the first derivative of the chord-length residual           *)
(* ==================================================================== *)

(** Compute Derive of P^2 + Q^2 + (L^2)*(P^2+Q^2) - d^2 the easy way:
    use auto_derive with the four basic integral derivatives as Hint
    Extern rules.  This avoids the manual is_derive_plus / is_derive_mult
    composition pain (Rplus vs plus unification). *)

(** Derive of P_int as a function of L. *)
Lemma Derive_P_int : forall L k0 k1,
  Derive (fun u => P_int u k0 k1) L = - T_int L k0 k1.
Proof. intros. apply is_derive_unique. apply is_derive_P. Qed.

Lemma Derive_Q_int : forall L k0 k1,
  Derive (fun u => Q_int u k0 k1) L = R_int L k0 k1.
Proof. intros. apply is_derive_unique. apply is_derive_Q. Qed.

(** ex_derive versions, for auto_derive to discharge. *)
Lemma ex_derive_P_int : forall L k0 k1, ex_derive (fun u => P_int u k0 k1) L.
Proof. intros. eexists. apply is_derive_P. Qed.

Lemma ex_derive_Q_int : forall L k0 k1, ex_derive (fun u => Q_int u k0 k1) L.
Proof. intros. eexists. apply is_derive_Q. Qed.

(** Now we register these as hints for auto_derive.
    Coquelicot's auto_derive uses the [derive] hint database; adding
    ex_derive instances + Derive_… computations is the recommended way. *)

(* (Hint setup not strictly needed for the Derive-based proof below.) *)

(** First main theorem (raw form, written with minus/plus/mult to match
    the conclusion shape of is_derive_minus/plus/mult exactly). *)
Theorem f_L_prime_eq_raw :
  forall (L k0 k1 d2 : R),
    is_derive (fun u => fL u k0 k1 d2) L
              (minus
                 (plus (mult (plus (mult 1 L) (mult L 1))
                             (plus (mult (P_int L k0 k1) (P_int L k0 k1))
                                   (mult (Q_int L k0 k1) (Q_int L k0 k1))))
                       (mult (mult L L)
                             (plus
                                (plus (mult (- T_int L k0 k1) (P_int L k0 k1))
                                      (mult (P_int L k0 k1) (- T_int L k0 k1)))
                                (plus (mult (R_int L k0 k1) (Q_int L k0 k1))
                                      (mult (Q_int L k0 k1) (R_int L k0 k1))))))
                 0).
Proof.
  intros L k0 k1 d2. unfold fL.
  apply (@is_derive_minus R_AbsRing R_NormedModule
           (fun u : R => mult (mult u u)
                              (plus (mult (P_int u k0 k1) (P_int u k0 k1))
                                    (mult (Q_int u k0 k1) (Q_int u k0 k1))))
           (fun _ : R => d2)
           L).
  - apply (@is_derive_mult R_AbsRing
             (fun u : R => mult u u)
             (fun u : R => plus (mult (P_int u k0 k1) (P_int u k0 k1))
                                (mult (Q_int u k0 k1) (Q_int u k0 k1)))
             L).
    + apply (@is_derive_mult R_AbsRing
               (fun u : R => u) (fun u : R => u) L).
      * apply (@is_derive_id R_AbsRing L).
      * apply (@is_derive_id R_AbsRing L).
      * apply Rmult_comm.
    + apply (@is_derive_plus R_AbsRing R_NormedModule
               (fun u : R => mult (P_int u k0 k1) (P_int u k0 k1))
               (fun u : R => mult (Q_int u k0 k1) (Q_int u k0 k1))
               L).
      * apply (@is_derive_mult R_AbsRing
                 (fun u : R => P_int u k0 k1)
                 (fun u : R => P_int u k0 k1) L).
        apply is_derive_P_at. apply is_derive_P_at. apply Rmult_comm.
      * apply (@is_derive_mult R_AbsRing
                 (fun u : R => Q_int u k0 k1)
                 (fun u : R => Q_int u k0 k1) L).
        apply is_derive_Q_at. apply is_derive_Q_at. apply Rmult_comm.
    + apply Rmult_comm.
  - apply (@is_derive_const R_AbsRing R_NormedModule d2 L).
Qed.

(** Differentiability of fL — follows immediately from f_L_prime_eq_raw. *)
Corollary ex_derive_fL : forall L k0 k1 d2, ex_derive (fun u => fL u k0 k1 d2) L.
Proof. intros. eexists. apply f_L_prime_eq_raw. Qed.

(** Corollary form using Derive.  The raw expression for Derive fL is what
    f_L_prime_eq_raw gives; this packages it. *)
Corollary Derive_fL_eq :
  forall L k0 k1 d2,
    Derive (fun u => fL u k0 k1 d2) L
    = minus
        (plus (mult (plus (mult 1 L) (mult L 1))
                    (plus (mult (P_int L k0 k1) (P_int L k0 k1))
                          (mult (Q_int L k0 k1) (Q_int L k0 k1))))
              (mult (mult L L)
                    (plus
                       (plus (mult (- T_int L k0 k1) (P_int L k0 k1))
                             (mult (P_int L k0 k1) (- T_int L k0 k1)))
                       (plus (mult (R_int L k0 k1) (Q_int L k0 k1))
                             (mult (Q_int L k0 k1) (R_int L k0 k1))))))
        0.
Proof. intros. apply is_derive_unique. apply f_L_prime_eq_raw. Qed.

(** Algebraic identity (paper form): the value above equals
    2*L*(P^2+Q^2) + 2*L^2*(Q*R - P*T).
    Proved in plain R, no ring-on-AbsRing needed. *)
Lemma fL_prime_value_simpl :
  forall (L : R) (P Q R T : R),
    minus
      (plus (mult (plus (mult 1 L) (mult L 1))
                  (plus (mult P P) (mult Q Q)))
            (mult (mult L L)
                  (plus
                     (plus (mult (- T) P) (mult P (- T)))
                     (plus (mult R Q) (mult Q R)))))
      0
    = 2 * L * (P * P + Q * Q) + 2 * L * L * (Q * R - P * T).
Proof.
  intros. unfold minus, plus, mult, opp; simpl.
  unfold Rminus; ring.
Qed.

(** Clean simplified form, obtained by transporting Derive_fL_eq through
    fL_prime_value_simpl. *)
Corollary Derive_fL :
  forall L k0 k1 d2,
    Derive (fun u => fL u k0 k1 d2) L
    = 2 * L * (P_int L k0 k1 * P_int L k0 k1
               + Q_int L k0 k1 * Q_int L k0 k1)
      + 2 * L * L * (Q_int L k0 k1 * R_int L k0 k1
                     - P_int L k0 k1 * T_int L k0 k1).
Proof.
  intros. rewrite Derive_fL_eq.
  apply fL_prime_value_simpl.
Qed.

(** Clean simplified form: f'(L) = 2 L (P^2+Q^2) + 2 L^2 (Q R - P T). *)
Theorem f_L_prime_eq :
  forall (L k0 k1 d2 : R),
    is_derive (fun u => fL u k0 k1 d2) L
              (2 * L * (P_int L k0 k1 * P_int L k0 k1
                        + Q_int L k0 k1 * Q_int L k0 k1)
               + 2 * L * L * (Q_int L k0 k1 * R_int L k0 k1
                              - P_int L k0 k1 * T_int L k0 k1)).
Proof.
  intros L k0 k1 d2.
  rewrite <- (Derive_fL L k0 k1 d2).
  apply Derive_correct. apply ex_derive_fL.
Qed.

(* ==================================================================== *)
(* 7. f''(L): the second derivative of the chord-length residual         *)
(* ==================================================================== *)

(** The f'(L) value, as a function of L. *)
Definition fpL (L k0 k1 : R) : R :=
  2 * L * (P_int L k0 k1 * P_int L k0 k1 + Q_int L k0 k1 * Q_int L k0 k1)
  + 2 * L * L * (Q_int L k0 k1 * R_int L k0 k1
                 - P_int L k0 k1 * T_int L k0 k1).

(** The second derivative of the chord-length residual.

    Proved by [auto_derive] (which produces nine [ex_derive] side
    conditions and a goal containing [Derive ...] placeholders for the
    four moment integrals), then rewriting those placeholders via the
    integral-level derivative lemmas [Derive_P_int], [Derive_Q_int],
    [is_derive_R], [is_derive_T], and closing by [ring]. *)
Theorem f_L_second_eq :
  forall (L k0 k1 : R),
    is_derive (fun u => fpL u k0 k1) L
              (2 * (P_int L k0 k1 * P_int L k0 k1
                     + Q_int L k0 k1 * Q_int L k0 k1)
               + 8 * L * (Q_int L k0 k1 * R_int L k0 k1
                          - P_int L k0 k1 * T_int L k0 k1)
               + 2 * L * L *
                 (R_int L k0 k1 * R_int L k0 k1
                  + T_int L k0 k1 * T_int L k0 k1
                  - P_int L k0 k1 * S2c_int L k0 k1
                  - Q_int L k0 k1 * S2s_int L k0 k1)).
Proof.
  intros L k0 k1.
  unfold fpL.
  auto_derive.
  - split. apply ex_derive_P_int.
    split. apply ex_derive_P_int.
    split. apply ex_derive_Q_int.
    split. apply ex_derive_Q_int.
    split. apply ex_derive_Q_int.
    split. eexists; apply is_derive_R.
    split. apply ex_derive_P_int.
    split. eexists; apply is_derive_T.
    exact I.
  - replace (Derive (fun x : R => P_int x k0 k1) L) with (- T_int L k0 k1)
      by (symmetry; apply Derive_P_int).
    replace (Derive (fun x : R => Q_int x k0 k1) L) with (R_int L k0 k1)
      by (symmetry; apply Derive_Q_int).
    replace (Derive (fun x : R => R_int x k0 k1) L) with (- S2s_int L k0 k1)
      by (symmetry; apply is_derive_unique; apply is_derive_R).
    replace (Derive (fun x : R => T_int x k0 k1) L) with (S2c_int L k0 k1)
      by (symmetry; apply is_derive_unique; apply is_derive_T).
    ring.
Qed.
