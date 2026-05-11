(** Moment-notation polishing for the Bertolazzi-Frego residual.

    Adds the named definitions f, X, Y to match the paper's notation, and
    proves the moment-form versions of the derivative identities:

        f'(A)  = X_2(2A, delta-A, phi0) - X_1(2A, delta-A, phi0)
        f''(A) = -(Y_4 - 2 Y_3 + Y_2)  (all at (2A, delta-A, phi0))

    These follow directly from f_prime_eq and f_second_eq in Clothoid.v by
    rewriting the integrand using the algebraic identities of section 1. *)

Require Import Reals.
Require Import Coquelicot.Coquelicot.
Require Import Clothoid.
Open Scope R_scope.

(** Helpers: bridge Coq Reals' [-] / Coquelicot's [opp] for RInt. *)
Lemma RInt_Ropp_R :
  forall (g : R -> R) (a b : R),
    ex_RInt g a b ->
    RInt (fun x => - g x) a b = - RInt g a b.
Proof.
  intros g a b H.
  change (- RInt g a b) with (opp (RInt g a b)).
  rewrite <- (RInt_opp g a b H). apply RInt_ext. intros. reflexivity.
Qed.

Lemma RInt_Rminus_R :
  forall (fr gr : R -> R) (a b : R),
    ex_RInt fr a b -> ex_RInt gr a b ->
    RInt (fun x => fr x - gr x) a b = RInt fr a b - RInt gr a b.
Proof.
  intros fr gr a b Hf Hg.
  change (RInt fr a b - RInt gr a b) with (minus (RInt fr a b) (RInt gr a b)).
  rewrite <- RInt_minus by assumption. apply RInt_ext. intros. reflexivity.
Qed.

(** Bertolazzi-Frego scalar residual. *)
Definition f (A delta phi0 : R) : R :=
  RInt (fun t => sin (phi A t delta phi0)) 0 1.

(** Generalized Fresnel moments, exactly as in the paper:
      X_k(a, b, c) = int_0^1 t^k * cos((a/2)*t^2 + b*t + c) dt
      Y_k(a, b, c) = int_0^1 t^k * sin((a/2)*t^2 + b*t + c) dt
    With (a, b, c) = (2A, delta-A, phi0), the phase reduces to
      (a/2)*t^2 + b*t + c
       = A*t^2 + (delta-A)*t + phi0
       = phi A t delta phi0 . *)
Definition X (k : nat) (a b c : R) : R :=
  RInt (fun t => t^k * cos ((a / 2) * t^2 + b * t + c)) 0 1.

Definition Y (k : nat) (a b c : R) : R :=
  RInt (fun t => t^k * sin ((a / 2) * t^2 + b * t + c)) 0 1.

(** Phase identity: with (a, b, c) = (2A, delta-A, phi0), the X/Y
    integrand reduces to t^k * trig(phi A t delta phi0). *)
Lemma phase_via_X :
  forall (A t delta phi0 : R),
    (2 * A / 2) * t^2 + (delta - A) * t + phi0 = phi A t delta phi0.
Proof. intros. unfold phi. field. Qed.

(** Algebraic moment identity needed for f': X_2 - X_1
    = int_0^1 (t^2 - t) * cos(phi). *)
Lemma X2_minus_X1_eq :
  forall (A delta phi0 : R),
    X 2 (2 * A) (delta - A) phi0 - X 1 (2 * A) (delta - A) phi0
    = RInt (fun t => (t * t - t) * cos (phi A t delta phi0)) 0 1.
Proof.
  intros A delta phi0. unfold X.
  rewrite <- RInt_Rminus_R.
  - apply RInt_ext. intros t _.
    rewrite phase_via_X. simpl. unfold Rminus. ring.
  - apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
    intros t _.
    apply (@ex_derive_continuous R_AbsRing R_NormedModule).
    eexists. auto_derive; trivial.
  - apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
    intros t _.
    apply (@ex_derive_continuous R_AbsRing R_NormedModule).
    eexists. auto_derive; trivial.
Qed.

(** Algebraic moment identity needed for f'':
    -(Y_4 - 2*Y_3 + Y_2) = int_0^1 -(t^2 - t)^2 * sin(phi). *)

(** Helper: an R-flavoured scalar-pull for RInt. *)
Lemma RInt_Rscal_R :
  forall (gr : R -> R) (k a b : R),
    ex_RInt gr a b ->
    k * RInt gr a b = RInt (fun x => k * gr x) a b.
Proof.
  intros gr k a b H.
  change (k * RInt gr a b) with (scal k (RInt gr a b)).
  rewrite <- (RInt_scal gr a b k H). apply RInt_ext. intros. reflexivity.
Qed.

(** Helper: addition for R-valued RInts. *)
Lemma RInt_Rplus_R :
  forall (fr gr : R -> R) (a b : R),
    ex_RInt fr a b -> ex_RInt gr a b ->
    RInt fr a b + RInt gr a b = RInt (fun x => fr x + gr x) a b.
Proof.
  intros fr gr a b Hf Hg.
  change (RInt fr a b + RInt gr a b) with (plus (RInt fr a b) (RInt gr a b)).
  rewrite <- RInt_plus by assumption. apply RInt_ext. intros. reflexivity.
Qed.

Lemma Y_pattern_eq :
  forall (A delta phi0 : R),
    - (Y 4 (2 * A) (delta - A) phi0
       - 2 * Y 3 (2 * A) (delta - A) phi0
       + Y 2 (2 * A) (delta - A) phi0)
    = RInt (fun t =>
              - ((t * t - t) * (t * t - t))
              * sin (phi A t delta phi0)) 0 1.
Proof.
  intros A delta phi0.
  (* Step 1: rewrite each Y in terms of phi (the phase identity). *)
  assert (HY4 : Y 4 (2 * A) (delta - A) phi0
                = RInt (fun t => t^4 * sin (phi A t delta phi0)) 0 1).
  { unfold Y. apply RInt_ext. intros t _.
    rewrite phase_via_X. reflexivity. }
  assert (HY3 : Y 3 (2 * A) (delta - A) phi0
                = RInt (fun t => t^3 * sin (phi A t delta phi0)) 0 1).
  { unfold Y. apply RInt_ext. intros t _.
    rewrite phase_via_X. reflexivity. }
  assert (HY2 : Y 2 (2 * A) (delta - A) phi0
                = RInt (fun t => t^2 * sin (phi A t delta phi0)) 0 1).
  { unfold Y. apply RInt_ext. intros t _.
    rewrite phase_via_X. reflexivity. }
  rewrite HY4, HY3, HY2.
  (* Step 2: integrability witnesses we'll need. *)
  assert (ex4 : ex_RInt (fun t => t^4 * sin (phi A t delta phi0)) 0 1).
  { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
    intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
    unfold phi. eexists. auto_derive; trivial. }
  assert (ex3 : ex_RInt (fun t => t^3 * sin (phi A t delta phi0)) 0 1).
  { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
    intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
    unfold phi. eexists. auto_derive; trivial. }
  assert (ex2 : ex_RInt (fun t => t^2 * sin (phi A t delta phi0)) 0 1).
  { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
    intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
    unfold phi. eexists. auto_derive; trivial. }
  assert (ex_target :
    ex_RInt (fun t => -((t * t - t) * (t * t - t))
                      * sin (phi A t delta phi0)) 0 1).
  { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
    intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
    unfold phi. eexists. auto_derive; trivial. }
  (* Step 3: combine the three Y-integrals into one big integral via
     linearity (RInt_Rminus_R, RInt_Rplus_R, RInt_Rscal_R, RInt_Ropp_R). *)
  (* Compute 2 * RInt(t^3 sin) as RInt(2 * t^3 sin) *)
  rewrite (RInt_Rscal_R _ 2 _ _ ex3).
  (* Now: - (RInt(t^4 sin) - RInt(2 * t^3 sin) + RInt(t^2 sin))
     Combine: RInt(t^4 sin) - RInt(2 * t^3 sin) = RInt(t^4 sin - 2 t^3 sin) *)
  rewrite <- (RInt_Rminus_R
                (fun t => t^4 * sin (phi A t delta phi0))
                (fun t => 2 * (t^3 * sin (phi A t delta phi0)))
                0 1 ex4).
  2: { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
       intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
       unfold phi. eexists. auto_derive; trivial. }
  (* Then add RInt(t^2 sin): RInt(...) + RInt(t^2 sin) = RInt(... + t^2 sin) *)
  rewrite RInt_Rplus_R.
  2: { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
       intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
       unfold phi. eexists. auto_derive; trivial. }
  2: { exact ex2. }
  (* Pull the leading minus into the integrand. *)
  rewrite <- RInt_Ropp_R.
  2: { apply (ex_RInt_continuous (V := R_CompleteNormedModule)).
       intros t _. apply (@ex_derive_continuous R_AbsRing R_NormedModule).
       unfold phi. eexists. auto_derive; trivial. }
  (* The two integrals are now equal by RInt_ext + ring. *)
  apply RInt_ext. intros t _. simpl. ring.
Qed.

(** Polishing lemma 1: f'(A) in moment notation. *)
Theorem f_prime_moment_notation :
  forall A delta phi0,
    is_derive (fun a => f a delta phi0) A
              (X 2 (2 * A) (delta - A) phi0
               - X 1 (2 * A) (delta - A) phi0).
Proof.
  intros A delta phi0.
  pose proof (f_prime_eq A delta phi0) as H.
  unfold f. rewrite X2_minus_X1_eq. exact H.
Qed.

(** Polishing lemma 2: f''(A) in moment notation.

    Derives directly from f_second_eq + Y_pattern_eq.

    [f_second_eq]   says  d/dA f'(A) = RInt (fun t => -(t^2-t)^2 * sin(phi)) 0 1
    [Y_pattern_eq]  says  RInt (...) = -(Y_4 - 2 Y_3 + Y_2). *)
Theorem f_second_moment_notation :
  forall A delta phi0,
    is_derive
      (fun a => RInt (fun t => (t * t - t) * cos (phi a t delta phi0)) 0 1)
      A
      (- (Y 4 (2 * A) (delta - A) phi0
          - 2 * Y 3 (2 * A) (delta - A) phi0
          + Y 2 (2 * A) (delta - A) phi0)).
Proof.
  intros A delta phi0.
  rewrite Y_pattern_eq.
  apply f_second_eq.
Qed.
