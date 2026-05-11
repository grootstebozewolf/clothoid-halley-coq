# Build the Coq proof file.
# Requires: Coq + Coquelicot in PATH (tested with Coq 8.13.1 + Coquelicot 3.x).

COQC ?= coqc
COQCHK ?= coqchk

all: Clothoid.vo

Clothoid.vo: Clothoid.v
	$(COQC) Clothoid.v

check: Clothoid.vo
	$(COQCHK) -silent -o Clothoid

audit: Clothoid.vo
	@echo "----- axioms required by f_prime_eq -----"
	@echo 'Require Import Clothoid. Print Assumptions f_prime_eq.' | $(COQC) -R . "" /dev/stdin
	@echo "----- axioms required by f_second_eq -----"
	@echo 'Require Import Clothoid. Print Assumptions f_second_eq.' | $(COQC) -R . "" /dev/stdin

clean:
	rm -f Clothoid.vo Clothoid.glob Clothoid.vos Clothoid.vok .Clothoid.aux

.PHONY: all check audit clean
