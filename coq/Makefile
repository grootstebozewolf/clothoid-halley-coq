# Build the Coq developments.
# Halley residual files: Coq + Coquelicot (tested 8.13.1 / 8.20.1 + Coquelicot 3.x).
# Distance.v / CompoundCurveKoc.v: stdlib Reals only (BSD-3-Clause; see LICENSE_BSD-3-Clause.txt).

COQC ?= coqc
COQCHK ?= coqchk

HALLEY_VO = Clothoid.vo Clothoid_L.vo ClothoidPolish.vo
KOC_VO    = Distance.vo CompoundCurveKoc.vo

all: $(HALLEY_VO) $(KOC_VO)

Clothoid.vo: Clothoid.v
	$(COQC) Clothoid.v

Clothoid_L.vo: Clothoid_L.v
	$(COQC) Clothoid_L.v

ClothoidPolish.vo: ClothoidPolish.v
	$(COQC) ClothoidPolish.v

Distance.vo: Distance.v
	$(COQC) Distance.v

CompoundCurveKoc.vo: CompoundCurveKoc.v Distance.vo
	$(COQC) CompoundCurveKoc.v

check: all
	$(COQCHK) -silent -o Clothoid
	$(COQCHK) -silent -o Clothoid_L
	$(COQCHK) -silent -o ClothoidPolish
	$(COQCHK) -silent -o Distance
	$(COQCHK) -silent -o CompoundCurveKoc

clean:
	rm -f *.vo *.glob *.vos *.vok .*.aux

.PHONY: all check clean
