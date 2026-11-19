# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Steven Van Overberghe
#
# genll -- exhaustive generation of locally linear graphs
#
#   make NAUTY=/path/to/nauty          build genll (N <= 32) and genllL (N <= 64)
#   make check NAUTY=/path/to/nauty    build, then compare against known counts
#   make clean
#
# NAUTY is a built nauty source directory: it must contain nauty.h and the archives
# nautyW1.a and nautyL1.a, which nauty's own "./configure && make" produces.
#
# Optional compile-time settings, e.g. "make STATS=1":
#   STATS=1     compile in per-candidate statistics (the >S and >K lines on stderr);
#               without it the counters are empty and cost nothing
#   ROUNDS=n    rounds of colour refinement in the last rung of the vertex
#               invariant (default 3)
#   VERIFY=1    check every incrementally computed invariant against its direct
#               definition (slow; for testing only)

NAUTY    ?= ../nauty
CXX      ?= g++
OPT      ?= -O3 -march=native -flto=auto
WARN     ?= -Wall -Wno-sign-compare

FLAGS    = -std=c++20 $(OPT) -DNDEBUG $(WARN) -I$(NAUTY) $(CXXFLAGS)
ifdef STATS
FLAGS   += -DGENLL_STATS=1
endif
ifdef ROUNDS
FLAGS   += -DREFINEMENT_ROUNDS=$(ROUNDS)
endif
ifdef VERIFY
FLAGS   += -DGENLL_VERIFY_KAPPA
endif

HEADERS  = nautyBitset.h graph.h

all: genll genllL

# 32-bit set words: N <= 32.  MAXN=WORDSIZE matches how nautyW1.a was compiled.
genll: genll.cpp $(HEADERS) $(NAUTY)/nautyW1.a
	$(CXX) $(FLAGS) -DMAXN=WORDSIZE -DWORDSIZE=32 -o $@ genll.cpp $(NAUTY)/nautyW1.a $(LDFLAGS)

# 64-bit set words: N <= 64.  MAXN=WORDSIZE matches how nautyL1.a was compiled.
genllL: genll.cpp $(HEADERS) $(NAUTY)/nautyL1.a
	$(CXX) $(FLAGS) -DMAXN=WORDSIZE -DWORDSIZE=64 -o $@ genll.cpp $(NAUTY)/nautyL1.a $(LDFLAGS)

$(NAUTY)/nautyW1.a $(NAUTY)/nautyL1.a:
	@echo "$@ not found." >&2
	@echo "Build nauty first (./configure && make in its directory) and point" >&2
	@echo "NAUTY at it, e.g.  make NAUTY=/path/to/nauty" >&2
	@exit 1

check: genll genllL
	./check.sh

clean:
	rm -f genll genllL

.PHONY: all check clean
