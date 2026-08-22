# Development targets for SCEFitting.jl. Every test target pins
# JULIA_NUM_THREADS=4: the suite REFUSES to run at one thread (the threaded-vs-
# serial gates are vacuous there), and CI pins the same value.

JULIA ?= julia
THREADS ?= 4
export JULIA_NUM_THREADS := $(THREADS)

.PHONY: test-all test-unit test-aqua test-jet test-oracle test-sunny test-glmnet \
        test-pin test-parity test-examples test-downstream test-ci ci-local docs \
        bench-setup bench-salcbasis bench-clusters bench-design-matrix bench-solver \
        bench-end-to-end bench-nd2fe14b

# ---- core suite (runtests.jl dispatches on TEST_MODE) -----------------------

test-all:
	TEST_MODE=all $(JULIA) --project -e 'using Pkg; Pkg.test()'

test-unit:
	TEST_MODE=unit $(JULIA) --project -e 'using Pkg; Pkg.test()'

test-aqua:
	TEST_MODE=aqua $(JULIA) --project -e 'using Pkg; Pkg.test()'

test-jet:
	TEST_MODE=jet $(JULIA) --project -e 'using Pkg; Pkg.test()'

# ---- separate environments (each carries the heavy / optional dependency) ---

# Pinned Magesty.jl as a numerical oracle. Needs a sibling Magesty.jl checkout;
# local-only (not in CI).
test-oracle:
	$(JULIA) --project=test/oracle -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=test/oracle test/oracle/runtests.jl

test-sunny:
	$(JULIA) --project=test/sunny -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=test/sunny test/sunny/runtests.jl

test-glmnet:
	$(JULIA) --project=test/glmnet -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=test/glmnet test/glmnet/runtests.jl

# Byte-level change detectors over the SALC chain (NOT correctness evidence —
# see test/pin/PIN.md). Both thread counts, every time: a pinned value that
# depends on the thread count is a race, not a pin.
test-pin:
	$(JULIA) --project=test/pin -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=test/pin -t 4 test/pin/runtests.jl
	$(JULIA) --project=test/pin -t 1 test/pin/runtests.jl

# Real-data parity of the moment channel vs the sibling SLCE.jl checkout.
# Needs SCE_PARITY_FEGE_DIR / SCE_PARITY_FEGE_POSCAR / SCE_PARITY_FERH_DIR;
# skips loudly without them. No CI job.
test-parity:
	$(JULIA) --project=test/parity -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=test/parity -t 4 test/parity/runtests.jl

# Every example self-gates with @asserts on a recovered coupling, so they have
# to RUN somewhere (CI's `runnable examples` job).
test-examples:
	$(JULIA) --project=examples -e 'using Pkg; Pkg.instantiate()'
	for f in examples/*.jl; do echo "== $$f"; $(JULIA) --project=examples "$$f" || exit 1; done

# The dependent's suite against THIS checkout (CI's `downstream` job). Needs a
# sibling SCEMonteCarlo.jl checkout whose Manifest develops ../SCEFitting.jl.
test-downstream:
	$(JULIA) --project=../SCEMonteCarlo.jl -e 'using Pkg; Pkg.test()'

# Strict Documenter build (warnonly = false, checkdocs = :public). Executes
# every @example block, so it is a test as much as a build.
docs:
	$(JULIA) --project=docs -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=docs docs/make.jl

# ---- CI parity ---------------------------------------------------------------

# Exactly the jobs GitHub Actions runs (minus `downstream`, which needs the
# sibling checkout — run test-downstream separately — and the oracle, which CI
# cannot reach). Run before a release or version bump.
test-ci: test-all test-sunny test-glmnet test-examples test-pin docs

# Cold-start reproduction of CI: no cached Manifest, juliaup `release` channel
# (matches CI's "1.12" while release resolves to 1.12.x). `juliaup add release`
# if missing.
ci-local:
	rm -f Manifest.toml
	TEST_MODE=all $(JULIA) +release --project -e 'using Pkg; Pkg.test()'
	$(MAKE) test-sunny test-glmnet test-examples test-pin docs

# ---- benchmarks (bench/README.md documents fixtures and arguments) ----------

bench-setup:
	$(JULIA) --project=bench -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'

bench-salcbasis:
	$(JULIA) --project=bench bench/bench_salcbasis.jl

bench-clusters:
	$(JULIA) --project=bench bench/bench_clusters.jl

bench-design-matrix:
	$(JULIA) --project=bench bench/bench_design_matrix.jl

bench-solver:
	$(JULIA) --project=bench bench/bench_solver.jl

bench-end-to-end:
	$(JULIA) --project=bench bench/bench_end_to_end.jl

bench-nd2fe14b:
	$(JULIA) --project=bench bench/bench_nd2fe14b.jl
