#!/usr/bin/env bash
# =============================================================================
# test_build.sh — Reproducible Cantera build & test script for Fedora + GCC 16
# =============================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

PASS=0
FAIL=0
step() {
    local n="$1"; shift
    echo ""
    echo "========================================"
    echo "  STEP $n: $*"
    echo "========================================"
}
pass() { echo "✅ STEP $1 PASSED"; ((PASS++)); }
fail() { echo "❌ STEP $1 FAILED"; ((FAIL++)); }

# ── Step 1: Clean repository ────────────────────────────────────────────────
step 1 "Clean repository"
git clean -fdx -e test_build.sh
git submodule deinit -f --all
pass 1

# ── Step 2: Update submodules ───────────────────────────────────────────────
step 2 "Update submodules"
git submodule update --init --recursive
pass 2

# ── Step 3: Install system dependencies ─────────────────────────────────────
step 3 "Install system dependencies (requires sudo)"
sudo dnf install -y \
    gcc-c++ python3 scons boost-devel hdf5-devel \
    python3-setuptools python3-jinja2 python3-devel python3-wheel
pass 3

# ── Step 4: Install Python dependencies via pip ─────────────────────────────
step 4 "Install Python dependencies via pip"
pip install numpy wheel cython ruamel.yaml pytest setuptools jinja2
pass 4

# ── Step 5: Apply GCC 16 patch (add #include <cstdint>) ────────────────────
step 5 "Apply GCC 16 patch to ext/yaml-cpp/src/emitterutils.cpp"
TARGET_FILE="ext/yaml-cpp/src/emitterutils.cpp"
if [ ! -f "$TARGET_FILE" ]; then
    echo "ERROR: $TARGET_FILE not found — submodule init may have failed."
    fail 5
else
    if grep -q '#include <cstdint>' "$TARGET_FILE"; then
        echo "Patch already applied — skipping."
    else
        # Insert #include <cstdint> after the last existing #include block
        # near the top of the file (around line 12).
        sed -i '/#include <sstream>/a #include <cstdint>' "$TARGET_FILE"
        echo "Injected #include <cstdint> into $TARGET_FILE"
    fi
    # Verify the patch is present
    if grep -q '#include <cstdint>' "$TARGET_FILE"; then
        pass 5
    else
        echo "ERROR: Patch verification failed."
        fail 5
    fi
fi

# ── Step 6: Compile core (disable doxygen to avoid Perl errors) ────────────
step 6 "Compile Cantera core (doxygen_docs=n)"
scons build doxygen_docs=n -j4
pass 6

# ── Step 7: Test core ──────────────────────────────────────────────────────
step 7 "Run Cantera test suite"
scons test
pass 7

# ── Step 8: Install Python module from wheel ────────────────────────────────
step 8 "Install Python module from built wheel"
pip install build/python/dist/*.whl
pass 8

# ── Step 9: Verify Cantera import ───────────────────────────────────────────
step 9 "Verify Cantera import and basic usage"
python3 -c "
import cantera as ct
gas = ct.Solution('gri30.yaml')
gas()
print('Cantera version:', ct.__version__)
print('GRI-Mech 3.0 loaded successfully with', gas.n_species, 'species')
"
pass 9

# ── Summary ─────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  BUILD SUMMARY"
echo "========================================"
TOTAL=$((PASS + FAIL))
echo "  Passed: $PASS / $TOTAL"
echo "  Failed: $FAIL / $TOTAL"
if [ "$FAIL" -eq 0 ]; then
    echo "  🎉 ALL STEPS PASSED — 100% SUCCESS"
    exit 0
else
    echo "  ⚠️  Some steps failed. Review output above."
    exit 1
fi
