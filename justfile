# Justfile for Matbench Discovery
# https://github.com/casey/just

set dotenv-load := false

# List available commands
default:
    @just --list

# Prepare a model submission by running all eval and plot scripts and checking the checklist
prepare-model-submission model_name overwrite="false":
    #!/usr/bin/env bash
    set -euo pipefail

    MODEL_INPUT="{{model_name}}"
    OVERWRITE_FLAG=""
    if [ "{{overwrite}}" = "true" ]; then
        OVERWRITE_FLAG="--overwrite"
    fi
    # Normalize: convert hyphens to underscores for enum lookup
    MODEL="${MODEL_INPUT//-/_}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Preparing model submission for: $MODEL"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    # Track checklist status
    CHECKLIST_PASS=0
    CHECKLIST_FAIL=0
    CHECKLIST_SKIP=0

    check_pass() {
        echo "  ✓ $1"
        CHECKLIST_PASS=$((CHECKLIST_PASS + 1))
    }

    check_fail() {
        echo "  ✗ $1"
        CHECKLIST_FAIL=$((CHECKLIST_FAIL + 1))
    }

    check_skip() {
        echo "  ○ $1 (skipped/optional)"
        CHECKLIST_SKIP=$((CHECKLIST_SKIP + 1))
    }

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "STEP 1: Checking PR checklist requirements"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # Get YAML path directly from the Model enum (handles name normalization)
    YAML_FILE=""
    if YAML_FILE=$(uv run python -c "from matbench_discovery.enums import Model; print(Model['$MODEL'].yaml_path)" 2>/dev/null); then
        if [ -f "$YAML_FILE" ]; then
            check_pass "Model YAML file exists: $YAML_FILE"
        else
            check_fail "Model YAML file not found at: $YAML_FILE"
            YAML_FILE=""
        fi
    else
        check_fail "Model '$MODEL' not found in Model enum. Add it to matbench_discovery/enums.py"
    fi

    # Check 2: Model is in the Model enum (already checked above, but confirm)
    if [ -n "$YAML_FILE" ]; then
        check_pass "Model '$MODEL' is registered in the Model enum"
    fi

    # Check model targets (E = energy only, EF/EFS = has forces for relaxation)
    MODEL_TARGETS=""
    if [ -n "$YAML_FILE" ]; then
        MODEL_TARGETS=$(grep -E "^targets:" "$YAML_FILE" 2>/dev/null | sed 's/targets:\s*//' | tr -d ' ')
    fi

    # Check 3-5: Check if prediction file URLs are in YAML
    if [ -n "$YAML_FILE" ]; then
        # Check for discovery pred_file_url
        if grep -q "pred_file_url:" "$YAML_FILE" 2>/dev/null; then
            check_pass "Prediction file URL found in YAML"
        else
            check_fail "No pred_file_url found in $YAML_FILE"
        fi

        # Check for geo_opt pred_file (relaxed structures) - skip for E-only models
        if [ "$MODEL_TARGETS" = "E" ]; then
            check_skip "Relaxed structures check skipped (model targets=$MODEL_TARGETS, no forces)"
        elif grep -A5 "geo_opt:" "$YAML_FILE" 2>/dev/null | grep -q "pred_file:"; then
            check_pass "Relaxed structures file reference found in YAML (geo_opt.pred_file)"
        else
            check_fail "Relaxed structures file (geo_opt.pred_file) not found in YAML"
        fi

        # Check for kappa/phonons predictions - skip for E-only models
        if [ "$MODEL_TARGETS" = "E" ]; then
            check_skip "Phonon predictions check skipped (model targets=$MODEL_TARGETS, no forces)"
        elif grep -A5 "phonons:" "$YAML_FILE" 2>/dev/null | grep -q "pred_file:"; then
            check_pass "Phonon predictions file reference found in YAML"
        else
            check_fail "Phonon predictions (phonons.kappa_103) not found in YAML"
        fi

        # Check for diatomics predictions - skip for E-only models
        if [ "$MODEL_TARGETS" = "E" ]; then
            check_skip "Diatomics predictions check skipped (model targets=$MODEL_TARGETS, no forces)"
        elif grep -A5 "diatomics:" "$YAML_FILE" 2>/dev/null | grep -q "pred_file:"; then
            check_pass "Diatomics predictions file reference found in YAML"
        else
            check_skip "Diatomics predictions not found in YAML"
        fi

        # Check 6: Test scripts exist for each task (discovery, kappa, diatomics)
        MODEL_DIR=$(dirname "$YAML_FILE")
        ARCH_NAME=$(basename "$MODEL_DIR")

        # Check for discovery test script
        if find "$MODEL_DIR" -name "test_*_discovery.py" -o -name "test_${ARCH_NAME}_discovery.py" 2>/dev/null | grep -q .; then
            DISCOVERY_SCRIPT=$(find "$MODEL_DIR" -name "test_*_discovery.py" 2>/dev/null | head -n1)
            check_pass "Discovery test script found: $(basename "$DISCOVERY_SCRIPT")"
        else
            check_fail "Discovery test script not found (test_*_discovery.py) in $MODEL_DIR"
        fi

        # Check for kappa test script - skip for E-only models
        if [ "$MODEL_TARGETS" = "E" ]; then
            check_skip "Kappa test script check skipped (model targets=$MODEL_TARGETS, no forces)"
        elif find "$MODEL_DIR" -name "test_*_kappa.py" -o -name "test_${ARCH_NAME}_kappa.py" 2>/dev/null | grep -q .; then
            KAPPA_SCRIPT=$(find "$MODEL_DIR" -name "test_*_kappa.py" 2>/dev/null | head -n1)
            check_pass "Kappa test script found: $(basename "$KAPPA_SCRIPT")"
        else
            check_fail "Kappa test script not found (test_*_kappa.py) in $MODEL_DIR"
        fi

        # Check for diatomics test script - skip for E-only models
        if [ "$MODEL_TARGETS" = "E" ]; then
            check_skip "Diatomics test script check skipped (model targets=$MODEL_TARGETS, no forces)"
        elif find "$MODEL_DIR" -name "test_*_diatomics.py" -o -name "test_${ARCH_NAME}_diatomics.py" 2>/dev/null | grep -q .; then
            DIATOMICS_SCRIPT=$(find "$MODEL_DIR" -name "test_*_diatomics.py" 2>/dev/null | head -n1)
            check_pass "Diatomics test script found: $(basename "$DIATOMICS_SCRIPT")"
        else
            check_skip "Diatomics test script not found (test_*_diatomics.py) in $MODEL_DIR"
        fi
    fi

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "STEP 2: Running evaluation scripts"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    echo ""
    echo ">> Running discovery metrics evaluation..."
    if uv run python scripts/evals/discovery.py --auto-download --models "$MODEL" $OVERWRITE_FLAG; then
        check_pass "Discovery metrics evaluation completed"
    else
        check_fail "Discovery metrics evaluation failed"
    fi

    echo ""
    echo ">> Running kappa (phonon) metrics evaluation..."
    if [ "$MODEL_TARGETS" = "E" ]; then
        check_skip "Kappa metrics skipped (model targets=$MODEL_TARGETS, no forces)"
    elif uv run python scripts/evals/kappa.py --auto-download --models "$MODEL" $OVERWRITE_FLAG; then
        check_pass "Kappa metrics evaluation completed"
    else
        check_fail "Kappa metrics evaluation failed"
    fi

    echo ""
    echo ">> Running geo_opt analysis..."
    if [ "$MODEL_TARGETS" = "E" ]; then
        check_skip "Geo-opt analysis skipped (model targets=$MODEL_TARGETS, no forces)"
    elif uv run python scripts/analyze_geo_opt.py --auto-download --models "$MODEL" $OVERWRITE_FLAG; then
        check_pass "Geo-opt analysis completed"
    else
        check_fail "Geo-opt analysis failed"
    fi

    echo ""
    echo ">> Running diatomic metrics evaluation..."
    if [ "$MODEL_TARGETS" = "E" ]; then
        check_skip "Diatomic metrics evaluation skipped (model targets=$MODEL_TARGETS, no forces)"
    elif uv run python scripts/evals/diatomic_metrics.py --auto-download --models "$MODEL" $OVERWRITE_FLAG; then
        check_pass "Diatomic metrics evaluation completed"
    else
        check_skip "Diatomic metrics evaluation skipped (no diatomic data or failed)"
    fi


    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "STEP 3: Generating required figures"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    echo ""
    echo ">> Running single model energy parity plot..."
    if uv run python scripts/model_figs/single_model_energy_parity.py --auto-download --models "$MODEL" --update-existing --no-show; then
        check_pass "Energy parity plot generated"
    else
        check_fail "Energy parity plot generation failed"
    fi

    echo ""
    echo ">> Running single model per-element errors plot..."
    if uv run python scripts/model_figs/single_model_per_element_errors.py --auto-download --models "$MODEL" --no-show; then
        check_pass "Per-element errors plot generated"
    else
        check_fail "Per-element errors plot generation failed"
    fi

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "SUMMARY"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  ✓ Passed:  $CHECKLIST_PASS"
    echo "  ✗ Failed:  $CHECKLIST_FAIL"
    echo "  ○ Skipped: $CHECKLIST_SKIP"
    echo ""

    if [ $CHECKLIST_FAIL -gt 0 ]; then
        echo "⚠️  Some checks failed. Please review the PR checklist:"
        echo "   https://github.com/janosh/matbench-discovery/blob/main/.github/pull_request_template.md"
        exit 1
    else
        echo "✅ All required checks passed!"
    fi
