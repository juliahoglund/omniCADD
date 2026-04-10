#!/bin/bash
# Debug script to test CI environment behavior
# This helps identify bash version differences and arithmetic behavior

echo "=== Environment Debug ==="
echo "Bash version: $BASH_VERSION"
echo "Shell: $SHELL"
type bash
echo ""

echo "=== Testing ((PASSED++)) with set -e ===="
set -e
set -o pipefail

PASSED=0
FAILED=0

test_pass() {
    echo "  test_pass() called"
    return 0
}

test_fail() {
    echo "  test_fail() called"
    return 1
}

echo "Starting: PASSED=$PASSED FAILED=$FAILED"
echo ""

echo "Test 1: Should increment PASSED from 0 to 1"
if test_pass; then
    echo "  Before: PASSED=$PASSED"
    ((PASSED++)) && echo "  Arithmetic succeeded" || echo "  Arithmetic returned false but script continued"
    echo "  After: PASSED=$PASSED"
else
    ((FAILED++))
    echo "  FAILED=$FAILED"
fi
echo ""

echo "Test 2: Should increment PASSED from 1 to 2"
if test_pass; then
    echo "  Before: PASSED=$PASSED"
    ((PASSED++)) && echo "  Arithmetic succeeded" || echo "  Arithmetic returned false but script continued"
    echo "  After: PASSED=$PASSED"
else
    ((FAILED++))
    echo "  FAILED=$FAILED"
fi
echo ""

echo "Test 3: Should increment FAILED from 0 to 1"
if test_fail; then
    ((PASSED++))
else
    echo "  Before: FAILED=$FAILED"
    ((FAILED++)) && echo "  Arithmetic succeeded" || echo "  Arithmetic returned false but script continued"
    echo "  After: FAILED=$FAILED"
fi
echo ""

echo "=== Final Results ==="
echo "PASSED=$PASSED (expected: 2)"
echo "FAILED=$FAILED (expected: 1)"
echo ""

if [ "$PASSED" -eq 2 ] && [ "$FAILED" -eq 1 ]; then
    echo "✓ All arithmetic operations worked correctly"
    exit 0
else
    echo "✗ Unexpected values!"
    exit 1
fi
