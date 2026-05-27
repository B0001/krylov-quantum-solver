#!/usr/bin/env bash
# ==============================================================================
# WORKSPACE SANITIZATION PIPELINE: clean_workspace.sh
# Purges residual caching blocks, build artifacts, and compiled bytecodes 
# to return the repository to a pristine, deployment-ready state.
# ==============================================================================

echo "================================================================================"
echo "INITIATING REPOSITORY SANITIZATION PASS"
echo "================================================================================"

# 1. Purge Python Bytecode & PyTest Cache Directories
echo "[STEP 1] Locating and destroying __pycache__ blocks and compiled bytecode (.pyc)..."
find . -type d -name "__pycache__" -exec rm -rf {} +
find . -type d -name ".pytest_cache" -exec rm -rf {} +
find . -type f -name "*.pyc" -delete
find . -type f -name "*.pyo" -delete

# 2. Erase Build & Distribution Artifacts
echo "[STEP 2] Removing package distributions and egg-info metadata bounds..."
rm -rf dist/
rm -rf build/
rm -rf *.egg-info/

# 3. Clean Sphinx Documentation Builds
echo "[STEP 3] Wiping generated documentation HTML/Markdown builds..."
rm -rf docs/build/

# 4. Remove Editor Backup & OS Temporary Files
echo "[STEP 4] Scrubbing trailing editor backup artifacts (*~) and OS tracking files..."
find . -type f -name "*~" -delete
find . -type f -name ".DS_Store" -delete

echo "--------------------------------------------------------------------------------"
echo "[SUCCESS] Workspace is now clean, lean, and strictly source-controlled."
echo "================================================================================"