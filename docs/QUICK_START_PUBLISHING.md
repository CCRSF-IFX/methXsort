# Quick Start: Publishing methXsort

## TL;DR - What You Need to Do

### 1️⃣ One-Time Setup (15 minutes)

```bash
# a. Create accounts (if you don't have them)
# - Test PyPI: https://test.pypi.org/account/register/
# - Production PyPI: https://pypi.org/account/register/

# b. Generate API tokens
# - Test PyPI: Account Settings → API tokens → Add token (scope: entire account)
# - Production PyPI: Account Settings → API tokens → Add token (scope: entire account)

# c. Add secrets to GitHub repository
# - Go to: https://github.com/CCRSF-IFX/methXsort/settings/secrets/actions
# - Click "New repository secret"
# - Add: TEST_PYPI_API_TOKEN (paste your Test PyPI token)
# - Add: PYPI_API_TOKEN (paste your Production PyPI token)
```

### 2️⃣ Test Your First Release (5 minutes)

```bash
cd /mnt/ccrsf-ifx/Software/github/methXsort

# Test build locally
python -m pip install --upgrade build twine
python -m build
twine check dist/*

# Push to Test PyPI via GitHub Actions
git tag v0.2.0-test
git push origin v0.2.0-test

# Wait 2-3 minutes, then test installation in a new environment
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ methXsort
methXsort --help
methXsort --version
```

### 3️⃣ Production Release (5 minutes)

```bash
# Go to GitHub → Releases → Draft new release
# URL: https://github.com/CCRSF-IFX/methXsort/releases/new
# - Tag: v0.2.0
# - Title: methXsort v0.2.0
# - Description: [Your release notes - see example below]
# - Publish release

# Wait 2-3 minutes, then test installation
pip install methXsort
methXsort --version
```

**Example Release Notes:**
```markdown
## methXsort v0.2.0

### 🚀 New Features
- Converted to proper Python package structure
- Added `methxsort` command-line tool
- Now installable via pip

### ⚙️ Changes
- Requires Python 3.12+
- Removed toolshed dependency (replaced with built-in functions)
- Updated to use pyproject.toml (modern packaging)

### 📦 Installation
```bash
pip install methXsort
```

### 🔧 Usage
```bash
methXsort --help
methXsort convert-ref reference.fa
```

See [README.md](https://github.com/CCRSF-IFX/methXsort#readme) for full documentation.
```

## That's It! 🎉

Your package is now on PyPI and anyone can install it with:
```bash
pip install methXsort
```

---

## Files You'll Need to Create

To enable PyPI publishing, you need to create GitHub Actions workflows:

### 1. **`.github/workflows/publish-test-pypi.yml`**
   - Triggers on: Tags matching `v*-test` (e.g., `v0.2.0-test`)
   - Publishes to: Test PyPI
   - Use for: Testing releases before production

### 2. **`.github/workflows/publish-pypi.yml`**
   - Triggers on: GitHub Releases
   - Publishes to: Production PyPI
   - Use for: Official releases

### 3. **`.github/workflows/tests.yml`** ✅ (Already exists)
   - Triggers on: Push to main/develop, Pull Requests
   - Tests: Package installation and CLI commands
   - Use for: Continuous Integration

---

## Quick Workflow Setup

Create these two files:

**`.github/workflows/publish-test-pypi.yml`:**
```yaml
name: Publish to Test PyPI

on:
  push:
    tags:
      - 'v*-test'  # Triggers on tags like v0.2.0-test

jobs:
  publish:
    name: Publish to Test PyPI
    runs-on: ubuntu-latest
    permissions:
      contents: read
      
    steps:
    - name: Checkout code
      uses: actions/checkout@v4
      
    - name: Set up Python
      uses: actions/setup-python@v5
      with:
        python-version: '3.12'
        
    - name: Install build dependencies
      run: |
        python -m pip install --upgrade pip
        pip install build twine
        
    - name: Build package
      run: python -m build
      
    - name: Check package
      run: twine check dist/*
      
    - name: Publish to Test PyPI
      env:
        TWINE_USERNAME: __token__
        TWINE_PASSWORD: ${{ secrets.TEST_PYPI_API_TOKEN }}
      run: |
        twine upload --repository testpypi dist/*
```

**`.github/workflows/publish-pypi.yml`:**
```yaml
name: Publish to PyPI

on:
  release:
    types: [published]  # Triggers when you publish a release on GitHub

jobs:
  publish:
    name: Publish to PyPI
    runs-on: ubuntu-latest
    permissions:
      contents: read
      
    steps:
    - name: Checkout code
      uses: actions/checkout@v4
      
    - name: Set up Python
      uses: actions/setup-python@v5
      with:
        python-version: '3.12'
        
    - name: Install build dependencies
      run: |
        python -m pip install --upgrade pip
        pip install build twine
        
    - name: Build package
      run: python -m build
      
    - name: Check package
      run: twine check dist/*
      
    - name: Publish to PyPI
      env:
        TWINE_USERNAME: __token__
        TWINE_PASSWORD: ${{ secrets.PYPI_API_TOKEN }}
      run: |
        twine upload dist/*
```

---

## methXsort-Specific Info

### Package Details
- **Package Name**: `methXsort`
- **CLI Command**: `methXsort` (capital X)
- **Current Version**: `0.2.0` (defined in `pyproject.toml` and `methxsort/__init__.py`)
- **Python Requirement**: `>=3.12`
- **Dependencies**: `pysam>=0.15.0`, `xengsort>=2.0.9`

### Subcommands Available
```bash
methXsort convert-ref          # Convert reference genome
methXsort convert-reads        # Convert FASTQ reads
methXsort filter-fastq-by-bam  # Extract reads from FASTQ
methXsort bbsplit              # Split reads with bbsplit
methXsort bbsplit-index        # Build bbsplit index
methXsort stat-split           # Generate statistics
methXsort xengsort-index       # Build xengsort index
methXsort xengsort-classify    # Classify reads
methXsort restore-fastq        # Restore original sequences
methXsort parse-xengsort-summary       # Parse xengsort summary
methXsort parse-xengsort-log-table     # Parse xengsort log
```

---

## Need More Details?

For comprehensive documentation, see:
- **INSTALLATION.md** - Installation instructions and requirements
- **README.md** - Usage examples and workflow diagrams
- **PYPROJECT_TOML_EXPLAINED.md** - Understanding the package configuration
- **PYTHON_312_MIGRATION.md** - Python version requirements
- **GITHUB_ACTIONS_FIX.md** - CI/CD troubleshooting

---

## Common Commands

```bash
# Build package locally
python -m build

# Check package
twine check dist/*

# Clean dist folder (before new build)
rm -rf dist/ build/ *.egg-info

# Test release to Test PyPI
git tag v0.2.0-test
git push origin v0.2.0-test

# Delete a tag (if you made a mistake)
git tag -d v0.2.0-test
git push origin --delete v0.2.0-test

# Production release (use GitHub UI - recommended)
# Go to: https://github.com/CCRSF-IFX/methXsort/releases/new
# Create release with tag v0.2.0

# Or use command line (not recommended, use GitHub UI instead)
git tag v0.2.0
git push origin v0.2.0
# Then create release on GitHub from this tag
```

---

## Pre-Release Checklist

Before creating a release:

- [ ] Update version in `methxsort/__init__.py`
- [ ] Update version in `pyproject.toml`
- [ ] Update CHANGELOG or release notes
- [ ] Test locally: `pip install -e . && methXsort --help`
- [ ] Run tests: GitHub Actions on latest commit
- [ ] Test build: `python -m build && twine check dist/*`
- [ ] Clean up: `rm -rf dist/ build/`

---

## Getting Help

- **Repository**: https://github.com/CCRSF-IFX/methXsort
- **Issues**: https://github.com/CCRSF-IFX/methXsort/issues
- **Email**: xies4@nih.gov
- **Documentation**: See README.md in the repository

---

## Troubleshooting

### "Package already exists" error
- Version already published to PyPI
- Increment version in `pyproject.toml` and `__init__.py`
- PyPI doesn't allow replacing versions (even if deleted)

### "Invalid token" error
- Check GitHub Secrets are correctly named: `TEST_PYPI_API_TOKEN`, `PYPI_API_TOKEN`
- Verify tokens are not expired
- Regenerate tokens if needed

### Build fails
- Check `pyproject.toml` syntax
- Verify all files are committed to git
- Clean build artifacts: `rm -rf dist/ build/ *.egg-info`

### Import fails after installation
- Check package structure in `[tool.setuptools.packages.find]`
- Verify `methxsort/__init__.py` exists
- Check Python version: `python --version` (must be 3.12+)
