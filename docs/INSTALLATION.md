# methXsort Installation Guide

## Package Structure

methXsort has been converted from a single script into a proper Python package with the following structure:

```
methXsort/
├── methxsort/              # Package directory
│   ├── __init__.py         # Package initialization
│   ├── __main__.py         # CLI entry point
│   └── core.py             # Core functionality
├── setup.py                # Package configuration
├── requirements.txt        # Dependencies
├── README.md              # Documentation
└── methXsort.py           # Original script (preserved for reference)
```

## Installation

### Option 1: Install in Development Mode (Recommended for Development)

From the project root directory:

```bash
pip install -e .
```

This installs the package in "editable" mode, meaning changes to the source code are immediately reflected without reinstalling.

### Option 2: Install from PyPI (When Published)

```bash
pip install methXsort
```

### Option 3: Install from Source

```bash
git clone <repository-url>
cd methXsort
pip install .
```

## Dependencies

- **Python >= 3.12** (Required)
- pysam >= 0.15.0
- xengsort >= 2.0.9

**Note:** The `toolshed` dependency has been removed and replaced with built-in Python functions for better compatibility with modern Python versions.

## Usage

After installation, you can use methXsort as a command:

```bash
# Show version
methXsort --version

# Show help
methXsort --help

# Run subcommands
methXsort convert-ref reference.fa
methXsort convert-reads --read R1.fastq.gz --read2 R2.fastq.gz
```

You can also run it as a Python module:

```bash
python -m methxsort --help
```

## Verification

To verify the installation:

```bash
# Check if methXsort is installed
pip list | grep methXsort

# Test the command
methXsort --version

# Test a subcommand help
methXsort convert-ref --help
```

## Uninstallation

```bash
pip uninstall methXsort
```

## Changes from Original Script

### Key Improvements:
1. **Package Structure**: Organized into a proper Python package
2. **Entry Point**: Can be invoked as `methxsort` command directly
3. **Removed toolshed Dependency**: Replaced with built-in functions for better compatibility
4. **Module Import**: Functions can be imported: `from methxsort.core import convert_fasta`
5. **Version Management**: Centralized version in `__init__.py`

### Backward Compatibility:
- The original `methXsort.py` script is preserved
- All subcommands and functionality remain unchanged
- Command-line interface is identical

## Development

For developers working on methXsort:

1. Clone the repository
2. Create a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Linux/Mac
   # or
   venv\Scripts\activate  # On Windows
   ```
3. Install in development mode:
   ```bash
   pip install -e .
   ```
4. Make changes to `methxsort/core.py` or `methxsort/__main__.py`
5. Changes are immediately available when running `methxsort`

## Troubleshooting

### ImportError: No module named 'methxsort'
- Ensure you've installed the package: `pip install -e .`
- Check you're in the correct Python environment

### Command not found: methxsort
- Ensure the package is installed: `pip list | grep methXsort`
- Check that your Python scripts directory is in PATH
- Try using: `python -m methxsort` instead

### Missing dependencies
- Install requirements: `pip install -r requirements.txt`
- Or let pip handle it: `pip install -e .`
