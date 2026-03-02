# RNAhub - RNA Homology Search Pipeline

RNAhub is a Python-based bioinformatics pipeline for discovering conserved structural RNAs through iterative homology searching using HMMER, followed by structural analysis with R-scape and INFERNAL. The pipeline processes FASTA sequences through multiple search iterations to build covariance models for RNA structure prediction.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Bootstrap and Environment Setup
- Install Python dependencies:
  - `python -m pip install icecream biopython` -- Required Python packages. NEVER CANCEL: Takes 2-3 minutes.
- Install system packages:
  - `sudo apt-get update && sudo apt-get install -y figlet hmmer infernal` -- NEVER CANCEL: Takes 3-5 minutes.
- Create local configuration:
  ```bash
  cat > config_local.py << 'EOF'
  # Local configuration overrides
  nhmmer = "/usr/bin/nhmmer"
  INFERNAL_PATH = '/usr/bin/'
  EASEL_PATH = '/usr/bin/'  # Warning: EASEL tools not available in Ubuntu packages
  RSCAPE_PATH = "/usr/bin/R-scape-NOT-INSTALLED"  # Not available in standard repos
  REPEAT_MASKER_PATH = '/usr/bin/RepeatMasker-NOT-INSTALLED'  # Not available
  RFAM_DB_PATH = "/dev/null/Rfam.cm"  # Database not included
  RFAM_FILES = "/dev/null/rfam"
  EOF
  ```

### Basic Testing and Validation
- Test environment setup:
  - `python rnahub.py --version` -- Shows version info, completes in <1 second
  - `python rnahub.py --help` -- Shows all command options, completes in <1 second
- Run dry-run test:
  - `python rnahub.py --dry --input example/u6.fa --db example/example.fa --job-folder /tmp/test --dev-skip-rscape --dev-skip-infernal` -- Shows commands without execution, completes in <1 second
- Test core functionality (partial):
  - `python rnahub.py --input example/u6.fa --db example/example.fa --job-folder /tmp/test --dev-skip-rscape --dev-skip-infernal --cpus 1` -- NEVER CANCEL: May take 2-5 minutes. Will fail at EASEL tools but tests nhmmer.

### Full Pipeline Requirements (External Tools)
CRITICAL: The complete pipeline requires external bioinformatics tools not available in standard package repositories:

**Available via apt packages:**
- nhmmer (HMMER 3.4) - Available via `apt install hmmer`
- INFERNAL tools (cmalign, cmbuild, cmsearch) - Available via `apt install infernal`

**NOT available in Ubuntu packages - require manual installation:**
- EASEL tools (esl-alimanip, esl-alistat, esl-reformat) - Must build from HMMER source
- R-scape - Must download and build from http://rivaslab.org/software.html
- RepeatMasker - Must download and install separately

**Required databases (not included):**
- Large genome databases for homology searches (typically 1-50 GB)
- Rfam covariance model database (Rfam.cm)

## Core Workflows

### Main Pipeline (rnahub.py)
- Basic usage: `python rnahub.py --input <fasta_file> --db <genome_database> --job-folder <output_dir>`
- With flanking sequences: `python rnahub.py --input <query.fa> --flanked <flanked.fa> --db <database> --job-folder <output>`
- Skip expensive steps during development:
  - `--dev-skip-search` - Skip all nhmmer iterations
  - `--dev-skip-rscape` - Skip R-scape structural analysis  
  - `--dev-skip-infernal` - Skip INFERNAL covariance model building
  - `--dry` - Show commands without execution

### Supporting Tools
- BLAST search: `python search_blast.py --db <database> <fasta_file>` - Requires BLAST+ installation and databases
- Rfam search: `python search_rfam.py <fasta_file>` - Requires Rfam database and installation

### Development and Testing Options
- Use `--dry` flag to see all commands without execution - ALWAYS use this first to validate command structure
- Use skip flags (`--dev-skip-*`) to test individual pipeline components
- Use small example files in `example/` directory for initial testing
- Set `--cpus 1` for testing to avoid resource conflicts

## Timing Expectations and Critical Warnings

### Short Operations (<1 minute)
- Help and version commands
- Dry runs and command validation
- Job folder creation and file copying

### Medium Operations (2-10 minutes)
- **NEVER CANCEL**: Python package installation (pip install)
- **NEVER CANCEL**: System package installation (apt install)
- **NEVER CANCEL**: Single nhmmer iteration with small databases and queries

### Long Operations (30+ minutes)
- **NEVER CANCEL**: Complete pipeline runs with large genome databases - Set timeout to 120+ minutes
- **NEVER CANCEL**: R-scape structural analysis on large alignments - Set timeout to 60+ minutes  
- **NEVER CANCEL**: Multiple nhmmer iterations with large databases - Set timeout to 180+ minutes

**CRITICAL**: Bioinformatics pipelines commonly run for hours. NEVER cancel commands that appear to hang - they may be processing large datasets. Always set explicit timeouts of 60+ minutes for build and analysis commands.

## Validation

### After Environment Setup
Always validate the installation works before making changes:
1. Run `python rnahub.py --version` - should complete without import errors
2. Run dry-run test with example data
3. Test partial pipeline with available tools
4. Check job folder creation and file copying works correctly

### Manual Testing Scenarios
- **Basic functionality**: Use example files to test job creation, file copying, and command generation
- **Tool availability**: Verify which external tools are available vs missing
- **Configuration**: Test config_local.py overrides work correctly
- **Error handling**: Verify graceful failure when external tools or databases are missing

### Common Validation Commands
Always run these before committing changes:
- `python -c "import icecream, Bio.SeqIO; print('Dependencies OK')"` - Verify Python dependencies
- `python rnahub.py --dry --input example/u6.fa --db example/example.fa --job-folder /tmp/validate` - Test command generation
- `ls /tmp/validate/u6/` - Verify job folder structure creation

## Repository Structure

### Key Files
- `rnahub.py` - Main pipeline script (770 lines)
- `config.py` - Default configuration with tool paths
- `config_local.py` - Local overrides (create this file)
- `search_blast.py` - BLAST search utility
- `search_rfam.py` - Rfam database search utility
- `test.sh` - Test script for development environment

### Important Directories
- `example/` - Test FASTA files for validation
- `opt/` - Expected location for manually installed tools (not present in sandboxed environment)
- `jobs/` - Default output directory for pipeline results

### Configuration Files
- `config.py` - Contains default paths to external tools and databases
- `config_local.py` - Override file for local installations (import order: config.py, then config_local.py)

## Common Tasks and Troubleshooting

### Installation Issues
- **Import errors**: Install missing Python packages with `python -m pip install <package>`
- **Tool not found**: Check if tool needs manual installation vs apt package
- **Permission errors**: Use `sudo` for system package installation

### Pipeline Failures
- **No hits found**: Expected behavior when query doesn't match database
- **File not found errors**: Often due to missing external tools - check tool availability
- **Database errors**: Requires specific bioinformatics databases not included in repository

### External Tool Status
The repository assumes a full bioinformatics environment with manually compiled tools. In a fresh environment:
- **Working**: Basic Python functionality, job management, command generation
- **Partial**: nhmmer and INFERNAL tools via apt packages
- **Missing**: EASEL tools, R-scape, RepeatMasker, large sequence databases

### Development Workflow
1. Always test with `--dry` flag first
2. Use development skip flags to test components individually  
3. Start with small example files before using large databases
4. Check job folder outputs to understand pipeline state
5. Use `--verbose` flag for detailed logging when debugging

## Example Command Outputs

### Repository Root Structure
```
.
├── README.md
├── config.py
├── config_local.py  # Create this
├── example/
│   ├── u6.fa
│   ├── example.fa
│   └── xrrna.fa
├── rnahub.py
├── search_blast.py
├── search_rfam.py
└── test.sh
```

### Successful Dry Run Output
```
python rnahub.py --dry --input example/u6.fa --db example/example.fa --job-folder /tmp/test
# Shows: version info, argument parsing, command generation without execution
# Expected runtime: <1 second
```

### Partial Pipeline Test
```
python rnahub.py --input example/u6.fa --db example/example.fa --job-folder /tmp/test --dev-skip-rscape --dev-skip-infernal
# Shows: nhmmer execution, tool failures at missing EASEL commands
# Expected runtime: 2-5 minutes
```