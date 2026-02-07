# FFT-DSP Toolkit: Repository Structure & .gitignore Guide

## What Should Be in GitHub ✓

### Source Code
- `src/*.c` - Implementation files
- `include/*.h` - Header files (public API)
- `tests/*.c` - Test cases and test framework
- `examples/*.c` - Example programs

### Configuration & Build
- `Makefile` - Build system
- `CMakeLists.txt` - CMake configuration
- `.github/workflows/*.yml` - CI/CD pipelines

### Documentation
- `*.md` files (README, architecture, API, etc.)
- `docs/diagrams/*.puml` - PlantUML source files (version control these!)
- `docs/images/` - Static images created by hand

### Project Files
- `.gitignore` - Repository configuration
- `LICENSE` - License file
- `.github/` - GitHub-specific config (issues, PR templates)

---

## What Should NOT Be in GitHub ✗

### Build Artifacts
```
build/                 # Compiled output
bin/                   # Executables
lib/                   # Libraries (.a, .so)
*.o                    # Object files
*.a                    # Static libraries
*.so, *.dylib, *.dll   # Shared libraries
CMakeFiles/            # CMake cache
cmake_install.cmake    # CMake generated
CMakeCache.txt         # CMake cache
```

**Reason**: Everyone rebuilds locally with their toolchain. Binary compatibility varies.

### Generated Files (Can Regenerate)
```
docs/diagrams/*.png    # Generated from .puml with plantuml
coverage_html/         # Generated from coverage.info
perf.data              # Profiling data (varies per machine)
*.gcda, *.gcno         # Coverage files
```

**Reason**: These are artifacts; keep the source (.puml, .c) instead.

**To Regenerate PNGs**:
```bash
cd docs/diagrams
python3 plantuml_render.py        # Online service
# OR
plantuml *.puml                   # Local installation
# OR
docker run --rm -v $(pwd):/d plantuml/plantuml /d/*.puml
```

### IDE & Editor Files
```
.vscode/               # VS Code settings (personal)
.idea/                 # IntelliJ settings
*.swp, *.swo           # Vim temp files
*~                     # Emacs backups
.DS_Store              # macOS
Thumbs.db              # Windows
```

**Reason**: Every developer has different IDE preferences. Use `.editorconfig` instead.

### Test & Coverage Output
```
coverage/              # Coverage reports
*.gcda, *.gcno         # GCC coverage files
coverage.info          # LCOV data
.pytest_cache/         # Python test cache
```

**Reason**: Regenerated on each test run. Only commit code, not reports.

### Profiling & Debug Output
```
perf.data              # Linux perf profiler output
perf.data.old          # Old perf data
*.prof                 # Python profiling
callgrind.out.*        # Valgrind profiling
```

**Reason**: Machine-specific; varies per run.

### System Files
```
.DS_Store              # macOS index
Thumbs.db              # Windows thumbnails
__pycache__/           # Python cache
*.pyc, *.pyo           # Python compiled
```

**Reason**: Not part of the project; system-generated.

### Logs, Temp Files
```
*.log                  # Build/test logs
*.tmp, *.bak           # Temporary files
output/                # Various output dirs
results/               # Test results
.cache/                # Cache directories
```

**Reason**: Ephemeral; regenerated on demand.

### Security Files
```
token.txt              # GitHub tokens
secrets.h              # API keys
*.key, *.pem           # Certificates
```

**Reason**: NEVER commit secrets. Use environment variables instead.

---

## Repository Size Analysis

Before committing, check:

```bash
# Size of repository
du -sh .git                      # Repository size
git ls-files | head -20          # Tracked files
git ls-files | wc -l             # Total tracked

# Find large files
find . -type f -size +1M         # > 1MB files
git rev-list --all --objects | sort -k2 | tail -10
```

### Size Goals
- **Target**: < 10 MB (source + docs only)
- **With PNG diagrams**: ~15-20 MB (acceptable)
- **With build artifacts**: 100+ MB (NEVER)

---

## Recommended Additions to .gitignore

Your current `.gitignore` is comprehensive. Key entries:

```gitignore
# Build artifacts (always)
build/
*.o
*.a
*.so

# Generated diagrams (regenerate as needed)
docs/diagrams/*.png

# IDE (everyone has different preferences)
.vscode/
.idea/

# Security (CRITICAL)
token.txt
*.key
secrets.h

# Test outputs (regenerate on test run)
coverage/
*.gcda
```

---

## Workflow: Before Committing

```bash
# 1. Check what will be committed
git status

# 2. See diffs
git diff

# 3. Verify no secrets
grep -r "password\|token\|key\|secret" src/ include/

# 4. Clean build artifacts
make clean
make distclean

# 5. Commit only source + docs
git add src/ include/ tests/ examples/ docs/ *.md Makefile CMakeLists.txt
git commit -m "Feature: Add windowing functions"
git push origin main
```

---

## When Diagrams Change

PlantUML source files (`.puml`) SHOULD be committed. PNG images can be:

- Regenerated locally on demand
- Kept in repo for convenience (15-20 MB is acceptable)
- Or excluded and regenerated by CI/CD

### Option A: Exclude PNG, Regenerate on Demand
```gitignore
docs/diagrams/*.png    # Exclude - regenerate locally
```

```bash
# Developer regenerates when needed:
cd docs/diagrams && python3 plantuml_render.py
```

### Option B: Include PNG for Convenience
```gitignore
# (Remove PNG exclusion line)
```

```bash
# Make sure to regenerate before committing:
python3 plantuml_render.py
git add docs/diagrams/*.png
```

**Recommendation**: **Exclude PNGs** (Option A). They're generated from `.puml` files; no need to version control them.

---

## CI/CD Consideration

GitHub Actions can regenerate diagrams automatically:

```yaml
# .github/workflows/generate-docs.yml
- name: Generate diagrams
  run: |
    cd docs/diagrams
    python3 plantuml_render.py
    git add *.png
    git commit -m "Auto-generated diagrams" || true
    git push
```

---

## Current Status

Your `.gitignore` now includes:
- ✓ Build artifacts
- ✓ Generated PNG files
- ✓ IDE/editor files
- ✓ Test coverage output
- ✓ Profiling data
- ✓ Security files (tokens, keys)
- ✓ System files
- ✓ Cache directories

### What's OK to Commit
- ✓ `*.puml` files (source diagrams)
- ✓ `*.md` documentation
- ✓ `*.c` and `*.h` source code
- ✓ Configuration files (Makefile, CMakeLists.txt, etc.)

### What Will Be Ignored
- ✗ `docs/diagrams/*.png` (regenerated)
- ✗ `build/`, `bin/`, `lib/` (build artifacts)
- ✗ `.vscode/`, `.idea/` (IDE settings)
- ✗ `*.log`, `*.tmp` (temporary)
- ✗ `token.txt` (security)

---

## Files to Clean Before First Push

```bash
cd ~/Documents/fft-dsp-toolkit

# Remove generated files NOT in .gitignore yet
rm -f docs/diagrams/*.png      # Regenerated from .puml
rm -rf build/                  # Build artifacts
rm -rf .vscode/                # IDE settings (if present)

# Verify what will be committed
git status
git ls-files | head -20
```

---

## See Also

- [GitHub .gitignore templ ates](https://github.com/github/gitignore)
- [C language .gitignore](https://github.com/github/gitignore/blob/main/C.gitignore)
- [Best practices](https://git-scm.com/docs/gitignore)

---

**Summary**: Commit source & docs, ignore build artifacts & generated files.
