# Files to Commit vs Skip

## Summary

**DO COMMIT** ✓ (Essential)
- Source code: `src/`, `include/`, `tests/`, `examples/`
- Configuration: `Makefile`, `CMakeLists.txt`, `.github/`
- Documentation: `*.md` files  
- Diagrams source: `*.puml` files
- Build system & project files

**DO NOT COMMIT** ✗ (Can regenerate)
- PNG diagrams: `docs/diagrams/*.png`
- Build artifacts: `build/` directory
- Test file: `TEST_PUSH.txt` (temporary test)

---

## Files by Category

### ✓ COMMIT (Source & Configuration)
```
.github/workflows/ci.yml          CI/CD pipeline config
Makefile                          Build system
CMakeLists.txt                    CMake configuration
.gitignore                        Repository rules
PROJECT_EXPANSION_PLAN.md         Development roadmap
README.md                         Project overview

src/fft.c                         FFT implementation
src/filter.c                      Filter implementation
src/dsp_utils.c                   Utility functions

include/fft.h                     FFT header
include/filter.h                  Filter header
include/dsp_utils.h               Utils header

tests/test_fft.c                  FFT tests
tests/test_filter.c               Filter tests
tests/test_framework.h            Test framework

examples/fft_demo.c               Example program
examples/filter_demo.c            Example program
```

### ✓ COMMIT (Documentation)
```
docs/API.md                       Complete API reference
docs/ARCHITECTURE.md              System architecture
docs/GITIGNORE_GUIDE.md           Repository guidelines
docs/README.md                    Docs overview
docs/dsp_references.md            Reference materials
docs/fft_notes.md                 FFT notes

docs/diagrams/*.puml              PlantUML source (VERSION CONTROL!)
docs/diagrams/README.md           Diagram documentation
docs/diagrams/plantuml_render.py  PNG generator script
docs/diagrams/render_diagrams.sh  Rendering helper
```

### ✗ DON'T COMMIT (Generated)
```
docs/diagrams/*.png               Generated from .puml (regenerate locally)
build/                            Compiled artifacts (make clean)
TEST_PUSH.txt                     Temporary test file (remove)
```

---

## Clean & Commit Steps

```bash
cd ~/Documents/fft-dsp-toolkit

# 1. Clean temporary test file
rm -f TEST_PUSH.txt

# 2. Clean build artifacts
make clean
rm -rf build

# 3. Don't add PNG files (regenerated)
# Already in .gitignore

# 4. Add new documentation
git add docs/API.md
git add docs/ARCHITECTURE.md
git add docs/GITIGNORE_GUIDE.md
git add docs/diagrams/*.puml
git add docs/diagrams/*.sh
git add docs/diagrams/*.py
git add docs/diagrams/README.md

# 5. Update .gitignore
git add .gitignore

# 6. View what will be committed
git status

# 7. Commit
git commit -m "docs: Add comprehensive visual documentation (10 PlantUML diagrams, API & architecture guides)"

# 8. Push
git push origin main
```

---

## File Size Summary

### Before Cleanup
- Repository Size: ~150 KB
- Largest items:
  - `.git/` folder: ~100+ KB
  - `build/` directory: 100 KB
  - `docs/diagrams/*.png`: ~110 KB (10 files total)

### After Cleanup (Recommended)
- Keep: Source code + PlantUML files
- Remove: PNG files, build artifacts, test files
- Expected: ~50-70 KB repository size

### With PNG Diagrams (Optional)
- If keeping PNG files: ~120-150 KB
- Pros: No need to regenerate locally
- Cons: More storage; less important to version control

---

## Recommendation

**Exclude PNG files from Git:**
- They're generated from `.puml` source
- Take up space
- Can be regenerated anytime with:
  ```bash
  cd docs/diagrams && python3 plantuml_render.py
  ```

**Include everything else:**
- PlantUML sources (version control the source!)
- All `.md` documentation
- Rendering scripts

---

## Next Steps

1. **Clean repo:**
   ```bash
   rm -f TEST_PUSH.txt
   make clean
   ```

2. **Verify .gitignore:**
   ```bash
   git status        # Should NOT show *.png files
   ```

3. **Add new files:**
   ```bash
   git add docs/ .gitignore
   ```

4. **Commit & push:**
   ```bash
   git commit -m "docs: Comprehensive visual & API documentation"
   git push origin main
   ```

5. **Verify on GitHub:**
   - Check that `.puml` files are there
   - Check that `.png` files are NOT there (in .gitignore)
   - Verify file structure is clean

---

**Checklist**: 
- [ ] Remove TEST_PUSH.txt
- [ ] Run make clean
- [ ] Verify .gitignore excludes *.png
- [ ] Add all new *.md and *.puml files
- [ ] Commit with clear message
- [ ] Verify GitHub reflects changes
