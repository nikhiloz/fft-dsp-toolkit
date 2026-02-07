# DSP Tutorial Suite: Visual Documentation

Complete PlantUML diagrams for the DSP Tutorial Suite project. All diagrams are automatically rendered to PNG for easy viewing.

## Diagram Index

### 1. **architecture.puml** - System Architecture
- Overview of the complete DSP Tutorial Suite
- Application layer, core library organization, and platform abstraction
- Shows how components interact and their relationships
- **Purpose**: Understanding the high-level system design

> [🔍 View System Architecture diagram →](architecture.png)

### 2. **signal_flow.puml** - Signal Processing Data Flow
- Typical signal processing pipeline from input to output
- Time domain and frequency domain operations
- Windowing, FFT, filtering, and post-processing stages
- **Purpose**: Understanding the data transformation pipeline

> [🔍 View Signal Processing Pipeline diagram →](signal_flow.png)

### 3. **modules.puml** - Module Dependencies
- Detailed module structure and interdependencies
- Shows how source files depend on each other
- Highlights modularity and separation of concerns
- **Purpose**: Understanding code organization and dependencies

> [🔍 View Module Dependencies diagram →](modules.png)

### 4. **fft_sequence.puml** - FFT Processing Sequence
- Step-by-step call sequence for FFT operations
- Shows interaction between application, buffers, and SIMD kernels
- Includes real-time streaming with overlap-add
- **Purpose**: Understanding runtime execution flow

> [🔍 View FFT Processing Sequence diagram →](fft_sequence.png)

### 5. **realtime_architecture.puml** - Real-Time Streaming Architecture
- Complete real-time audio processing pipeline
- Ring buffers, synchronization, and priority scheduling
- Low-latency design for audio/sensor applications
- **Purpose**: Understanding real-time constraints and design patterns

> [🔍 View Real-Time Streaming Architecture diagram →](realtime_architecture.png)

### 6. **optimization_roadmap.puml** - Performance Optimization Strategy
- 5-stage optimization approach (Compiler → Algorithm → SIMD → Threading → Platform)
- Expected speedups at each stage
- Total potential for 16-128x improvement
- **Purpose**: Understanding performance enhancement strategy

> [🔍 View Optimization Roadmap diagram →](optimization_roadmap.png)

### 7. **api_reference.puml** - Public API Structure
- All public functions organized by module
- Function signatures at a glance
- Relationship between modules
- **Purpose**: Quick API reference and module overview

> [🔍 View Public API Reference diagram →](api_reference.png)

### 8. **roadmap.puml** - Project Development Roadmap
- 6-phase project timeline
- Current status and planned milestones
- Duration and key deliverables for each phase
- **Purpose**: Project planning and progress tracking

> [🔍 View Development Roadmap diagram →](roadmap.png)

### 9. **benchmarks.puml** - Performance Benchmarks
- Latency comparison: dsp_core vs competitors
- Shows baseline, SIMD-optimized, and real-time kernel versions
- Comparison with FFTW, GSL, NumPy, MATLAB, Eigen
- **Purpose**: Performance goals and competitive analysis

> [🔍 View Performance Benchmarks diagram →](benchmarks.png)

### 10. **use_cases.puml** - Primary Use Cases
- Audio processing, embedded systems, research applications
- Real-time applications categories
- Actor roles and their interactions
- **Purpose**: Understanding target market and applications

> [🔍 View Use Cases diagram →](use_cases.png)

## Rendering Diagrams

### Option 1: Online (No Installation)
Visit [PlantUML Online Editor](http://planttext.com/) or use the VS Code extension:
```bash
code --install-extension jebbs.plantuml
```

Then open any `.puml` file and press `Alt+D` to preview.

### Option 2: Local Rendering (Linux/macOS)
Install PlantUML and Graphviz:
```bash
# Ubuntu/Debian
sudo apt-get install plantuml graphviz

# macOS
brew install plantuml graphviz

# Then render all diagrams:
cd docs/diagrams
for f in *.puml; do plantuml "$f"; done
```

### Option 3: Docker (Cross-Platform)
```bash
docker run --rm -v $(pwd):/workspace plantuml/plantuml \
  reference/diagrams/*.puml -o /workspace/docs/diagrams
```

## Directory Structure

```
reference/diagrams/
├── *.puml          (PlantUML source files)
├── *.png           (Generated PNG images)
└── README.md       (This file)
```

## Usage in Documentation

All diagrams are embedded in the main project documentation:
- `reference/ARCHITECTURE.md` - System design and architecture
- `reference/PERFORMANCE.md` - Optimization strategies and benchmarks
- `reference/API.md` - Public API reference
- `reference/REALTIME.md` - Real-time programming guide

## Updating Diagrams

1. Edit the `.puml` file with your changes
2. Regenerate PNG using PlantUML
3. Update corresponding documentation
4. Commit both `.puml` and `.png` files

## Quality Standards

- **Clarity**: Each diagram focuses on one aspect
- **Completeness**: Show relationships and dependencies
- **Consistency**: Use standard shapes and colors
- **Maintainability**: Keep text concise and readable
- **Accessibility**: Provide alt-text and descriptions

## Related Files

- [ARCHITECTURE.md](../ARCHITECTURE.md) - Detailed design document
- [API.md](../API.md) - Function reference

---

**Note**: All `.png` files are generated from `.puml` sources. 
Edit the `.puml` files directly; regenerate PNG on changes.
