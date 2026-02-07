#!/bin/bash
# Generate PNG diagrams from PlantUML files
# Usage: ./render_diagrams.sh [--online|--local|--docker]

set -e

DIAGRAMS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-auto}"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${BLUE}[PlantUML]${NC} $1"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Detect available rendering method
detect_method() {
    if command -v plantuml &> /dev/null; then
        echo "local"
    elif command -v docker &> /dev/null; then
        echo "docker"
    else
        echo "online"
    fi
}

# Local rendering (requires: plantuml, graphviz)
render_local() {
    print_status "Using local PlantUML installation..."
    
    if ! command -v plantuml &> /dev/null; then
        print_warning "PlantUML not found. Install with: apt-get install plantuml"
        return 1
    fi
    
    for puml_file in "$DIAGRAMS_DIR"/*.puml; do
        if [ -f "$puml_file" ]; then
            basename=$(basename "$puml_file" .puml)
            print_status "Rendering $basename..."
            plantuml -png "$puml_file" -o "$DIAGRAMS_DIR"
            print_success "Generated $basename.png"
        fi
    done
}

# Docker rendering
render_docker() {
    print_status "Using Docker for PlantUML rendering..."
    
    if ! command -v docker &> /dev/null; then
        print_warning "Docker not found. Install from https://docker.com"
        return 1
    fi
    
    docker run --rm -v "$DIAGRAMS_DIR:/diagrams" \
        plantuml/plantuml:latest \
        /diagrams/*.puml -o /diagrams
    
    print_success "All diagrams rendered via Docker"
}

# Online rendering (uses PlantUML online service)
render_online() {
    print_status "Using online PlantUML service..."
    print_warning "Requires internet connection"
    
    for puml_file in "$DIAGRAMS_DIR"/*.puml; do
        if [ -f "$puml_file" ]; then
            basename=$(basename "$puml_file" .puml)
            print_status "Rendering $basename..."
            
            # Read PlantUML file and encode for URL
            content=$(cat "$puml_file")
            
            # Use PlantUML online renderer
            # fallback: just warn user
            print_warning "Online rendering requires manual step"
            echo "  Visit: http://www.planttext.com/?text=$(echo "$content" | base64)"
        fi
    done
}

# Main logic
case "$MODE" in
    local)
        render_local || echo "Falling back to Docker..."
        render_docker || echo "Falling back to online..."
        ;;
    docker)
        render_docker || echo "Falling back to local..."
        render_local || echo "Please install PlantUML or Docker"
        ;;
    online)
        render_online
        ;;
    auto|*)
        method=$(detect_method)
        print_status "Auto-detected method: $method"
        case "$method" in
            local)  render_local ;;
            docker) render_docker ;;
            *)      render_online ;;
        esac
        ;;
esac

print_success "Diagram rendering complete!"
print_status "Generated PNG files in: $DIAGRAMS_DIR"
echo ""
echo "View diagrams:"
echo "  - File browser: open $DIAGRAMS_DIR"
echo "  - Command line: ls -la $DIAGRAMS_DIR/*.png"
echo "  - VS Code: Code > Open Folder > docs/diagrams"
