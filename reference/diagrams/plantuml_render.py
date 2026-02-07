#!/usr/bin/env python3
"""
PlantUML to PNG converter using online service.
Generates PNG diagrams from .puml files.
"""

import os
import sys
import urllib.parse
import urllib.request
import base64
from pathlib import Path

def plantuml_to_png(puml_file, output_file):
    """Convert PlantUML file to PNG using online service."""
    
    # Read PlantUML source
    with open(puml_file, 'r') as f:
        puml_content = f.read()
    
    # Remove @startuml/@enduml markers for encoding
    lines = puml_content.split('\n')
    filtered = [l for l in lines if not l.strip().startswith('@')]
    filtered_content = '\n'.join(filtered)
    
    # Encode for PlantUML URL
    # Method: compress with zlib and encode in base64url format
    import zlib
    
    compressed = zlib.compress(filtered_content.encode('utf-8'))
    # PlantUML uses base64url encoding (but with normal base64 chars)
    encoded = base64.b64encode(compressed).decode('ascii')
    
    # Build PlantUML rendering URL
    url = f"https://www.plantuml.com/plantuml/png/{encoded}"
    
    print(f"  Rendering {Path(puml_file).name}...")
    
    try:
        # Download PNG from PlantUML service
        response = urllib.request.urlopen(url, timeout=10)
        png_data = response.read()
        
        # Save PNG file
        with open(output_file, 'wb') as f:
            f.write(png_data)
        
        size_kb = len(png_data) / 1024
        print(f"    ✓ Generated {Path(output_file).name} ({size_kb:.1f} KB)")
        return True
        
    except Exception as e:
        print(f"    ✗ Error: {e}")
        # Create placeholder PNG
        create_placeholder_png(output_file, Path(puml_file).stem)
        return False

def create_placeholder_png(output_file, name):
    """Create a simple placeholder PNG."""
    # Minimal PNG file (1x1 white pixel)
    # PNG header + IHDR + IDAT + IEND chunks
    png_minimal = bytes([
        0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A,  # PNG signature
        0x00, 0x00, 0x00, 0x0D,                           # IHDR chunk size
        0x49, 0x48, 0x44, 0x52,                           # IHDR
        0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,  # 1x1 dimensions
        0x08, 0x02, 0x00, 0x00, 0x00,                     # 8-bit RGB
        0x90, 0x77, 0x53, 0xDE,                           # CRC
        0x00, 0x00, 0x00, 0x0C,                           # IDAT chunk size
        0x49, 0x44, 0x41, 0x54,                           # IDAT
        0x08, 0x99, 0x63, 0xF8, 0xFF, 0xFF, 0x3F, 0x00,  # Compressed white data
        0x00, 0x05, 0x00, 0x01,
        0x0D, 0x0A, 0x2D, 0xB4,                           # CRC
        0x00, 0x00, 0x00, 0x00,                           # IEND chunk size
        0x49, 0x45, 0x4E, 0x44,                           # IEND
        0xAE, 0x42, 0x60, 0x82,                           # CRC
    ])
    
    with open(output_file, 'wb') as f:
        f.write(png_minimal)
    
    print(f"    ℹ Created placeholder: {Path(output_file).name}")

def main():
    diagrams_dir = Path(__file__).parent
    
    puml_files = sorted(diagrams_dir.glob('*.puml'))
    
    if not puml_files:
        print("No .puml files found in current directory")
        return 1
    
    print(f"Found {len(puml_files)} PlantUML files")
    print("Rendering diagrams (using online PlantUML service)...\n")
    
    success_count = 0
    for puml_file in puml_files:
        output_file = puml_file.with_suffix('.png')
        if plantuml_to_png(puml_file, output_file):
            success_count += 1
    
    print(f"\n✓ Generated {success_count}/{len(puml_files)} diagrams")
    print(f"Output directory: {diagrams_dir}")
    
    # List generated files
    png_files = sorted(diagrams_dir.glob('*.png'))
    if png_files:
        print(f"\nGenerated PNG files:")
        for png in png_files:
            size = png.stat().st_size / 1024
            print(f"  • {png.name} ({size:.1f} KB)")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
