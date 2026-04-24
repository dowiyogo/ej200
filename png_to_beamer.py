#!/usr/bin/env python3
"""
Script para generar automáticamente una presentación Beamer de LaTeX
desde imágenes PNG en un directorio.

Uso:
    python3 png_to_beamer.py [directorio_entrada] [archivo_salida.tex]
    
Ejemplo:
    python3 png_to_beamer.py fpt_png/ presentacion_timing_detector.tex
"""

import os
import sys
from pathlib import Path
from datetime import datetime
import re

def extract_caption(filename):
    """
    Extrae información del nombre del archivo y crea un caption informativo.
    
    Patrones soportados:
    - muon_-650mm_left.png → Muón @ x = -650 mm (End-Left SiPMs, IDs 0-7)
    - muon_+390mm_right.png → Muón @ x = +390 mm (End-Right SiPMs, IDs 8-15)
    - muon_+0mm_top.png → Muón @ x = 0 mm (Top SiPMs, IDs 16-35)
    """
    
    # Patrón: muon_[posición]_[cara]
    match = re.match(r'muon_([+-]?\d+)mm_(\w+)\.png', filename)
    
    if match:
        posicion = match.group(1)
        cara = match.group(2).lower()
        
        # Mapeo de caras a información de SiPMs
        cara_info = {
            'left': 'End-Left SiPMs (IDs 0–7)',
            'right': 'End-Right SiPMs (IDs 8–15)',
            'top': 'Top SiPMs (IDs 16–35)'
        }
        
        sipm_info = cara_info.get(cara, cara.capitalize())
        
        # Crear caption informativo
        if posicion.startswith('+'):
            caption = f"Muón @ x = +{posicion[1:]} mm — {sipm_info}"
        else:
            caption = f"Muón @ x = {posicion} mm — {sipm_info}"
        
        return caption
    
    # Si no coincide el patrón, devolver nombre limpio
    return filename.replace("_", "\\_").replace(".png", "")

def generate_beamer_presentation(input_dir, output_file, title="Timing Detector Analysis", author="SHiP Experiment"):
    """
    Genera una presentación Beamer desde imágenes PNG.
    
    Args:
        input_dir: Directorio con imágenes PNG
        output_file: Archivo .tex de salida
        title: Título de la presentación
        author: Autor de la presentación
    """
    
    # Validar directorio
    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: '{input_dir}' no es un directorio válido")
        sys.exit(1)
    
    # Obtener lista de PNGs ordenados
    png_files = sorted(input_path.glob("*.png"))
    
    if not png_files:
        print(f"Error: No hay archivos PNG en '{input_dir}'")
        sys.exit(1)
    
    print(f"Encontrados {len(png_files)} archivos PNG")
    
    # Generar LaTeX
    latex_content = []
    
    # Preámbulo
    latex_content.append(r"""\documentclass{beamer}
\usepackage[utf8]{inputenc}
\usepackage[spanish]{babel}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subcaption}

\usetheme{Madrid}
\usecolortheme{default}

% Configuración del título
""")
    
    latex_content.append(f"\\title{{{title}}}\n")
    latex_content.append(f"\\author{{{author}}}\n")
    latex_content.append(f"\\date{{\\today}}\n")
    
    latex_content.append(r"""
\AtBeginSection[]
{
  \begin{frame}
    \frametitle{Contenido}
    \tableofcontents[currentsection]
  \end{frame}
}

\begin{document}

\begin{frame}
  \titlepage
\end{frame}

\begin{frame}
  \frametitle{Contenido}
  \tableofcontents
\end{frame}

\section{Análisis de Fotones}

""")
    
    # Crear diapositivas para cada PNG
    for idx, png_file in enumerate(png_files, 1):
        filename = png_file.name
        
        # Extraer información del nombre: muon_[posición]_[cara].png
        caption = extract_caption(filename)
        
        # Obtener ruta relativa del archivo
        rel_path = os.path.relpath(png_file, Path(output_file).parent)
        
        latex_content.append(f"""\\begin{{frame}}[fragile]
  \\frametitle{{Análisis {idx}}}
  \\begin{{center}}
    \\includegraphics[width=0.90\\textwidth,height=0.80\\textheight,keepaspectratio]{{{rel_path}}}
    
    \\vspace{{0.1cm}}
    \\small \\texttt{{{caption}}}
  \\end{{center}}
\\end{{frame}}

""")
    
    # Cierre del documento
    latex_content.append(r"""\end{document}
""")
    
    # Escribir archivo
    output_path = Path(output_file)
    output_path.write_text("".join(latex_content), encoding='utf-8')
    
    print(f"\n✓ Presentación generada: {output_file}")
    print(f"  - {len(png_files)} diapositivas creadas")
    print(f"\nPara compilar:")
    print(f"  pdflatex -interaction=nonstopmode {output_file}")
    print(f"  # O si prefieres xelatex/lualatex para mejor soporte UTF-8:")
    print(f"  xelatex {output_file}")

def main():
    if len(sys.argv) < 2:
        print("Uso: python3 png_to_beamer.py <directorio_entrada> [archivo_salida.tex]")
        print("\nEjemplo:")
        print("  python3 png_to_beamer.py fpt_png/ presentacion.tex")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "presentacion_timing_detector.tex"
    
    generate_beamer_presentation(input_dir, output_file)

if __name__ == "__main__":
    main()
