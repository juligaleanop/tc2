#!/usr/bin/env bash
set -euo pipefail
 
shopt -s nullglob          # si no hay .pdf, el loop no ejecuta nada
pdfs=( ./*.pdf )
 
if [[ ${#pdfs[@]} -eq 0 ]]; then
    echo "No se encontraron archivos .pdf en el directorio actual."
    exit 1
fi
 
for pdf in "${pdfs[@]}"; do
    base="${pdf%.pdf}"          # ej: "./documento" a partir de "./documento.pdf"
    nombre="${base#./}"         # ej: "documento" (sin el ./)
 
    echo "── Procesando: $nombre.pdf"
 
    # 1. pdfcrop: sobreescribe el mismo archivo
    pdfcrop --margins 10 "$nombre.pdf" "$nombre.pdf"
 
    # 2. pdftoppm: genera <nombre>-1.png (para la primera página)
    pdftoppm -png "$nombre.pdf" "$nombre"
 
    # 3. Renombrar <nombre>-1.png → <nombre>.png
    if [[ -f "${nombre}-1.png" ]]; then
        mv "${nombre}-1.png" "${nombre}.png"
        echo "   ✓ Imagen guardada: $nombre.png"
    else
        echo "   ✗ No se encontró ${nombre}-1.png — verificá si el PDF tiene páginas."
    fi
done
 
echo ""
echo "Listo. Se procesaron ${#pdfs[@]} archivo(s)."
