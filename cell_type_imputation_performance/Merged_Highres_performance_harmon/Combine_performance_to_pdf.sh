
<<comment

cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/cell_type_imputation_performance/Merged_Highres_performance_harmon/

chmod +x Combine_performance_to_pdf.sh

./Combine_performance_to_pdf.sh

comment

#!/usr/bin/env bash
set -euo pipefail

# Combine all PNGs in this folder into a single bookmarked PDF with filename annotations.
# Usage:
#   ./Combine_performance_to_pdf.sh                # creates Combined.pdf
#   ./Combine_performance_to_pdf.sh output.pdf     # custom output name

out_pdf="${1:-Combined.pdf}"
tmp_pdf="$(mktemp -t combined_pngs.XXXXXX).pdf"
bookmarks_file="$(mktemp -t pdf_bookmarks.XXXXXX).txt"
annot_dir="annotated"

# --- dependencies ---
for bin in convert pdftk; do
  command -v "$bin" >/dev/null 2>&1 || {
    echo "Error: '$bin' is required but not found in PATH." >&2
    exit 1
  }
done

# --- collect PNGs (case-insensitive), natural sort like Fig_S1, Fig_S2, ... ---
# If you need strict lexicographic order, replace 'sort -zV' with 'sort -z'
mapfile -d '' files < <(find . -maxdepth 1 -type f \( -iname '*.png' \) -print0 | sort -zV)

if (( ${#files[@]} == 0 )); then
  echo "No PNG files found in $(pwd)" >&2
  exit 1
fi

# --- annotate each image with its base filename into ./annotated ---
mkdir -p "$annot_dir"

# Increase density if you want sharper PDFs (e.g., -density 300 before input)
# Adjust -pointsize as desired.
for f in "${files[@]}"; do
  base="$(basename "$f")"
  convert "$f" \
    -gravity north -background white -splice 0x100 \
    -gravity north -pointsize 70 -annotate +0+40 "$base" \
    "$annot_dir/$base"
done

# --- build the combined PDF from annotated images, preserving order ---
# Pass files explicitly to avoid glob reordering.
# You can tweak density/compression here if needed.
# shellcheck disable=SC2206
annot_list=()
for f in "${files[@]}"; do
  annot_list+=("$annot_dir/$(basename "$f")")
done

# Create temporary combined PDF (no bookmarks yet)
convert "${annot_list[@]}" "$tmp_pdf"

# --- create PDF bookmarks matching the page order ---
# Use UTF-8 info to safely handle non-ASCII names.
: > "$bookmarks_file"
page=1
for f in "${files[@]}"; do
  title="$(basename "$f")"
  {
    echo "BookmarkBegin"
    echo "BookmarkTitle: $title"
    echo "BookmarkLevel: 1"
    echo "BookmarkPageNumber: $page"
  } >> "$bookmarks_file"
  ((page++))
done

pdftk "$tmp_pdf" update_info_utf8 "$bookmarks_file" output "$out_pdf"

# --- cleanup ---
rm -f "$tmp_pdf" "$bookmarks_file"

echo "✅ Created '$out_pdf' with $((${#files[@]})) pages and bookmarks."
echo "   Annotated images are in './$annot_dir'. You can delete that folder if not needed."

