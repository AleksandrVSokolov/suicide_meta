
<<comment

cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/

chmod +x Prepare_supplementary_figures.sh

./Prepare_supplementary_figures.sh prepared_supplementary_figures

comment
#!/bin/bash

# Usage: ./Prepare_supplementary_figures.sh /path/to/input_dir

if [ $# -lt 1 ]; then
  echo "Usage: $0 <input_dir>"
  exit 1
fi

input_dir="$1"
output_dir="$input_dir/annotated"

mkdir -p "$output_dir"

shopt -s nullglob
# Collect PNG files, extract S<number> for sorting
mapfile -t imagefiles < <(
  for f in "$input_dir"/*.png; do
    if [[ $(basename "$f") =~ _S([0-9]+) ]]; then
      num=${BASH_REMATCH[1]}
      printf "%04d %s\n" "$num" "$f"
    fi
  done | sort -n | cut -d' ' -f2-
)

if [ ${#imagefiles[@]} -eq 0 ]; then
  echo "No PNG files with _S<number> found in $input_dir"
  exit 1
fi

# Extract start and end numbers
start_num=$(basename "${imagefiles[0]}" | sed -E 's/.*_S([0-9]+).*/\1/')
end_num=$(basename "${imagefiles[-1]}" | sed -E 's/.*_S([0-9]+).*/\1/')

output_pdf="$input_dir/Fig_S${start_num}_S${end_num}.pdf"

# Annotate files in sorted order (with extra space for label)
annotated_files=()
for file in "${imagefiles[@]}"; do
  filename=$(basename "$file")
  annotated_file="$output_dir/$filename"
  
  # Get image dimensions
  width=$(identify -format "%w" "$file")
  height=$(identify -format "%h" "$file")
  
  # Add 100px of extra white space at the top, then label it
  convert "$file" -background white -gravity north \
    -splice 0x100 -pointsize 70 -annotate +0+30 "$filename" "$annotated_file"
  
  annotated_files+=("$annotated_file")
done

# Combine into a PDF
combined_pdf="$input_dir/Annotated_images_tmp.pdf"
convert "${annotated_files[@]}" "$combined_pdf"

# Create bookmarks.txt
bookmarks="$input_dir/bookmarks.txt"
rm -f "$bookmarks"
page=1
for file in "${imagefiles[@]}"; do
  filename=$(basename "$file")
  echo "BookmarkBegin" >> "$bookmarks"
  echo "BookmarkTitle: $filename" >> "$bookmarks"
  echo "BookmarkLevel: 1" >> "$bookmarks"
  echo "BookmarkPageNumber: $page" >> "$bookmarks"
  page=$((page + 1))
done

# Final PDF with bookmarks
pdftk "$combined_pdf" update_info "$bookmarks" output "$output_pdf"

# Cleanup
rm "$combined_pdf" "$bookmarks"

echo "Done! Output saved to $output_pdf"
