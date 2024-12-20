
<<comment

cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/

chmod +x Prepare_supplementary_figures.sh

./Prepare_supplementary_figures.sh

comment
#!/bin/bash

mkdir -p annotated

declare -a imagefiles
imagefiles=(Fig_S1.png Fig_S2.png Fig_S3.png Fig_S4.png)

for file in "${imagefiles[@]}"; do
  convert "$file" -gravity north -pointsize 70 -annotate +0+10 "$(basename "$file")" "annotated/$file"
done

# Combine annotated PNGs into a PDF
convert annotated/* Fig_S1_S4_tmp.pdf

# Create the bookmarks.txt file
rm -f bookmarks.txt
page=1


for file in "${imagefiles[@]}"; do
  echo "BookmarkBegin" >> bookmarks.txt
  echo "BookmarkTitle: $file" >> bookmarks.txt
  echo "BookmarkLevel: 1" >> bookmarks.txt
  echo "BookmarkPageNumber: $page" >> bookmarks.txt
  page=$((page + 1))
done

pdftk Fig_S1_S4_tmp.pdf update_info bookmarks.txt output Fig_S1_S4.pdf

rm  Fig_S1_S4_tmp.pdf
rm bookmarks.txt