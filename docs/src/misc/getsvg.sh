
#!/bin/bash

# Set the path where PDF files are located
input_path=$(pwd)"/tikz"

# Check if the path exists
if [ ! -d "$input_path" ]; then
  echo "Error: Path does not exist."
  exit 1
fi

# Loop through each PDF file in the specified path
for pdf_file in "$input_path"/*.pdf; do
  if [ -e "$pdf_file" ]; then
    # Get the filename without extension
    filename=$(basename -- "$pdf_file")
    filename_no_ext="${filename%.*}"
    path="$input_path/$filename_no_ext"

    # Convert PDF to SVG 
    inkscape --without-gui --export-plain-svg="$path"".svg" "$path"".pdf"
    # Invert colors
    # inkscape --without-gui --export-type=svg --export-filename="$path""-dark.svg" --batch-process --actions="EditSelectAll;EditInvert;FileSave;FileQuit" "$path"".svg"
    cp "$path"".svg" "$path""-dark.svg"
    sed -i 's/fill:#000000/fill:#ffffff/g' "$path""-dark.svg"
    sed -i 's/stroke:#000000/stroke:#ffffff/g' "$path""-dark.svg"

  fi
done

echo "Conversion complete."
