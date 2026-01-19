#!/bin/bash

URL="https://datahub.assets.cbioportal.org/msk_chord_2024.tar.gz"

if [ -e "../cbioportal" ]; then
    echo "Error: Destination already exists." >&2
    echo "Please remove or rename it before running this script." >&2
    exit 1
fi

cd ..
curl -L -o "msk_chord_2024.tar.gz" "$URL"

if [ ! -s "msk_chord_2024.tar.gz" ]; then
    echo "Download failed or file is empty."
    exit 1
fi

tar -xzf "msk_chord_2024.tar.gz"
rm -f "msk_chord_2024.tar.gz"
mv msk_chord_2024 cbioportal
cd cbioportal
echo "Download and extraction complete. Files are in $(pwd)"
