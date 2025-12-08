#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

patterns=(
  ".snakemake/"
  "logs/"
  "../rmds/*.Rds"
  "../rmds/*.html"
  "../helpers/*.html" # Accidental renderings
  "../rmds/*_cache/"
  "../rmds/*_files/"
  "../rmds/figures/"
  "../rmds/._*"
  "../output/"
)

items_to_remove=()
for pat in "${patterns[@]}"; do
  for entry in $pat; do
    [[ -e $entry ]] || continue
    items_to_remove+=("$entry")
  done
done

if (( ${#items_to_remove[@]} )); then
  echo "Will remove:"
  printf '  %s\n' "${items_to_remove[@]}"
  read -p "Remove ALL listed items? [y/N] " ans
  if [[ $ans =~ ^[Yy]$ ]]; then
    for item in "${items_to_remove[@]}"; do
      rm -rf -- "$item"
      echo "Removed: $item"
    done
    echo "Done."
  else
    echo "Aborted."
  fi
else
  echo "Nothing to remove."
fi
