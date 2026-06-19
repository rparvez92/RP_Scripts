#!/usr/bin/env bash

set -u

repo_root=$(git rev-parse --show-toplevel 2>/dev/null) || {
  echo "Run this script from inside the RP_Scripts repository." >&2
  exit 1
}

config=${DATA_LINKS_CONFIG:-"$repo_root/.data-links.local.tsv"}
mode=${1:-check}
force=${2:-}

usage() {
  cat <<'EOF'
Usage:
  tools/manage_data_links.sh check
  tools/manage_data_links.sh apply [--force]

The default manifest is .data-links.local.tsv. Override it with:
  DATA_LINKS_CONFIG=/path/to/manifest.tsv tools/manage_data_links.sh check

Each non-comment line in the manifest contains:
  repository-relative link path<TAB>absolute target path
EOF
}

if [[ "$mode" != "check" && "$mode" != "apply" ]]; then
  usage >&2
  exit 2
fi

if [[ "$force" != "" && "$force" != "--force" ]]; then
  usage >&2
  exit 2
fi

if [[ ! -f "$config" ]]; then
  echo "Manifest not found: $config" >&2
  echo "Copy data-links.example.tsv to .data-links.local.tsv and edit it for this machine." >&2
  exit 1
fi

checked=0
healthy=0
changed=0
problems=0

while IFS=$'\t' read -r link_path target extra; do
  [[ -z "${link_path// }" || "$link_path" == \#* ]] && continue

  checked=$((checked + 1))

  if [[ -n "${extra:-}" || -z "${target:-}" ]]; then
    echo "INVALID  $link_path (expected exactly two tab-separated fields)" >&2
    problems=$((problems + 1))
    continue
  fi

  if [[ "$link_path" = /* || "$link_path" == *".."* ]]; then
    echo "INVALID  $link_path (link path must stay inside the repository)" >&2
    problems=$((problems + 1))
    continue
  fi

  if [[ "$target" != /* ]]; then
    echo "INVALID  $link_path (target must be an absolute path)" >&2
    problems=$((problems + 1))
    continue
  fi

  destination="$repo_root/$link_path"

  if [[ -L "$destination" ]]; then
    current=$(readlink "$destination")
    if [[ "$current" == "$target" ]]; then
      if [[ -e "$target" ]]; then
        echo "OK       $link_path -> $target"
        healthy=$((healthy + 1))
      else
        echo "BROKEN   $link_path -> $target" >&2
        problems=$((problems + 1))
      fi
      continue
    fi

    if [[ "$mode" == "apply" && "$force" == "--force" ]]; then
      unlink "$destination"
    else
      echo "MISMATCH $link_path -> $current"
      echo "         expected -> $target"
      problems=$((problems + 1))
      continue
    fi
  elif [[ -e "$destination" ]]; then
    echo "BLOCKED  $link_path exists and is not a symlink" >&2
    problems=$((problems + 1))
    continue
  fi

  if [[ "$mode" == "check" ]]; then
    echo "MISSING  $link_path -> $target"
    problems=$((problems + 1))
    continue
  fi

  if [[ ! -e "$target" ]]; then
    echo "WARNING  target does not currently exist: $target" >&2
  fi

  mkdir -p "$(dirname "$destination")"
  ln -s "$target" "$destination"
  echo "CREATED  $link_path -> $target"
  changed=$((changed + 1))
done < "$config"

echo
echo "Checked: $checked; healthy: $healthy; changed: $changed; problems: $problems"

(( problems == 0 ))
