# Machine-specific data links

Replay and simulation data do not belong in Git. Their locations differ between
the local Mac/T7 setup, cdaq, and ifarm, so each checkout should maintain its
own ignored `.data-links.local.tsv` manifest.

## Configure a machine

From the repository root:

```bash
cp data-links.example.tsv .data-links.local.tsv
```

Edit the second, tab-separated column for that machine, then run:

```bash
tools/manage_data_links.sh apply
tools/manage_data_links.sh check
```

The script creates missing links but will not replace an existing link unless
you explicitly use:

```bash
tools/manage_data_links.sh apply --force
```

It never replaces a real file or directory.

## Recommended workflow

- Keep code and lightweight selected results in Git.
- Keep ROOT files, SIMC outputs, replay reports, logs, and compiled executables
  outside Git.
- Use the same repository branch on the local machine, cdaq, and ifarm.
- Give each machine its own `.data-links.local.tsv`; never commit that file.
- Run `tools/manage_data_links.sh check` after cloning, pulling, or remounting
  external storage.

The repository-relative names used by the analysis macros remain unchanged.
Only their machine-specific targets differ.
