# Changelog

## 0.4.1 - 2026-06-09

- expand the documentation with quickstart, input-file, validation, and migration pages
- add a simple Makie-based vibrational spectrum example to the docs
- clarify CLI usage, common failure cases, and the `0.4` structure-input transition

## 0.4.0 - 2026-06-09

- add `read_structure(...)` as the primary public structure reader for restart and XYZ inputs
- align the CLI and documentation on neutral `structure` input naming instead of restart-specific terminology
- define XYZ plus moldescriptor compatibility explicitly: only single-molecule moldescriptor files are accepted because XYZ files do not encode molecule-type assignments
