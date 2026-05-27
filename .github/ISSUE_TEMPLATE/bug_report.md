---
name: Bug report
about: Report a defect in the pipeline
title: '[bug] '
labels: bug
assignees: ''
---

## Summary
<!-- One-line description of the problem. -->

## Steps to reproduce
1. `nextflow run main.nf -profile ... --sample_sheet ...`
2. ...

## Expected behavior
<!-- What you expected to happen. -->

## Actual behavior
<!-- What actually happened. Include the failing process name and last few lines of the .command.log. -->

## Environment
- Pipeline version (from `nextflow run main.nf --version`):
- Nextflow version (`nextflow -version`):
- Profile(s) used (e.g. `conda,docker,gandalf`):
- Host/OS (e.g. gandalf RHEL 9):
- Sample sheet (redact PHI):

## Logs and outputs
<!--
Attach (or paste in fenced blocks):
- .nextflow.log (or relevant tail)
- the failing task's .command.log and .command.err
- pipeline_info/trace.txt if available

Do NOT paste real patient identifiers. Sequencing IDs only.
-->

## Workaround tried
<!-- If you tried -resume, --skip_*, etc., note what worked or didn't. -->
