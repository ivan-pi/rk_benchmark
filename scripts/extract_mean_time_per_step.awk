#!/usr/bin/awk -f

# Extract approach label and mean "us/step" from rk_benchmark text output.
#
# Usage from repository root:
#   ./build/rk_benchmark | ./scripts/extract_mean_time_per_step.awk
#
# Output format:
#   <approach label>\t<mean us/step>
# This tab-separated output is consumed by gnuplot using:
#   using 2:xtic(1)
#
# IMPORTANT:
# This parser is intentionally tied to the fixed-width table emitted by
# `rk_benchmark.f90` in the main benchmark-driver summary. If that table layout
# changes, the column offsets below must be updated too.
#
# Match only the six benchmark data rows:
# - ^[[:space:]]*    : optional leading spaces
# - [1-6]\.          : row index "1." through "6."
# - length($0) >= 80 : guard to ensure fixed-width columns are present
#   (this excludes shorter "Penalty" rows from the overhead-summary section)
/^[[:space:]]*[1-6]\./ && length($0) >= 80 {
  # The label starts at column 5 and spans 30 chars in the benchmark table.
  label = substr($0, 5, 30)
  # Trim leading/trailing whitespace from the fixed-width label field.
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", label)
  # The mean us/step value starts at column 69 and spans 12 chars.
  us_field = substr($0, 69, 12)
  # Require a numeric us/step field; skip anything that does not match.
  # This is a second safeguard against parsing non-summary rows.
  if (us_field !~ /^[[:space:]]*[0-9]+(\.[0-9]+)?([Ee][+-]?[0-9]+)?[[:space:]]*$/) {
    next
  }
  # `+ 0` coerces the extracted numeric substring to numeric form.
  us = us_field + 0
  print label "\t" us
}
