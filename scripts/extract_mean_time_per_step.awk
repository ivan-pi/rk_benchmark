# Extract approach label and mean "us/step" from rk_benchmark text output.
#
# Usage from repository root:
#   ./build/rk_benchmark | awk -f scripts/extract_mean_time_per_step.awk
#
# Output format:
#   <approach label>\t<mean us/step>
# This tab-separated output is consumed by gnuplot using:
#   using 2:xtic(1)
#
# Match only the six benchmark data rows:
# - ^[[:space:]]*    : optional leading spaces
# - [1-6]\.          : row index "1." through "6."
# - length($0) >= 80 : guard to ensure fixed-width columns are present
/^[[:space:]]*[1-6]\./ && length($0) >= 80 {
  # The label starts at column 5 and spans 30 chars in the benchmark table.
  label = substr($0, 5, 30)
  # Trim leading/trailing whitespace from the fixed-width label field.
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", label)
  # The mean us/step value starts at column 69 and spans 12 chars.
  # `+ 0` coerces the extracted substring to numeric form.
  us = substr($0, 69, 12) + 0
  print label "\t" us
}
