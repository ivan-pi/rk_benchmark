/^[[:space:]]*[1-6]\./ && length($0) >= 80 {
  label = substr($0, 5, 30)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", label)
  us = substr($0, 69, 12) + 0
  print label "\t" us
}
