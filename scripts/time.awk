#!/usr/bin/awk -f

$1 ~ /^[0-9]+/ {
    total += $NF;
}

END {
    printf "Total time: %8.1f sec\n", total;
}
