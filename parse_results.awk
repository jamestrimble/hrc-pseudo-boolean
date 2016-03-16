/* NAME/ {
    varnames[$3] = $4;
}

/^v/ {
    for (i=2; i<=NF; i++) {
        if (!match($i, "-")) {
            gsub("x", "", $i);
            print varnames[$i];
        }
    }
}
