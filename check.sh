#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Steven Van Overberghe
#
# Run by "make check".  Compares the numbers of graphs genll and genllL generate with
# the known counts in counts.txt, then checks that the parts of a res/mod split add up
# to the whole.

cd "$(dirname "$0")" || exit 1

for prog in genll genllL; do
    [ -x "./$prog" ] || { echo "./$prog not built -- run make first" >&2; exit 1; }
done

checks=0
failures=0

# The number of graphs a run generates, read from its ">Z" summary line.
count() {
    "$@" -u 2>&1 >/dev/null | sed -n 's/^>Z \([0-9]*\) graphs.*/\1/p'
}

# verify EXPECTED DESCRIPTION ACTUAL
verify() {
    checks=$((checks + 1))
    if [ "$3" != "$1" ]; then
        echo "FAIL  $2: got ${3:-nothing}, expected $1"
        failures=$((failures + 1))
    fi
}

# --- the known counts, with both executables ------------------------------------
for prog in ./genll ./genllL; do
    while read -r n ll mll c4 mc4 c5; do
        case $n in ''|\#*) continue ;; esac
        [ "$ll"  = - ] || verify "$ll"  "$prog $n"       "$(count "$prog" "$n")"
        [ "$mll" = - ] || verify "$mll" "$prog -M $n"    "$(count "$prog" -M "$n")"
        [ "$c4"  = - ] || verify "$c4"  "$prog -f $n"    "$(count "$prog" -f "$n")"
        [ "$mc4" = - ] || verify "$mc4" "$prog -f -M $n" "$(count "$prog" -f -M "$n")"
        [ "$c5"  = - ] || verify "$c5"  "$prog -5 $n"    "$(count "$prog" -5 "$n")"
        [ "$c5"  = - ] || verify "$c5"  "$prog -f -5 $n" "$(count "$prog" -f -5 "$n")"
    done < counts.txt
    echo "$prog: known counts checked"
done

# --- res/mod: the parts of a split add up to the whole ---------------------------
for options in "" "-f" "-M" "-5"; do
    case $options in
        -f) n=22 ;;
        -5) n=23 ;;
        *)  n=18 ;;
    esac
    whole=$(count ./genll $options $n)
    sum=0
    for r in 0 1 2 3 4 5 6; do
        sum=$((sum + $(count ./genll $options "$r/7" $n)))
    done
    verify "$whole" "res/mod: sum of 7 parts of genll $options $n" "$sum"
done
echo "res/mod splits checked"

if [ "$failures" -eq 0 ]; then
    echo "all $checks checks passed"
else
    echo "$failures of $checks checks FAILED"
    exit 1
fi
