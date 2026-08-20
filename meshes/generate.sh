#!/usr/bin/env bash
set -euo pipefail

TRIANGLE_BIN="${TRIANGLE_BIN:-triangle}"

A=0.4096;"${TRIANGLE_BIN}" -pqa${A}Djz square.poly
A=0.1024;"${TRIANGLE_BIN}" -prqa${A}Djz square.1
A=0.0512;"${TRIANGLE_BIN}" -prqa${A}Djz square.2
## Agressive
A=0.0016;"${TRIANGLE_BIN}" -prqa${A}Djz square.3
A=0.0001;"${TRIANGLE_BIN}" -prqa${A}Djz square.4
A=0.00000625; "${TRIANGLE_BIN}" -prqa${A}Djz square.5
A=0.000000390625; "${TRIANGLE_BIN}" -pqa${A}Djz square.6


mkdir -p dirt

mv square.[1567].* dirt/
#
cp precious/square.[234].* dirt/
#
rm -f square.[234].*


cp precious/square.msh dirt/


## For mg
cp square.poly squaremg.poly
A=0.00640000;"${TRIANGLE_BIN}" -pqa${A}Djz squaremg.poly
A=0.00160000;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.1
A=0.00040000;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.2
A=0.00010000;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.3
A=0.00002500;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.4
A=0.00000600;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.5
A=0.00000150;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.6
A=0.00000035;"${TRIANGLE_BIN}" -prqa${A}Djz squaremg.7

mv squaremg.[12345678].* dirt/
rm -f squaremg.poly
