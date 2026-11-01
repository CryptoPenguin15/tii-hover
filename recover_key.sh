#!/usr/bin/env bash

if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 <challenge-number>"
    exit 1
fi

PARAM="$1"

echo get
curl -fL -o "data/pk_McEliece_${PARAM}.txt" \
  "https://raw.githubusercontent.com/ElenaKirshanova/tii_decoding_challenge/main/public_keyRec/pk_McEliece_${PARAM}.txt"

echo
echo retrieve
time python hover.py "data/pk_McEliece_${PARAM}.txt" --degree auto --checkpoint-dir checkpoints --output "results/my_${PARAM}.json"

echo
echo verify
time python tii_cpu.py verify "data/pk_McEliece_${PARAM}.txt" "results/my_${PARAM}.json"
