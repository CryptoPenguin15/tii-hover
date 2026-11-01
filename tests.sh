python recover_key.py --demo --method exhaustive --output demo/key.json
python tii_cpu.py verify demo/demo_public_key.txt demo/key.json

python recover_key.py data/pk_McEliece_100.txt --method algebraic --output recovered_100.json --timeout 300

python recover_key.py demo/demo_public_key.txt --method exhaustive --checkpoint demo/search.json --max-candidates 3
python recover_key.py demo/demo_public_key.txt --method exhaustive --checkpoint demo/search.json --resume --output demo/key.json

python tii_cpu.py inspect data/pk_McEliece_100.txt
sage load_tii.sage data/pk_McEliece_100.txt
python tii_cpu.py preflight data/pk_McEliece_100.txt
python -m unittest discover -s tests -v
