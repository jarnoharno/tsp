#/bin/bash
for f in input/*; do ./tsp $f > out/`basename $f`; done
