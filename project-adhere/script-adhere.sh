#!/bin/bash

sed -e '/W        W/d' -e '/ION  /d' adhesion-tests/membrane-start.gro > adhesion-tests/temp.gro
newcount=$(($(sed -n '$=' adhesion-tests/temp.gro)-3))
echo "membrane without water" > adhesion-tests/membrane-start-no-solvent.gro
echo $newcount >> adhesion-tests/membrane-start-no-solvent.gro
sed -e '1d' -e '2d' adhesion-tests/temp.gro >> adhesion-tests/membrane-start-no-solvent.gro
rm adhesion-tests/temp.gro
