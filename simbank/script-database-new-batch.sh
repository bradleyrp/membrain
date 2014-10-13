#!/usr/bin/bash

echo -e "y\ny\ny\ny\n" | ./script-mesosims-database-setup.py
./script-mesosims-database-update.py -c "structure,update_kappas" -b 1200 -e 2200 -n 2017
./script-mesosims-database-update.py -c update_kappas -n 2017
./script-coupling-adv-batch.py -n 2017 --calcs batchcalc
./script-mesosims-database-update.py -c update_datarefs
cd ../
./script-coupling-adv-batch.py -n 2017 -c batchcalc
cd simbank
./script-coupling-adv-batch.py -c megaplot
