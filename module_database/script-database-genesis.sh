#!/usr/bin/bash

'''Run this script only once to prepare a database for module_database.'''

#---parameters
centerdir=/home/rpb/worker/project-postgresql-central

#---start a new server
mkdir /home/rpb/worker/project-postgresql-central
initdb -D $centerdir
cd $centerdir
pg_ctl -D $centerdir -l logfile start

#---create a simbank to interface with membrain codes
createdb membrain_simbank
