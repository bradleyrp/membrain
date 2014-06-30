#!/usr/bin/bash

#---parameters
centerdir=/home/rpb/worker/project-postgresql-central

#---start a new server
mkdir /home/rpb/worker/project-postgresql-central
initdb -D $centerdir
cd $centerdir
pg_ctl -D $centerdir -l logfile start

#---create a simbank to interface with membrain codes
createdb membrain_simbank
