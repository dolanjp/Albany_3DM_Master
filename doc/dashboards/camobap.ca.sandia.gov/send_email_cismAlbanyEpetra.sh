#!/bin/bash

#source $1 

TTT=`grep "(Failed)" /home/ikalash/nightlyCDash/nightly_log_cismAlbanyEpetra.txt -c`
TTTT=`grep "(Not Run)" /home/ikalash/nightlyCDash/nightly_log_cismAlbanyEpetra.txt -c`
TTTTT=`grep "(Timeout)" /home/ikalash/nightlyCDash/nightly_log_cismAlbanyEpetra.txt -c`
TT=`grep "...   Passed" /home/ikalash/nightlyCDash/nightly_log_cismAlbanyEpetra.txt -c`


#/bin/mail -s "Albany ($ALBANY_BRANCH): $TTT" "albany-regression@software.sandia.gov" < $ALBOUTDIR/albany_runtests.out
/bin/mail -s "IKTCismAlbanyEpetra, camobap.ca.sandia.gov: $TT tests passed, $TTT tests failed, $TTTT tests not run, $TTTTT timeouts" "ikalash@sandia.gov, lbertag@sandia.gov, mperego@sandia.gov" < /home/ikalash/nightlyCDash/results_cismAlbanyEpetra
