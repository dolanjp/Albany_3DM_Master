
grep "Test   #" nightly_log_blakeAlbanySerial.txt >& results0
grep "Test  #" nightly_log_blakeAlbanySerial.txt >& results1
cat results0 results1 >& results11
grep "Test #" nightly_log_blakeAlbanySerial.txt >& results0
cat results11 results0 >& results1
echo "" >> results1
grep " tests failed" nightly_log_blakeAlbanySerial.txt >& results2
cat results1 results2 >& results3
grep "Total Test" nightly_log_blakeAlbanySerial.txt >& results4
cat results3 results4 >& results5
grep "(Failed)" nightly_log_blakeAlbanySerial.txt >& results6
cat results5 results6 >& results_serial
echo "" >> results_serial 
echo "The Albany CDash site can be accessed here: https://my.cdash.org/index.php?project=Albany" >> results_serial
echo "" >> results_serial
rm results0 results1 results11 results2 results3 results4
bash send_email_serial.sh
