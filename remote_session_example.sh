echo "Running NLSPMC 2D ELDOR simulations" > ~/logfile.txt
echo "start at" >> ~/logfile.txt
date >> ~/logfile.txt
echo "read 95Gtest2.run" | ./nlspmc
echo "finish at" >> ~/logfile.txt
date >> ~/logfile.txt
cat ~/logfile.txt | mail -s "job done" sc974@cornell.edu

