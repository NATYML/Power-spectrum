#rm ./pk_fit/pk*
make clean
make
#./fit.out ../../output_files/PS_TSC_1e11.dat > pk
#./fit.out ./PS_TSC_1e11.dat 
./fit.out  ./pk/PS_0.50_kf.dat

