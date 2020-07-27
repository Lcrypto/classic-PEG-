# Compile
cd ../classic_PEG
make all
mv ./MainPEG ../reference_result_linux/MainPEG
make clean
cd ../reference_result_linux

# Run
./MainPEG -numM 252 -numN 504 -codeName 20.alist -degFileName DenEvl_20.deg