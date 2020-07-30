# Compile
cd ../classic_PEG
make all
mv ./MainPEG ../reference_result_linux/MainPEG
make clean
cd ../reference_result_linux

# Run
./MainPEG -numM 256 -numN 512 -codeName 22.alist -degFileName DenEvl_20.deg