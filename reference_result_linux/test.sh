# Compile
cd ../classic_PEG
make all
mv ./MainPEG ../reference_result_linux/MainPEG
make clean
cd ../reference_result_linux

# Run
# It's not performant enough to create large matrices btw, even m=2000 is abit much
# some of the matrices are not proper fit
codeLen=100

## declare an array variable
declare -a arr=("DenEvl_cr_0.50_thr_0.1071" \
"DenEvl_cr_0.55_thr_0.0904" "DenEvl_cr_0.60_thr_0.0766" \
"DenEvl_cr_0.65_thr_0.0633" "DenEvl_cr_0.70_thr_0.0504" \
"DenEvl_cr_0.75_thr_0.0392" "DenEvl_cr_0.80_thr_0.0298" \
"DenEvl_cr_0.85_thr_0.0199" "DenEvl_cr_0.90_thr_0.0109")

declare -a arrNumerators=("50" \
"55" "60" \
"65" "70" \
"75" "80" \
"85" "90")

declare -a arrDenominators=("100" \
"100" "100" \
"100" "100" \
"100" "100" \
"100" "100")

# get length of an array
arraylength=${#arr[@]}

## now loop through the above array
for (( i=1; i<${arraylength}+1; i++ ));
do
    codeRateNumerator=${arrNumerators[$i-1]}
    codeRateDenominator=${arrDenominators[$i-1]}
    msgLen=$(python -c "print(($codeLen * 1.0)/($codeRateDenominator)*($codeRateNumerator))")

    printf "Running $msgLen $codeLen $codeRateNumerator/$codeRateDenominator"

    ./MainPEG -numM $msgLen -numN $codeLen \
    -codeName "${arr[$i-1]}.alist" \
    -degFileName "DenEvls/arXiv:0901.2140/${arr[$i-1]}.deg"
done