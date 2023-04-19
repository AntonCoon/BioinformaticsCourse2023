OUT_DIR="pollux_output"
INPUT_FILENAME="_1.fastq"
INPUT_FASTQ="input/$INPUT_FILENAME"
OUT_TRIM_FASTQ="$OUT_DIR/$INPUT_FILENAME.trimmed"
OUT_FASTQ="$OUT_DIR/$INPUT_FILENAME.trimmed.corrected"
META="input/meta.txt"

mkdir -p input

echo "Running simulation python script"
python3 simulation.py $INPUT_FASTQ $META
echo "Read results are in $INPUT_FASTQ"
echo

# If default parameters are used then all sequences are recognized as wrong
# So parameters are tuned
echo "Running trimmomatic"
java -jar lib/trimmomatic-0.39.jar SE -phred33 $INPUT_FASTQ $OUT_TRIM_FASTQ \
  LEADING:1 TRAILING:1 SLIDINGWINDOW:5:2 MINLEN:36
echo

echo "Analyze trimming software"
python3 analyze_trimming.py $INPUT_FASTQ $OUT_TRIM_FASTQ
echo

# !!! Script is copied from
# https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/scripts/wrappers/run.pollux.sh
echo "Running error correcting software"
./run.pollux.sh $OUT_TRIM_FASTQ _ $OUT_DIR 25
echo "Result is in $OUT_FASTQ"
echo

echo "Analyze error correcting software"
python3 analyze_correcting.py $OUT_TRIM_FASTQ $OUT_FASTQ $META
echo
