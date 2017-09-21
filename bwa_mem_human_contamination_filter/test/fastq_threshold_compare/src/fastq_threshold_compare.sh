#!/bin/bash
# fastq_threshold_compare 0.0.1


# Absolute value
 abs() {
   (( $(bc <<< "scale=7; $1 < 0") )) && bc <<< "scale=7; $1 * -1" || echo $1
 }

# ratio max check. $1-test $2-max
check-max-ratio() {
    per_inp=$(float-to-percent $1)
    per_max=$(float-to-percent $2)
    (( $(bc <<<"scale=7; $1 < $2") )) && echo "PASS: Filtered down to ${per_inp}%, fine" || echo "FAIL: ratio remaining after filter, ${per_inp}%, too high. Expected ${per_max}%"
}

# turn $1 from % to decimal
percent-to-float() {
    bc <<<"scale=7; $1 / 100"
}

# turn $1 from % to decimal
float-to-percent() {
    bc <<<"scale=7; $1 * 100"
}

###############################################################
#
# function verifies line counts between 2 files are equal
#
# function arguments:
#       1. file 1
#       2. file 2
#
# returns string starting with PASS or FAIL sentinel. Assumes
# Both files are gzipped or not.
#
###############################################################
verify-line-counts() {
    if [[ $1 == *\.gz ]] && [[ $2 == *\.gz ]]; then
        param_one_count=$(<"${1}" zcat | wc -l)
        param_two_count=$(<"${2}" zcat | wc -l)
    else
        param_one_count=$(<"${1}" wc -l)
        param_two_count=$(<"${2}" wc -l)        
    fi

    if (( param_one_count == param_two_count )); then
        echo "PASS: $1 and $2 are equal in lenght"
    else
        echo "FAIL: $1 and $2 are not equal in lenght"
    fi
}


###############################################################
#
# function calculates percentage reads removed within threshold
#
# function arguments:
#       1. original fastq
#       2. decontaminated fastq
#
# returns ratio of reads removed. Assumes all inputs are gzipped
#
###############################################################
calculate-reads-removed() {
    if [[ $1 == *\.gz ]] && [[ $2 == *\.gz ]]; then
        ori_count=$(<"${1}" zcat | wc -l)
        decontaminate_count=$(<"${2}" zcat | wc -l)
    else
        ori_count=$(wc -l <"$1")
        decontaminate_count=$(wc -l <"$2")      
    fi

    bc <<<"scale=7; $decontaminate_count/$ori_count"
}


###############################################################
#
# function verifies an input is within a +/- tolerance
#
# function arguments:
#       1. Input int/float
#       2. Expected int/float
#       3. Tolerance int/float
#
# returns a string 
#
###############################################################
threshold-verify() {
    difference=$(bc <<<"scale=7; $1 - $2")
    difference=$(abs $difference)
    thresh=$(abs $3)

    perc_diff=$(float-to-percent $difference)
    perc_exp=$(float-to-percent $2)
    (( $(bc <<< "$difference > $thresh") )) && echo "FAIL: Deviated ${perc_diff}% is outside +/- ${thresh}% window of ${perc_exp}%" || echo "PASS: Deviated ${perc_diff}% is within +/- ${thresh}% window of ${perc_exp}%"
}

main() {

    echo "Value of ori_fastq: '$ori_fastq'"
    echo "Value of ori_fastq2: '$ori_fastq2'"
    echo "Value of new_fastq: '$new_fastq'"
    echo "Value of new_fastq2: '$new_fastq2'"
    echo "Value of expected_filter_rate: '${expected_filter_rate}'%"
    echo "Value of tolerance: '${tolerance}'%"
    echo "Value of max_filter: '${max_filter}'%"

    dx-download-all-inputs
    expected_filter_rate=$(percent-to-float $expected_filter_rate)
    tolerance=$(percent-to-float $tolerance)
    max_filter=$(percent-to-float $max_filter)


    declare -a fastq_file_validation  # Logs validation process

    fastq_file_validation+=("$(fastQValidator --file "${ori_fastq_path}" 1>&2 && echo "PASS: ${ori_fastq_name} valid FASTQ " || echo "FAIL: ${ori_fastq_name} invalid FASTQ")")
    if [ -n "$ori_fastq2" ]
    then
        fastq_file_validation+=("$(fastQValidator --file "${ori_fastq2_path}" 1>&2 && echo "PASS: ${ori_fastq2_name} valid FASTQ" || echo "FAIL: ${ori_fastq2_name} invalid FASTQ")")
        fastq_file_validation+=("$(verify-line-counts "${ori_fastq_path}" "${ori_fastq2_path}")")
    fi

    fastq_file_validation+=("$(fastQValidator --file "${new_fastq_path}" 1>&2 && echo "PASS: ${new_fastq_name} valid FASTQ" || echo "FAIL: ${new_fastq_name} invalid FASTQ")")
    if [ -n "$new_fastq2" ]
    then
        fastq_file_validation+=("$(fastQValidator --file "${new_fastq2_path}" 1>&2 && echo "PASS: ${new_fastq2_name} valid FASTQ" || echo "FAIL: ${new_fastq2_name} invalid FASTQ")")
        fastq_file_validation+=("$(verify-line-counts "${new_fastq_path}" "${new_fastq2_path}")")
        if [ -z "$ori_fastq2" ]; then
            fastq_file_validation+=("FAIL: Must provide both paired end files")
        fi
    fi

    ratio=$(calculate-reads-removed "${ori_fastq_path}" "${new_fastq_path}")
    fastq_file_validation+=("$(check-max-ratio $ratio $max_filter)")
    fastq_file_validation+=("$(threshold-verify "${ratio}" "${expected_filter_rate}" "${tolerance}")")

    # Print results and pass of fail job
    sentinel=0
    for msg in "${fastq_file_validation[@]}"; do
        echo $msg
        if [[ $msg == FAIL* ]]; then
            sentinel=1
        fi
    done
    exit $sentinel
}
