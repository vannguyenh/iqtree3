#!/bin/bash
expected_column="$1"

WD="test_scripts/test_data"
input_file="${WD}/expected_memory.tsv"
reported_file="${WD}/time_log.tsv"
selected_columns_file="${WD}/selected_columns.tsv"
cut -f1,2,3,"$expected_column" "$input_file" > "${selected_columns_file}"
final_file="${WD}/combined_with_reported.tsv"

# assuming the reported file and expected file have the same order of commands
cut -f"$expected_column" "$reported_file" > /tmp/reported_column.tsv
paste "$selected_columns_file" /tmp/reported_column.txt > "$final_file"


fail_count=0

echo -e "Command\tExpected\tThreshold\tActual\tExceededBy"

# Skip header
tail -n +2 "$input_file" | while IFS=$'\t' read -r command threshold expected reported; do
    # Compute allowed = expected + threshold
    allowed=$(echo "$expected + $threshold" | bc -l)
    # Check if reported > allowed
    is_exceed=$(echo "$reported > $allowed" | bc -l)
    if [ "$is_exceed" -eq 1 ]; then
        diff=$(echo "$reported - $expected" | bc -l)
        echo -e "$command\t$expected\t$threshold\t$reported\t$diff"
        ((fail_count++))
    fi
done

echo
if [ "$fail_count" -eq 0 ]; then
    echo "✅ All runtime checks passed."
else
    echo "❌ $fail_count checks failed."
    exit 1
fi
