#!/bin/sh
testing_manifest='testing/manifest_'
error_manifest='testing/manifest_errors'

#no errors
echo **"No issues**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_two.tsv testing/samples_two.tsv
cat 'testing/manifest_clean.txt'
echo

#sample duplicates
echo
echo "**Sample col has duplicates**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_two.tsv testing/samples_sampledups.tsv
cat $error_manifest*
echo

#barcode duplicates
echo
echo **"Barcode col has duplicates**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_one.tsv testing/samples_barcodedups.tsv
cat $error_manifest*
echo

#non-alpha numeric
echo
echo "**Barcodes and adaptors have non-alpha numeric values**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_one.tsv testing/samples_nonalpha.tsv
cat $error_manifest*
echo

#filename duplicates
echo
echo "**File_name col has duplicates**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_filedups.tsv testing/samples_two.tsv
cat $error_manifest*
echo

#filename extension incorrect
echo
echo "**File_name col has non fast.gz extensions**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_ext.tsv testing/samples_one.tsv
cat $error_manifest*
echo

#multiplex names don't match
echo
echo **"Multiplex col doesn't match between multiplex and sample files**"
python workflow/scripts/manifest_check.py $testing_manifest testing/multiplex_extrasamples.tsv testing/samples_extrasamples.tsv
cat $error_manifest*
echo
echo

rm $error_manifest*
rm 'testing/manifest_clean.txt'