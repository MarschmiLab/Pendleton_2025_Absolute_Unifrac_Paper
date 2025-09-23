for file in $(cat accessions.txt)
do 
fastq-dump --split-files --gzip $file
done

