### Prepare the transcriptome indices for Salmon's mapping-based mode of transcript quantification
Salmon's documentation can be found here: https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode.
Specific information on how to build high-quality salmon indeces for the mapping-based mode can be found here: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/.

To use this mode, we need to build a salmon index for the transcriptome of our interest (the human hg38 and mouse mm39 transcriptome in our case). To do this, we first need to get the transcript and genome sequences:
```
# hg38 sequences
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
mv gencode.v41.transcripts.fa.gz resources/hg38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa.gz resources/hg38

# mm39 sequences
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.transcripts.fa.gz
mv gencode.vM31.transcripts.fa.gz resources/mm39
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz
mv GRCm39.primary_assembly.genome.fa.gz resources/mm39
```

Salmon indexing requires the names of the genome targets and the concatenated transcriptome and genome reference

#### Building the human (hg38) index

```
# Genome targets
grep "^>" <(gunzip -c resources/hg38/GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > resources/hg38/decoys.txt
sed -i.bak -e 's/>//g' resources/hg38/decoys.txt

# Concatenating the transcriptome and genome reference
cat resources/hg38/gencode.v41.transcripts.fa.gz resources/hg38/GRCh38.primary_assembly.genome.fa.gz > resources/hg38/gentrome.fa.gz
```

Now we can build the salmon index:
```
conda activate salmon2
salmon index -t resources/hg38/gentrome.fa.gz -d resources/hg38/decoys.txt -p 12 -i resources/hg38/salmon_index --gencode -k 25
````

### Integrate the human genome (hg38) + virus sequences for related project (case examples):

````
grep "^>" <(cat hg38_with_viruses.fa) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat transcripts_with_viruses.fa hg38_with_viruses.fa > gentrome.fa
gzip gentrome.fa

conda activate salmon2
salmon index -t gentrome.fa.gz -d decoys.txt -p 40 -i salmon_hg38_virus --gencode -k 25
````
