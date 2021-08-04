#!/bin/bash
#
#TOUGHBUGS - written by Brad Hart Nov 2020
#This pipeline is intended to simplify the interrogation of short read WGS while utilising the latest AMR databases from NCBI and CARD.
#required software - ncbi amrfinder, card rgi, prokka, ariba, abricate
#
## TO DO - ADD BACMET
#
toughbugs="$( cd -P "$( dirname "$0" )" && pwd )"
samps=
r4check=(*.fna)
r0check=(*.fasta)
r1check=(*1.fastq.gz)
r2check=(*2.fastq.gz)
work=TOUGHBUGS
cpus=$(nproc)
parp=( 1.Acinetobacter_baumannii 2.Campylobacter 3.Enterococcus_faecalis 4.Enterococcus_faecium 5.Escherichia 6.Klebsiella 7.Salmonella 8.Staphylococcus_aureus 9.Staphylococcus_pseudintermedius 10.Streptococcus_agalactiae 11.Streptococcus_pneumoniae 12.Streptococcus_pyogenes 13.Vibrio_cholerae )
bases=( card megares argannot plasmidfinder resfinder srst2_argannot vfdb_core vfdb_full virulencefinder )
#
echo -e "\nAre the sequences any of the following organisms?:\n"
echo -e "1.Acinetobacter_baumannii\n2.Campylobacter\n3.Enterococcus_faecalis\n4.Enterococcus_faecium\n5.Escherichia\n6.Klebsiella\n7.Salmonella\n8.Staphylococcus_aureus\n9.Staphylococcus_pseudintermedius\n10.Streptococcus_agalactiae\n11.Streptococcus_pneumoniae\n12.Streptococcus_pyogenes\n13.Vibrio_cholerae\n"
echo
read -p "Yes/No? (y/n): " yehnah
if [[ "$yehnah" =~ ^([yY][eE][sS]|[yY])+$ ]]; then
echo
read -p "Input number matching your organism: " choice
toot=$(($choice - 1))
organism=$(echo ${parp[$toot]} | sed 's/.*\.//g')
fi
#
echo -e "Samples\t""Read1\t""Read2\t""Description" > samples_metadata.txt
if test -n "$(shopt -s nullglob; echo ${r0check})"; then
for i in ${r0check[*]}; do
echo -e "${i%.fasta}""\t""GENOME""\t""${i}""\t""Reads from ${i%.fasta}" >> samples_metadata.txt
done
echo -e "\nfasta files added to samples_metadata.txt"
else
echo -e "\nno fasta files detected. Ignore if intentional."
fi
if test -n "$(shopt -s nullglob; echo ${r4check})"; then
for i in ${r4check[*]}; do
echo -e "${i%.fna}""\t""GENOME""\t""${i}""\t""Reads from ${i%.fna}" >> samples_metadata.txt
done
echo -e "\nfna files added to samples_metadata.txt"
else
echo -e "\nno fna files detected. Ignore if intentional."
fi
if [[ -f "$r1check" && -f "$r2check" ]] ; then
for i in $(ls *1.fastq.gz); do
  echo -e "$(echo $i | sed "s/_.*//")""\t""$i""\t""$(echo $i | sed "s/1.fastq.gz/2.fastq.gz/")""\t""Reads from $(echo $i | sed "s/_.*//")" >> samples_metadata.txt;
done
echo -e "\nfastq.gz files added to samples_metadata.txt"
else
echo -e "\nMissing or incomplete fastq.gz file set with naming syntax of xxxx1.fastq.gz and xxxx2.fastq.gz"
echo "Ignore this message if this is intentional."
fi
lines=$(cat samples_metadata.txt | tail -n +2 | wc -l)
if [ "$lines" -gt 0 ]; then
  echo -e "\nsamples_metadata.txt file created, type: cat samples_metadata.txt to check entries"
  else
  echo -e "\nNO SEQUENCES DETECTED - ENSURE SEQUENCES ARE PRESENT IN FOLDER YOU RAN THIS SCRIPT FROM\n"
  rm -rf samples_metadata.txt
  exit 1
fi
#
samps=$(cat samples_metadata.txt | tail -n +2 | awk '{print $1}' | sort -u)
#
## CREATE WORKING DIRECTORY AND VARIABLE
mkdir -p $work
if [ ! -d $work ]; then
        echo -e "\nERROR: $OUTPUT could not be created in the selected location. Please check\n"
        exit 1
fi
#
for i in $samps; do
  mkdir -p "$work"/"$i"
  if [ ! -d "$work"/"$i" ]; then
        echo -e "\nERROR: $work/$i could not be created in the selected location. Please check\n"
        exit 1
  fi
done
#
for i in $samps; do
  cp ${i}* $work/
done
#
for i in $samps; do
if [ -e ${work}/${i}.faa -a -e ${work}/${i}.gff ]; then
 echo -e "Both ${i}.faa and gff files present, continuing\n"
 else
 echo -e "oh no, no faa or gff file, time for some prokka!\n"
 sleep 1 
 mkdir -p ${work}/${i}/annotation/
  if [[ "$yehnah" =~ ^([yY][eE][sS]|[yY])+$ ]]; then
    echo -e "\n Running prokka on ${i} to create .gff and .faa files"
    prokka $work/$i\.fasta --outdir $work/$i/annotation/ --prefix $i --locustag $i --genus $organism --cpus $cpus --force > /dev/null 2>&1
    cp ${work}/${i}/annotation/${i}.{faa,gff} $work/
  else
    echo -e "\n Running prokka on ${i} to create .gff and .faa files"
    prokka $work/$i\.fasta --outdir $work/$i/annotation/ --prefix $i --locustag $i --cpus $cpus --force > /dev/null 2>&1
    cp ${work}/${i}/annotation/${i}.{faa,gff} $work/
  fi
fi
done
## RGI Card resistance checking
echo -e "Beginning RGI Card Resistance Database search\n"
 sleep 1
for i in $samps; do
  mkdir -p $work/$i/card_resistance
 rgi main -i $work/${i}.fasta -o "$work"/"$i"/card_resistance/"$i"_card_fasta -t contig -a BLAST -n $cpus --clean
 rgi main -i $work/${i}.faa -o "$work"/"$i"/card_resistance/"$i"_card_faa -t protein -a DIAMOND -n $cpus --clean
done
#
#Amend gff file format for amrfinderplus
#
for i in $samps; do
  perl -pe '/^##FASTA/ && exit; s/(\W)Name=/$1OldName=/i; s/ID=([^;]+)/ID=$1;Name=$1/' ${work}/${i}.gff > ${work}/${i}amr.gff
done
#
#amrfinderplus run
#
echo -e "\n Running NCBI AMRfinder Plus"
 sleep 1
for i in $samps; do
  echo -e "\nRunning NCBI AMRfinder Plus on ${i}"
  mkdir -p ${work}/${i}/amrfinderplus/
if [[ "$yehnah" =~ ^([yY][eE][sS]|[yY])+$ ]]; then
  amrfinder -p ${work}/${i}.faa -g ${work}/${i}amr.gff -n ${work}/${i}.fasta -O $organism --plus --protein_output ${work}/${i}/amrfinderplus/${i}amrfprot.fa --nucleotide_output ${work}/${i}/amrfinderplus/${i}amrfnuc.fa --mutation_all ${work}/${i}/amrfinderplus/${i}_amrfmutations --name ${i} --output ${work}/${i}/amrfinderplus/${i}_amrf_out.txt --threads $cpus
else
  echo -e "\nRunning NCBI AMRfinder Plus on ${i}"
  amrfinder -p ${work}/${i}.faa -g ${work}/${i}amr.gff -n ${work}/${i}.fasta --plus --protein_output ${work}/${i}/amrfinderplus/${i}amrfprot.fa --nucleotide_output ${work}/${i}/amrfinderplus/${i}amrfnuc.fa --name ${i} --output ${work}/${i}/amrfinderplus/${i}_amrf_out.txt --threads $cpus
fi
done
#
#ARIBA analysis
#
echo -e "\n Running ARIBA analysis\n"
for i in ${samps[*]}; do
echo $i >> samps.txt
done
#
while read line; do
mkdir -p ${work}/${line}/ariba
for i in ${bases[*]}; do
echo -e "\n ${line} -  ${i} Analysis\n"
ariba run --threads ${cpus} /home/harbj019/dbs/ariba1/out.${i}.prepareref ${work}/${line}*1.fastq.gz ${work}/${line}*2.fastq.gz ${work}/${line}/ariba/${i}
done
done < samps.txt

#
#Abricate analysis
#
echo -e "\n Running Abricate Analysis\n"
sleep 1
while read line; do
mkdir -p ${work}/${line}/abricate
db=( resfinder card argannot vfdb plasmidfinder ecoh ncbi ecoli_vf megares )
for i in ${db[@]}; do
abricate --db ${i} --minid 85 --mincov 85 --nopath ${work}/${line}.fasta > ${work}/${line}/abricate/${line}_${i}.tab 
done
done < samps.txt
#
echo -e "\n Abricate Analysis Complete\n"
#
rm -f ${work}/*.{faa,gff,gz,fasta,fa,fasta,fna} > /dev/null 2>&1
for i in ${samps}; do
rm -f ${work}/${i}/annotation > /dev/null 2>&1
done
rm samples_metadata.txt
rm samps.txt
echo -e "\n All Done!\n"
