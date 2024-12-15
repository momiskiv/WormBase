#!/usr/bin/bash
#SBATCH --job-name=hmmsearch_nm471
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=8gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nm471@student.le.ac.uk
#SBATCH --export=NONE

# setting files/directories
out=/home/n/nm471/
hmm=/home/n/nm471/PFam_HMM
dbase=/home/n/nm471/wormbase

# executable from ALICE
hmmsearch=/cm/shared/spack/opt/spack/linux-rocky9-x86_64_v3/gcc-12.3.0/hmmer-3.3.2-6xbbubuci3m4xefi3yqu6thg5loc4jju/bin/hmmsearch

# module needed for using HPC installed software
module load gcc/12.3.0-yxgv2bl
module load openmpi/4.1.5-fzc7xdf
module load hmmer/3.3.2-6xbbubu
s
# Execute the job code
$hmmsearch --tblout $out/PF00890_ancylostoma_caninum.out -E 0.1 --noali $hmm/PF00890.txt $dbase/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF00890_allodiplogaster_sudhausi.out -E 0.1 --noali $hmm/PF00890.txt $dbase/allodiplogaster_sudhausi.PRJEB48369.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF00890_oesophagostomum_dentatum.out -E 0.1 --noali $hmm/PF00890.txt $dbase/oesophagostomum_dentatum.PRJNA72579.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF00171_ancylostoma_caninum.out -E 0.1 --noali $hmm/PF00171.txt $dbase/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF00171_allodiplogaster_sudhausi.out -E 0.1 --noali $hmm/PF00171.txt $dbase/allodiplogaster_sudhausi.PRJEB48369.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF00171_oesophagostomum_dentatum.out -E 0.1 --noali $hmm/PF00171.txt $dbase/oesophagostomum_dentatum.PRJNA72579.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF01127_ancylostoma_caninum.out -E 0.1 --noali $hmm/PF01127.txt $dbase/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF01127_allodiplogaster_sudhausi.out -E 0.1 --noali $hmm/PF01127.txt $dbase/allodiplogaster_sudhausi.PRJEB48369.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF01127_oesophagostomum_dentatum.out -E 0.1 --noali $hmm/PF01127.txt $dbase/oesophagostomum_dentatum.PRJNA72579.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF03937_ancylostoma_caninum.out -E 0.1 --noali $hmm/PF03937.txt $dbase/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF03937_allodiplogaster_sudhausi.out -E 0.1 --noali $hmm/PF03937.txt $dbase/allodiplogaster_sudhausi.PRJEB48369.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF03937_oesophagostomum_dentatum.out -E 0.1 --noali $hmm/PF03937.txt $dbase/oesophagostomum_dentatum.PRJNA72579.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF14290_ancylostoma_caninum.out -E 0.1 --noali $hmm/PF14290.txt $dbase/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF14290_allodiplogaster_sudhausi.out -E 0.1 --noali $hmm/PF14290.txt $dbase/allodiplogaster_sudhausi.PRJEB48369.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF14290_oesophagostomum_dentatum.out -E 0.1 --noali $hmm/PF14290.txt $dbase/oesophagostomum_dentatum.PRJNA72579.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF05328_ancylostoma_caninum.out -E 0.1 --noali $hmm/PF05328.txt $dbase/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF05328_allodiplogaster_sudhausi.out -E 0.1 --noali $hmm/PF05328.txt $dbase/allodiplogaster_sudhausi.PRJEB48369.WBPS19.protein.fa.gz
$hmmsearch --tblout $out/PF05328_oesophagostomum_dentatum.out -E 0.1 --noali $hmm/PF05328.txt $dbase/oesophagostomum_dentatum.PRJNA72579.WBPS19.protein.fa.gz
