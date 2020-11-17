mkdir -p E13_HL
cd E13_HL
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227225/suppl/GSM4227225%5FE13matrix%2Emtx%2Egz -O matrix.mtx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227225/suppl/GSM4227225%5FE13genes%2Etsv%2Egz -O features.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227225/suppl/GSM4227225%5FE13barcodes%2Etsv%2Egz -O barcodes.tsv.gz

mkdir -p ../E15_HL
cd ../E15_HL
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227226/suppl/GSM4227226%5FE15barcodes%2Etsv%2Egz -O barcodes.tsv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227226/suppl/GSM4227226%5FE15genes%2Etsv%2Egz -O features.tsv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227226/suppl/GSM4227226%5FE15matrix%2Emtx%2Egz -O matrix.mtx.gz

mkdir -p ../E11_HL
cd ../E11_HL
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227224/suppl/GSM4227224%5FE11barcodes%2Etsv%2Egz -O barcodes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227224/suppl/GSM4227224%5FE11genes%2Etsv%2Egz -O features.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4227nnn/GSM4227224/suppl/GSM4227224%5FE11matrix%2Emtx%2Egz -O matrix.mtx.gz
