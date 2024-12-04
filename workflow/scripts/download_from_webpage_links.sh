#
##########################################################################################
## Download all of the raw data
##########################################################################################
## setting the output directory
#outdir="results/raw/"
#
## links are coming from: 
## https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/raw_datasets/training_dataset/
#wget --no-check-certificate \
#     --no-directories \
#     --recursive \
#     --no-parent \
#     --directory-prefix $outdir \
#     https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/raw_datasets/training_dataset/
#
#
##########################################################################################
## Download all of the batch corrected data
##########################################################################################
## setting the output directory
#outdir="results/batch_corrected/"
#
## links are coming from: 
## https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_allData_batchCorrected_TSV/
#wget --no-check-certificate \
#     --no-directories \
#     --recursive \
#     --no-parent \
#     --directory-prefix $outdir \
#     https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_allData_batchCorrected_TSV/


#########################################################################################
# Download all of the batch corrected data
#########################################################################################
# setting the output directory
outdir="results/harmonized_new/"

# links are coming from: 
# https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/
wget --no-check-certificate \
     --no-directories \
     --recursive \
     --no-parent \
     --directory-prefix $outdir \
     https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/
