# Source this file to get a working CIAO + partial-CALDB environment:
#     source scripts/setup_ciao.sh
# (CIAO 4.18 was installed with conda without the 2.3 GB caldb_main package;
#  the calibration files needed for each ObsID were fetched with
#  download_obsid_caldb into $ASCDS_INSTALL/CALDB instead.)
eval "$(conda shell.bash hook 2>/dev/null || conda shell.zsh hook)"
conda activate ciao-4.18
export CALDB=$ASCDS_INSTALL/CALDB
export CALDBCONFIG=$CALDB/software/tools/caldb.config
export CALDBALIAS=$CALDB/software/tools/alias_config.fits
