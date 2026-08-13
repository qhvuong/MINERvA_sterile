

# htgettoken -a htvaultprod.fnal.gov -i minerva
htgettoken -i minerva --vaultserver htvaultprod.fnal.gov
export BEARER_TOKEN_FILE=/run/user/`id -u`/bt_u`id -u`
#export BEARER_TOKEN_FILE=/tmp/bt_token_minerva_Analysis_`id -u`
