#! /usr/bin/env bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --configfile)
    CONFIGFILE="$2"
    shift # past argument
    shift # past value
    ;;
    --snakefile)
    SNAKEFILE="$2"
    shift # past argument
    shift # past value
    ;;
    --database)
    DATABASE="$2"
    shift # past argument
    shift # past value
    ;;
    --species)
    SPECIES="$2"
    shift # past argument
    shift # past value
    ;;
    --config)
    CONFIG="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SHADOW_DIR="$( paste -d'/' <(echo "$DIR") <(paste -d'_' <( echo ".$SPECIES" ) <( echo "$DATABASE" )) )"


if [ -z "$DATABASE" ]
then
    echo "You need to provide the database name: refseq, ensembl or mytranscripts"
    exit 1
fi

if [ -z "$SPECIES" ]
then
    echo "You need to provide the species name, ie: hsapiens, mmusculus..."
    exit 1
fi

if [ -z "$SNAKEFILE" ]
then
    SNAKEFILE="config/${DATABASE}/${SPECIES}/Snakefile.smk"
fi

if [ -z "$CONFIGFILE" ]
then
    CONFIGFILE="config/${DATABASE}/${SPECIES}/config.yaml"
fi

if [ ! -f $SNAKEFILE ] || [ ! -f $CONFIGFILE ]
then
    echo "The snakefile or configfile requested do not exist. Please, make sure that you are using a supported species."
    exit 1
fi

exec snakemake -pr --use-conda --conda-frontend conda --snakefile $SNAKEFILE --configfile $CONFIGFILE --config db=$DATABASE $CONFIG --directory $DIR --cores 8 all --rerun-incomplete #--unlock #-n
