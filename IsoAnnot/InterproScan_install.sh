!/bin/bash
# Title: InterProScan installation script for IsoAnnot
# Author: Darío González
# Description: This file installs InterProScan within IsoAnnot


Variables
INTERPROSCAN="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz"
INTERPROSCAN_HASH="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz.md5"
INTERPROSCAN_DOWNLOAD=$(basename $INTERPROSCAN)  # name of the downloaded file
INTERPROSCAN_HASH_DOWNLOAD=$(basename $INTERPROSCAN_HASH)  # name of the hash file

SIGNALP_DOWNLOAD_PAGE="https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=4.1&packageversion=4.1g&platform=Linux"
TMHMM_DOWNLOAD_PAGE="https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=tmhmm&version=2.0c&packageversion=2.0c&platform=Linux"

# Download InterProScan database
echo "         ###################################################"
echo "         ## InterProScan installation script for IsoAnnot ##"
echo "         ###################################################"
echo 
echo "Please choose one option of the following: "
echo "  a) I want to INSTALL InterProScan from scratch."
echo "  b) I moved the location of the IsoAnnot folder and want to MODIFY InterProScan installation PATHS."
echo
read -p "Enter option a or b: " execution_option


if [ $execution_option == "a" ]; then
        
    echo "# Downloading InterProScan #"

    MAX_RETRIES=3
    TRY=0

    while [ $TRY -lt $MAX_RETRIES ]
    do
        # Download the files
        echo "# Downloading the file #"
        echo "This may take a while"
        wget ${INTERPROSCAN}
        wget ${INTERPROSCAN_HASH}

        # Calculate the downloaded file's MD5 hash
        echo "# Checking integrity of the download #"
        EXPECTED_HASH=$(awk '{print $1}' $INTERPROSCAN_HASH_DOWNLOAD)
        DOWNLOADED_HASH=$(md5sum $INTERPROSCAN_DOWNLOAD | awk '{print $1}')

        # Retry download if needed
        if [ $DOWNLOADED_HASH = $EXPECTED_HASH ]
        then
            echo "File downloaded successfully and verified!"
            break
        else
            echo "The downloaded file is corrupted!"
            echo "Retrying..."
            TRY += 1
            rm -f $INTERPROSCAN_DOWNLOAD
            rm -f $INTERPROSCAN_HASH_DOWNLOAD
        fi
    done

    # Exit if maximum download retries were reached
    if [ $TRY -ge $MAX_RETRIES ]; then
        echo "Failed to download and verify the file after $TRY attempts."
        echo "Aborting."
        exit 1
    fi

    # Extract the files directly to the folder IsoAnnot needs
    echo "# Decompressing the download #"
    tar -pxvzf interproscan-*-bit.tar.gz --one-top-level=software/interproscan --strip-components 1

    # Run the setup file of InterProScan
    echo "## Configuring InterProScan ##"
    cd software/interproscan/
    python3 setup.py -f interproscan.properties
    cd ../..

    # Install the propietary databases in InterProScan
    echo
    echo "###############################################################"
    echo "## INSTALL LICENSED DATABASES                                ##"
    echo "##                                                           ##"
    echo "## IsoAnnot uses 2 licensed databases from InterProScan.     ##"
    echo "## You can download them for free for academic purposes.     ##"
    echo "##                                                           ##"
    echo "## This script helps you automate the installation for each  ##"
    echo "## database, but you will need to register and accept the    ##"
    echo "## licenses on your own.                                     ##"
    echo "##                                                           ##"
    echo "## IF YOU DON'T INSTALL THIS DATABASES ISOANNOT WON'T WORK!! ##"
    echo "###############################################################"
    echo
    read -p "Continue with the installation? [yes/no]: " ANSWER

    # Exit the program if the user requested it
    if [[ $ANSWER =~ ^[Nn][Oo]$ ]]
    then
        echo "Installation cancelled."
        exit 1
    fi

    # Create a temporal folder to store the downloaded files
    mkdir -p software/temp

    # Print SignalP manual installation message
    echo
    echo "# Installing SignalP #"
    echo "Follow the link bellow to the download page for SignalP, read and accept the license."
    echo "You will receive an email with the download link for the software when you accept."
    echo 
    echo "$SIGNALP_DOWNLOAD_PAGE"
    echo 

    ANSWER="No"  # set the variable for input loop

    while [[ $ANSWER =~ ^[Nn][Oo]$ ]]
    do
        read -p "Paste the link from the email here: " LINK
        read -p "Check the link is correct. Do you confirm the link? [yes/no]: " ANSWER
    done

    # Download and decompress the program
    echo "Downloading file"
    wget -P software/temp -nd  -np -r -nH -l1 -A "signalp*.tar.gz" $LINK

    tar -xvzf software/temp/signalp-* -C software/interproscan/bin/signalp/4.1/ --strip-components=1

    # Move files to InterProScan folder
    echo "Moving required files to InterProScan file tree"
    mv software/temp/signalp-*/* software/interproscan/bin/signalp/4.1/

    # Modify the `signalp` binary file so it can execute correctly
    echo "Configuring SignalP"
    PATH_TO_SIGNALP=$(realpath software/interproscan/bin/signalp/4.1)
    echo $PATH_TO_SIGNALP
    mv software/interproscan/bin/signalp/4.1/signalp software/interproscan/bin/signalp/4.1/signalp_backup
    awk -v var="    \$ENV{SIGNALP} = '$PATH_TO_SIGNALP'" '{ if (NR == 13) print var; else print $0}' software/interproscan/bin/signalp/4.1/signalp_backup > software/interproscan/bin/signalp/4.1/signalp


    # Add lines to InterProScan config file to activate SignalP 
    echo "Adding SignalP to InterProScan"

    echo -en "\n\n" >> software/interproscan/interproscan.properties
    echo "# SignalP" >> software/interproscan/interproscan.properties
    echo "signalp_euk.signature.library.release=4.1" >> software/interproscan/interproscan.properties
    echo "signalp_gram_positive.signature.library.release=4.1" >> software/interproscan/interproscan.properties
    echo "signalp_gram_negative.signature.library.release=4.1" >> software/interproscan/interproscan.properties
    echo "binary.signalp.path=bin/signalp/4.1/signalp" >> software/interproscan/interproscan.properties
    echo "signalp.perl.library.dir=bin/signalp/4.1/lib" >> software/interproscan/interproscan.properties

    echo "Finished installing SignalP"

    Install TMHMM database
    Print SignalP manual installation message
    echo
    echo "# Installing TMHMM #"
    echo "Follow the link bellow to the download page for TMHMM, read and accept the license."
    echo "You will receive an email with the download link for the software when you accept."
    echo 
    echo "$TMHMM_DOWNLOAD_PAGE"
    echo 

    ANSWER="No"  # set the variable for input loop

    while [[ $ANSWER =~ ^[Nn][Oo]$ ]]
    do
        read -p "Paste the link from the email here: " LINK
        read -p "Check the link is correct. Do you confirm the link? [yes/no]: " ANSWER
    done

    # Download and decompress the program
    echo "Downloading file"
    # wget -P software/temp $LINK
    wget -P software/temp -nd  -np -r -nH -l1 -A "tmhmm*.tar.gz" $LINK

    tar -xvzf software/temp/tmhmm-*.tar.gz -C software/temp/

    # Move files to InterProScan folder
    echo "Moving required files to InterProScan file tree"
    mv software/temp/tmhmm*/bin/* software/interproscan/bin/tmhmm/2.0c/
    mv software/temp/tmhmm*/lib/* software/interproscan/data/tmhmm/2.0c/

    # Add lines to InterProScan config file to activate TMHMM 
    echo "Adding TMHMM to InterProScan"

    echo -en "\n\n" >> software/interproscan/interproscan.properties
    echo "# TMHMM" >> software/interproscan/interproscan.properties
    echo "tmhmm.signature.library.release=2.0c" >> software/interproscan/interproscan.properties
    echo "binary.tmhmm.path=bin/tmhmm/2.0c/decodeanhmm.Linux_x86_64" >> software/interproscan/interproscan.properties
    echo "tmhmm.model.path=data/tmhmm/2.0c/TMHMM2.0.model" >> software/interproscan/interproscan.properties

    echo "Finished installing TMHMM"

    echo "# Removing leftover files #"
    rm -rf software/temp

    echo -en "\n\n"
    echo "Your installation of InterProScan is complete!"

elif [ $execution_option == "b" ]; then
    # Modify the `signalp` binary file so it can execute correctly
    echo "Configuring SignalP"
    PATH_TO_SIGNALP=$(realpath software/interproscan/bin/signalp/4.1)
    mv software/interproscan/bin/signalp/4.1/signalp software/interproscan/bin/signalp/4.1/signalp_backup
    awk -v var="    \$ENV{SIGNALP} = '$PATH_TO_SIGNALP'" '{ if (NR == 13) print var; else print $0}' software/interproscan/bin/signalp/4.1/signalp_backup > software/interproscan/bin/signalp/4.1/signalp

else
    echo "Try again and choose option a or b."
    exit 1

fi 